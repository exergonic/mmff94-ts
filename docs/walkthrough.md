# mmff94-ts — End-to-End Walkthrough

**What happens when you call `calc_energy(typed_mol)` or `optimize_lbfgs(typed_mol)`?**

This document traces the full pipeline: from a raw SDF string to a per-atom energy
decomposition or a minimized geometry. It is the companion to the AGENTS.md project
guide — this one explains *how the pieces fit together* rather than *why the project
exists*.

---

## 1. Entry point — the public API

```typescript
import { parse_sdf, assign_atom_types, calc_energy } from 'mmff94-ts';

const mol = parse_sdf(sdf_text);          // Step 1
const typed = assign_atom_types(mol);     // Step 2
const energy = calc_energy(typed);        // Step 3
console.log(energy.total);                // kcal/mol
```

Every function is pure — no hidden state, no `new ForceField()` — just data flowing
through a pipeline. The output of each step is the input of the next.

---

## 2. The data model — `src/types.ts`

Everything starts and ends with plain objects. No classes, no prototypes, no surprises.

```
Molecule                      TypedMolecule
┌──────────────┐              ┌──────────────────┐
│ atoms[]      │  assign_atom │ atoms[]          │
│   index      │  _types()    │ bonds[]          │
│   element    │ ───────────→ │ atom_types[]     │  ← MMFF94 type per atom
│   x, y, z    │              │ partial_charges? │  ← from BCI (optional)
│ bonds[]      │              └──────────────────┘
│   atom1      │                        │
│   atom2      │               calc_energy()
│   bond_order │                        ▼
│   name?      │              EnergyComponents
└──────────────┘              ┌──────────────────┐
                              │ total            │
                              │ bond_stretch     │
                              │ angle_bend       │
                              │ stretch_bend     │
                              │ torsion          │
                              │ van_der_waals    │
                              │ electrostatic    │
                              │ out_of_plane     │
                              └──────────────────┘
```

A `TypedMolecule` extends `Molecule` with one additional array — `atom_types`, which
maps each atom index to an MMFF94 integer type number (1–99). This type number is
the key that unlocks the right parameters for every energy term. Getting it wrong
means every subsequent calculation is wrong — which is why atom typing is the single
hardest piece.

Partial charges (`partial_charges[]`) are attached by `assign_bci_charges()`,
which returns a copy of the molecule with the charges on it — pure like the rest
of the pipeline (the input is untouched). They are reused by the gradient
calculation and stay valid through geometry optimization because MMFF94 charges
are geometry-independent. The energy terms also compute the charges on demand if
given a bare typed molecule, so you can skip the BCI step entirely if you only
need vdW + bonds + angles.

The top-level functions take this further: `calc_energy()`, `calc_gradient()`,
and the optimizers accept a bare `Molecule` straight from `parse_sdf()` and run
typing + charging on demand (`prepare.ts`). The simple path returns exactly what
the explicit path returns — the same full per-term breakdown, and for the
optimizers the same result with the typed/charged molecule at the minimum.

---

## 3. Input — `src/sdf.ts`

The `parse_sdf()` function reads a standard V2000 MOL block (the format that PubChem,
RDKit, and OpenBabel all emit). It does not use a tokenizer or a parser combinator —
it reads fixed-width columns, because that is what the V2000 spec requires.

```
Line 1:     molecule name (optional)
Line 2-3:   program/comment headers (ignored)
Line 4:     counts line — cols 0-2 = #atoms, cols 3-5 = #bonds
Lines 5+:   atom block — one line per atom
              cols 0-9:  x coordinate (Å)
              cols 10-19: y coordinate (Å)
              cols 20-29: z coordinate (Å)
              cols 31-33: element symbol (right-justified)
Next N:     bond block — one line per bond
              cols 0-2:  atom1 index (1-based)
              cols 3-5:  atom2 index (1-based)
              cols 6-8:  bond order (1=single, 2=double, 3=triple)
```

The function returns `{ atoms: [], bonds: [] }` on any parse error, so the caller
can fall through gracefully without a try/catch. It validates element symbols
against a known list and rejects obviously wrong structures (>999 atoms, negative
indices).

```typescript
// Internal helper: reads a fixed-width float from a string
function parse_float(line: string, start: number, end: number): number {
  const fragment = line.substring(start, end).trim();
  if (fragment === '') return NaN;
  return parseFloat(fragment);
}
```

---

## 4. Geometry primitives — `src/utils/vector.ts`

Every energy term needs at least one of these: distance, angle, or dihedral angle.
These are pure functions on `[x, y, z]` tuples.

```
Type:  Vec3 = [number, number, number]
```

| Function | What it computes | Used by |
|---|---|---|
| `distance(a, b)` | Euclidean distance \|a − b\| | bond-stretch, vdW, electrostatic |
| `angle_in_radians(b, a, c)` | Angle a-b-c at vertex b | angle-bend, stretch-bend |
| `dihedral_angle(i, j, k, l)` | Torsion angle i-j-k-l | torsion |

The dihedral angle uses the IUPAC sign convention (right-hand rule about j→k):

```
τ = atan2( (n₁ × n₂) · v̂₂ ,  n₁ · n₂ )
```

where n₁ and n₂ are the normals of planes (i,j,k) and (j,k,l), and v₂ is the
central bond vector j→k. This gives τ = 0° when i-j and k-l are eclipsed (cis),
and τ = ±180° when they are staggered (trans).

---

## 5. Atom typing — `src/mmff94/assign-atom-types.ts`

This is the hardest single piece. Every atom type is an integer (1–99) that selects
the correct row from every parameter table. A wrong type means wrong bond lengths,
wrong angles, wrong torsions, wrong vdW radii — everything.

The decision tree follows Halgren 1996 (J. Comput. Chem. 17, 520–552) and considers,
in order:

1. **Element** — C, N, O, H, S, P, F, Cl, Br, I, Si, etc.
2. **Coordination number** — how many immediate bonded neighbors the input
   declares. The SDF/MOL parser reads **explicit atoms only**: a structure
   written without explicit hydrogens (common from sketchers and some
   converters) types silently wrong — every bare carbon falls to generic
   type 1, no H types exist, and the energies are meaningless. Feed this
   library structures with all hydrogens present.
3. **Bond orders to neighbors** — single, double, triple, aromatic
4. **Neighbor elements** — e.g., C=O carbonyl vs. C=C alkene
5. **Neighbor types** — e.g., carbonyl C (type 3) vs. aromatic C (type 37)
6. **Ring membership** — is the atom in a small ring (3-4 membered)?
7. **Formal charge** — positively/negatively charged variants of neutral types

The implementation is a flat series of conditions, each with a comment explaining
the chemical reason for the branch:

```typescript
// Carbon with 3 neighbors and one double bond: alkene (sp²) or carbonyl (sp²)?
// If the double bond goes to oxygen → type 3 (carbonyl C).
// If the double bond goes to carbon  → type 2 (alkene sp² C).
```

The current implementation covers the common organic elements: carbon (sp³, sp²,
carbonyl, acetylenic, cyclopropyl, cyclobutyl, aromatic), hydrogen (bonded to C,
O, N, S), oxygen (divalent, carbonyl), nitrogen (amine, imine, amide, aromatic),
sulfur (thiol/sulfide, thiocarbonyl, sulfoxide, sulfone), halogens, silicon, and
phosphorus. Ring detection uses iterative leaf-stripping; aromaticity is detected
via the Kekulé pattern — every ring atom carries exactly one ring double bond,
except a 5-ring may have exactly one N/O/S with none (the lone-pair donor:
pyrrole's N, furan's O, thiophene's S); V2000 aromatic bonds (order 4) settle
the ring directly. Fused rings, chorded cage rings, and saturated rings
(triazines, γ-pyrones) are rejected by the pattern.

`assign_bci_charges()` lives in its own module, `src/mmff94/charges.ts` — it is
the electrostatics model, not atom typing. The Bond Charge Increment model: each
bond contributes a fixed increment to both of its atoms; the partial charge on an
atom is the sum of all its bond contributions, plus the formal-charge correction
of part V eq. (15). Charged types carry a primary formal charge q⁰ (e.g. +1 on
quaternary N, −0.5 on a carboxylate oxygen, +1/3 on a guanidinium N); a NEGATIVE
q⁰ shares half of itself with the bonded atoms — q_i = (1 − α_i·crd_i)·q⁰_i +
Σ_k α_k·q⁰_k + Σ w_ik, where α (fcadj, part V Table III) is the sharing factor and
the neighbor sum uses the NEIGHBOR's α. Type 32 is environment-dependent: −0.5
on a carboxylate oxygen, 0 on sulfone/nitro/nitrate oxygens. The suite's
reference per-atom charges (the `.mmd` pchg column) pin the model:
138/140 typing-exact molecules reproduce them to < 1e-3 (`charges-suite.test.ts`);
the two thiosulfinate anions are excluded (BatchMin's dative adjustment). The
fixture logs' per-atom charges are pinned in `charges.test.ts`.

---

## 6. Parameters — `src/mmff94/parameters/`

### 6.1 Source of truth

The numeric parameters come from **OpenBabel**'s `.par` files (which are themselves
machine-readable transcriptions of Halgren's published tables). The extraction
script (`scripts/extract-mmff94-par.py`) reads these `.par` files and writes TypeScript
literal objects. The script is a format converter — it implements zero MMFF94 logic.

### 6.2 Key format and the class column

The first column of every parameter table is the parameter **class**: class 0
is the general entry, and higher classes hold context-specific alternatives —
conjugated single bonds (bond class 1), angles with BT-flagged bonds (classes
1/2), 3-ring angles (3/5/6), 4-ring angles (4/7/8), and the torsion classes
(1/2/4/5). The class is NOT a priority ladder — it is a chemical question:
which class applies is decided from the molecule's bonding, not by falling
through entries. For example, a C=C double bond and a conjugated C-C single
bond both have type pair (2,2) but different force constants:

```
0   2    2     9.505     1.333   C94    <- class 0: C=C double bond
1   2    2     5.310     1.430   C94    <- class 1: conjugated single bond (BTij = 1)
```

The extraction preserves these as separate keys by including the class:

```typescript
BOND_PARAMS = {
  "0-2-2": { k_b: 9.505, r0: 1.333 },   // class 0: C=C double
  "1-2-2": { k_b: 5.310, r0: 1.430 },   // class 1: conjugated C-C
  ...
}
```

### 6.3 The class system — `parameter-classes.ts`

`lookup_param()` is only the class-0 fallback. The real resolution is the class
system in `src/mmff94/parameters/parameter-classes.ts`:

1. **BTij** (bond-type flag, part V p. 620) — 1 for a single bond between two
   sbmb-flagged types (conjugated dienes, C(=O)–C(ar)), 1 for a non-ring bond
   between two aromatic types (biphenyl), 0 for everything else including
   aromatic ring bonds (BatchMin flags ring bonds aromatic; Kekulé input files
   don't, so an in-ring BFS check reproduces it).
2. **ATijk** (angle) = BT(i,j) + BT(j,k) with ring overrides → classes 1/2,
   3/5/6 (3-ring), 4/7/8 (4-ring). **TTijkl** (torsion) → classes 1/2/4/5.
   **STijk** (stretch-bend) is a remap of ATijk (GetStrBndType).
3. **Step-down chain** (part I p. 513) — within the class: exact terminal types
   first, then the EqLvl3/4/5 equivalence levels of the terminals (from
   `mmffdef.par`, extracted into `atom-type-properties.ts`).
4. **Empirical rules** — when the chain misses: the part-V rules in
   `empirical.ts` — bonds (eqs. 18-19), angles (the θ₀ protocol + eq. 20),
   torsions (rules a–h), and the
   default-fsb table (`default-stretch-bend.ts`) for stretch-bend.

The terms import these helpers from the parameters barrel; each energy term
stays a pure evaluation of geometry × resolved parameters.

### 6.4 Parameter tables

| File | Content | Keyed by | Entry count |
|---|---|---|---|
| `bond.ts` | Bond stretch: k_b (mdyn/Å), r₀ (Å) | `"c-t1-t2"` | ~500 |
| `angle.ts` | Angle bend: k_a (mdyn·Å/rad²), θ₀ (°) | `"c-t1-t2-t3"` | ~2350 |
| `stretch-bend.ts` | Stretch-bend: k_sb_IJK, k_sb_KJI | `"c-t1-t2-t3"` | ~290 |
| `torsion.ts` | Torsion: up to 3 Fourier terms (V_n, γ_n, n) | `"c-t1-t2-t3-t4"` | ~940 |
| `van-der-waals.ts` | VDW per-atom: R*, α_i, N_i, G_i, DA flag | atom type number | ~95 |
| `bci.ts` | Bond charge increments | `"c-t1-t2"` + per-atom defaults | ~600 |
| `out-of-plane.ts` | OOP bending: k_oop | `"t1-t2-t3-t4"` | ~120 |
| `parameters/atom-types.ts` | Type definitions: symbol, element, valence | atom type number | ~95 |
| `atom-type-properties.ts` | Per-type flags: crd, val, pilp, mltb, arom, lin, sbmb + EqLvl3/4/5 | atom type number | ~95 |
| `default-stretch-bend.ts` | Element-row default k_sb values (mmffdfsb.par) | `"row-row-row"` | 30 |
| `parameter-classes.ts` | BTij/ATijk/TTijkl/STijk class selection + class-scoped resolution (hand-written) | — | — |
| `empirical.ts` | Part-V empirical rules: bonds (eqs. 18-19), angles (θ₀ + eq. 20), torsions (rules a–h) (hand-written) | — | — |

(The `c` in the key format is the class column; `lookup.ts` provides
`lookup_param()` for class-0 wildcard fallback.)

---

## 7. Energy terms — `src/mmff94/energy/`

Each energy term is a standalone function in its own file. Every function has the
same signature:

```typescript
function calc_term_energy(molecule: TypedMolecule): number
```

They are computed independently and summed by `total.ts`. The functions do NOT
apply 1-4 scaling — that is handled by the orchestrator so that each term is
testable in isolation.

### 7.1 Bond stretching — `bond-stretch.ts`

```
E_bond = 143.9325 · (k_b / 2) · (r − r₀)² · [1 + cs · (r − r₀) + 7/12 · cs² · (r − r₀)²]
```

where `cs = −2 Å⁻¹` is the cubic stretch constant (Halgren1996, eq. (2)).

For every bond in the molecule:
1. Resolve k_b and r₀ by the pair of atom types — class-aware: a BTij = 1
   bond (conjugated single bond) uses the class-1 entry when one exists
   (e.g. a diene's central C–C uses `1-2-2`, not the C=C double-bond
   `0-2-2`)
2. Compute the current bond length r from atomic coordinates
3. Compute the harmonic term (leading quadratic) and the anharmonic correction
   (cubic and quartic terms from the Morse expansion), then multiply them:
   ```
   harmonic     = 143.9325 · (k_b / 2) · (r − r₀)²
   anharmonic   = 1 + cs · (r − r₀) + 7/12 · cs² · (r − r₀)²
   E           += harmonic · anharmonic
   ```

The factor 143.9325 converts from mdyn/Å to kcal/mol/Å². The ½ is the
standard harmonic oscillator prefactor — the parameter table stores the
full force constant k_b, not k_b/2 (see note in source). The cubic expansion
through fourth order approximates a Morse potential with alpha = 2 Å⁻¹.

### 7.2 Angle bending — `angle-bend.ts`

Halgren1996, eq. (3) — for conventional bond angles:

```
E_angle = 0.043844 · (k_a / 2) · Δθ² · (1 + cb · Δθ)
```

where `Δθ = θ − θ₀`, and `cb = −0.007 deg⁻¹` is the cubic bend constant.

For linear centers — atom types with the `lin` flag in `atom-type-properties.ts`
(e.g. sp carbon) — eq. (3) is replaced by eq. (4):

```
E_angle = 143.9325 · k_a · (1 + cos θ)
```

which avoids the singularity of the cubic expansion at θ = 180°. Stretch-bend
is omitted for these angles (the same flag drives both terms).

For every angle i-j-k (where j is the central atom):

1. Resolve k_a, θ₀ and the linear flag via `angle_parameters()`: the angle
   class ATijk (BT sum, with 3/5/6 for 3-rings and 4/7/8 for 4-rings) selects
   the class-scoped entry — small rings have their own θ₀ ≈ 60°/90° — then the
   EqLvl3/4/5 step-down chain, then the part-II empirical rules (θ₀ from the
   central atom's coordination number, k_a from the Z/C element tables)
2. Compute the current angle θ
3. If `linear`: use the cosine form (eq. 4)
4. Otherwise: compute the harmonic term and the anharmonic correction:
   ```
   harmonic     = 0.043844 · (k_a / 2) · Δθ²
   anharmonic   = 1 + cb · Δθ
   E           += harmonic · anharmonic
   ```

The factor 0.043844 converts from mdyn·Å/rad² to kcal/mol/deg². The ½
prefactor follows the same convention as bond stretch — the parameter table
stores the full k_a. The cubic correction with cb = −0.007 breaks the
symmetry of the harmonic well (compression stiffer than extension).

### 7.3 Stretch-bend cross term — `stretch-bend.ts`

```
E_sb = 2.51210 · [k_sb_IJK · (r_IJ − r_IJ0) + k_sb_KJI · (r_KJ − r_KJ0)] · (θ − θ₀)
```

This is a CLASS II force field term — it couples bond stretching with angle bending.
Most simpler force fields (UFF, GAFF) omit this term, but MMFF94 includes it
because bond lengths and angles are physically coupled: when an H-C-H angle in
methane closes, the C-H bonds shorten slightly.

For each angle i-j-k:
1. Resolve the stretch-bend class: STijk is a REMAP of the angle class
   (1↔2 split by which side carries the BT flag, 2→3, 3→5, 4→4, 5→6/7,
   6→8, 7→9/10, 8→11), so ring/BT angles look up different keys than the
   angle term's
2. Look up k_sb_IJK and k_sb_KJI from `STRETCH_BEND_PARAMS` by the class and
   the type triplet; if the lookup misses, the default-fsb table
   (`default-stretch-bend.ts`, element-row keyed) supplies small F values —
   BatchMin still evaluates the angle (a stretched Si–C bond around an
   untyped angle can carry several kcal/mol)
3. Look up r₀ for bonds i-j and k-j via the class-aware `bond_parameters()`
   (each bond uses its OWN sorted type pair, not the angle's terminals)
4. Take θ₀ from the same `angle_parameters()` call the angle term uses, so
   ring and BT-flagged angles share the reference angle; skip the angle if
   the center is linear (eq. 4)
5. Compute current distances and angle
6. Accumulate: `E += 2.51210 * (k_IJK * dr_IJ + k_KJI * dr_KJ) * dθ`

Two separate k_sb constants are used because asymmetric environments (e.g., C-C-O
vs O-C-C) have different coupling strengths for the two sides of the angle. Many
entries have k_sb_IJK = k_sb_KJI (symmetric angles), but the parameter table
stores both independently.

### 7.4 Torsion (dihedral) — `torsion.ts`

```
E_tors = Σ (V_n / 2) · [1 + cos(n · τ − γ_n)]   for n = 1, 2, 3
```

A Fourier series. Each dihedral type has up to three terms with different
periodicities. The parameter table stores V₁, V₂, V₃ with fixed phase shifts:
γ₁ = 0°, γ₂ = 180°, γ₃ = 0° (by convention). The V_n values are barrier heights
in kcal/mol.

The function evaluates every dihedral i-j-k-l around every bond — including
double bonds: an alkene's C=C torsion is real (V₂ ≈ 12 kcal/mol, holding the
alkene planar; the class-0 entries for double bonds have V₁ = V₃ = 0). Only
torsions with no substituent on either central atom are skipped by construction.

Parameters come from the torsion class TTijkl (part IV p. 609): 1 = central
bond BT-flagged (conjugated diene — V₂ ≈ 1.8), 2 = terminal bond BT-flagged,
4 = all four atoms in a 4-ring, 5 = non-aromatic 5-ring with an sp3 carbon,
else 0. Within the class, the step-down chain runs in the order-canonical
direction only — the par file stores each entry in ONE direction (decided by
an order index), and consulting the other direction first would let wildcard
defaults like `*-1-1-*` steal exact reversed entries (an H-C-C-C dihedral must
resolve to `0-1-1-1-5`, not the generic `0-0-1-1-0`). If the chain misses
entirely, the part-V empirical rules apply (`empirical.ts`, rules
a–h: V₂ from the U parameters, V₃ from the V parameters, negative V₂ for O/S
central pairs, with some combinations skipping the torsion).

### 7.5 Van der Waals — `van-der-waals.ts`

```
E_vdw = ε_ij · [ (1.07·R* / (r + 0.07·R*))⁷ · (1.12·R*⁷ / (r⁷ + 0.12·R*⁷) − 2) ]
```

This is a **buffered 14-7** potential — NOT the Lennard-Jones 12-6 used by UFF and
GAFF. The "buffer" terms (0.07·R* and 0.12·R*⁷) eliminate the singularity at r = 0
that plagues the standard LJ potential, giving a finite repulsive wall. This makes
MMFF94 more numerically stable during optimization when atoms may approach very
closely.

Combination rules for mixed-atom pairs:
- **A donor atom present (DA flag = 1, e.g. H–N, H–O)**: arithmetic mean
  R*_ij = 0.5 × (R*_i + R*_j). If the OTHER atom is an acceptor (DA = 2,
  e.g. carbonyl O, amine N), the pair is a hydrogen bond: ε is halved and
  R* is scaled by 0.8 (ε always uses the unscaled R*).
- **No donor**: Waldman–Hagler combination,
  R*_ij = 0.5 × (R*_i + R*_j) × [1 + 0.2 × (1 − e^(−12·Δ²))] with
  Δ = (R*_i − R*_j)/(R*_i + R*_j).
- ε_ij = 181.16 × G_i × G_j × α_i × α_j / (√(α_i/N_i) + √(α_j/N_j)) (Slater–Kirkwood)

The donor flags come from the DA column of `mmffvdw.par` (extracted into the
vdW table). At r = R*, the expression simplifies to E = −ε — the well depth —
because both buffer fractions become exactly 1. This property is verified in
the test suite.

### 7.6 Electrostatic — `electrostatic.ts`

```
E_elec = 332.0716 · q_i · q_j / (D · (R_ij + S))
```

Coulomb's law with partial charges from the BCI model (`charges.ts`). The
factor 332.0716 converts from e²/Å to kcal/mol. Eq. (6) of part III adds
**S = 0.05 Å, the electrostatic buffering constant**, to every distance —
MMFF94 never evaluates the bare 1/r form.

- The dielectric D = 1.0 in vacuo (the MMFF94 default; D = r is the
  distance-dependent alternative).
- Only pairs separated by **three or more bonds** are evaluated: 1-2 and
  1-3 pairs are excluded (the same pair list as the van der Waals term —
  this is why ammonia's electrostatic energy is zero).
- 1-4 pairs (exactly three bonds apart) are scaled by **0.75** — the
  scaling lives inside this term because the term functions return
  totals, not pair lists.

Consumes the `partial_charges[]` attached by `assign_bci_charges()` (the term
also computes them on demand if they are absent).

Validated against the reference logs (per-atom charges AND energies) and
against BatchMin: 759/761 suite molecules match the electrostatic
component at |Δ| <= 1e-4 — the two exclusions (AN11A, DOZNIP) are the
delocalized-anion reference anomalies (see VALIDATION.md).

### 7.7 Out-of-plane bending — `out-of-plane.ts`

```
E_oop = 0.043844 · (k_oop / 2) · χ²
```

MMFF94 uses a dedicated oop term rather than the **improper torsion** approach
used by UFF and GAFF. These are NOT the same thing:

- **Improper torsion** (UFF/GAFF): treats the out-of-plane deformation as a
  dihedral rotation, which couples it with the true torsion term. This is
  physically incorrect — out-of-plane motion and torsional rotation are independent
  degrees of freedom.

- **MMFF94 oop term**: a pure out-of-plane bending that measures the deformation
  angle χ at the central atom, defined as the angle between the vector from the
  central atom to any one substituent and the plane defined by the other two
  substituents and the central atom.

Applied to every tri-coordinate center (planar or pyramidal) with exactly three
bonded neighbors. The sign of k_oop encodes real chemistry: zero for amine N
(pyramidalization comes from angle-bend reference values), negative for amide N
(MMFF94 gives pyramidal amide nitrogen deliberately). Validated against
BatchMin: 761/761 suite molecules at |Δ| <= 1e-4.

---

## 8. Total energy and 1-4 scaling — `total.ts`

```typescript
export function calc_energy(molecule: TypedMolecule): EnergyComponents {
  const bond_stretch  = calc_bond_stretch_energy(molecule);
  const angle_bend    = calc_angle_bend_energy(molecule);
  const stretch_bend  = calc_stretch_bend_energy(molecule);
  const torsion       = calc_torsion_energy(molecule);
  const van_der_waals = calc_vdw_energy(molecule);
  const electrostatic = calc_electrostatic_energy(molecule);
  const out_of_plane  = calc_oop_energy(molecule);

  // 1-4 scaling: electrostatic ×0.75, applied INSIDE the electrostatic
  // term (it is the only scaled term). vdW is NOT scaled at 1-4
  // (Halgren 1996, p. 496). The terms return totals, not pair lists,
  // so total.ts cannot rescale individual pairs.

  const total = bond_stretch + angle_bend + stretch_bend +
                torsion + van_der_waals + electrostatic + out_of_plane;

  return { total, bond_stretch, angle_bend, stretch_bend,
           torsion, van_der_waals, electrostatic, out_of_plane };
}
```

### Why 1-4 scaling?

Atoms separated by exactly three bonds (1-4 pairs) interact both through the
torsion term AND through direct vdW/electrostatic pairs. Without scaling, these
interactions would be double-counted. The MMFF94 specification mandates:

| Interaction | 1-4 scale factor | Reason |
|---|---|---|
| Van der Waals | **none** | Halgren 1996 p. 496: "1,4-vdW interactions are not differentially scaled in MMFF94" — the ×0.5 of MM2/GAFF is a different convention |
| Electrostatic | 0.75 | BCI charges are parameterized with this scaling built in |

The scaling is applied inside the electrostatic term — it is the only scaled
term, and the term functions return totals (not pair lists), so `total.ts`
cannot rescale individual pairs. Each term stays testable in isolation; the
electrostatic term's own tests pin the ×0.75 against the reference logs.

**Current status**: implemented. The total is a plain sum of the seven terms;
the ×0.75 lives inside the electrostatic term (see §7.6).

---

## 9. Gradients — `src/mmff94/gradient/`

The gradient is the derivative of the total energy with respect to every atomic
coordinate:

```
∇E = [dE/dx₁, dE/dy₁, dE/dz₁, dE/dx₂, ...]
```

Each energy term contributes its own analytical derivative. The gradient files
mirror the energy/ layout: one file per term, each exporting a function that
returns the gradient contribution for that term.

The analytical gradients are derived from the same functional forms as the energy
terms. Two rules keep them honest:

1. **Same computational path** — every derivative goes through the chain rule on
   the exact expression the energy term evaluates (the shared helpers in
   `gradient/derivatives.ts` mirror `utils/vector.ts`'s normalization order and
   handedness). A mathematically equal but differently-ordered derivative would
   still agree with finite differences, but floating-point rounding would show up
   as test noise.
2. **Shared resolution** — parameter lookup and pair enumeration are exported
   helpers (`stretch_bend_angle_terms`, `torsion_terms`, `vdw_pair_parameters`,
   `oop_force_constant`, `is_1_4_pair`) used by BOTH the energy term and its
   gradient, so the two cannot disagree about which parameters or which pairs
   apply. The 1-4 electrostatic ×0.75 is applied by the same `is_1_4_pair` in
   both.

**Current status**: complete. All seven term gradients plus the total are
cross-checked against central finite differences in `tests/gradient.test.ts`
(δ = 10⁻⁶ Å, relative error < 10⁻⁵; worst observed error 8×10⁻⁸). The gradient is
the negative of the force on each atom: F_i = −∇_i E.

---

## 10. Optimization — `src/optimize/`

Geometry optimization finds the nearest local minimum of the energy surface by
iteratively adjusting atomic coordinates. Two algorithms are planned:

### 10.1 L-BFGS (primary — implemented)

Limited-memory Broyden–Fletcher–Goldfarb–Shanno — a quasi-Newton method that
builds an approximate inverse Hessian from a limited history (m = 10 by default)
of position and gradient changes, needing only the energy and its gradient at
each step — exactly what `calc_energy`/`calc_gradient` provide. The implementation
follows Nocedal & Wright (2nd ed.): Algorithm 7.5 (two-loop recursion with the
initial-Hessian scaling γ = sᵀy/yᵀy), Algorithm 3.5 (strong-Wolfe line search,
c1 = 1e-4, c2 = 0.9) and Algorithm 3.6 (the zoom with cubic interpolation).

Two robustness details matter on MMFF94's stiff surfaces (see the source
comments for the full why):

- The first line-search trial is α₀ = 1/γ (capped at a 2 Å physical step):
  without the γ compensation, a tiny γ (~10⁻³–10⁻⁵ on stiff surfaces) lets the
  α = 1 trial satisfy the Wolfe conditions trivially and the optimizer accepts
  γ-sized steps — convergent in theory, glacial in practice.
- History pairs with steps below the noise floor (≤ 10⁻⁴ Å) are discarded, and a
  non-descent two-loop direction falls back to steepest descent with a history
  reset: wall-limited steps on vdW canyons would otherwise poison the direction
  with ill-scaled curvature information.

### 10.2 Steepest descent (fallback)

Simple gradient descent with Armijo line search: x ← x − α·∇E with α
halved until E decreases by at least c₁·α·‖∇E‖² (c₁ = 1e-4). The first
trial is capped in physical space (no atom moves more than 2 Å), the same
wall-guard as L-BFGS's first trial. Used when the L-BFGS line search fails or
as a robustness fallback.

Both optimizers stop when the maximum |gradient| component falls below
`gradient_tolerance` (default 0.05 kcal/mol/Å).

**Current status**: on the user-authored `*_non-optimized.sdf` fixtures
(ethane, butane, water — genuinely strained starting geometries), both
L-BFGS and steepest descent converge at the spec, max|g| < 0.05, with
the final energies inside the per-fixture windows and the final bond
lengths at the expected values. `tests/optimization.test.ts` covers
both.

### 10.3 Considered and deferred: DIIS (GDIIS)

Direct inversion of the iterative subspace (Császár & Pulay, 1984) was
considered as a third optimizer. Deferred as superfluous: it assembles
the same gradient history L-BFGS already exploits (its extrapolation is
effectively a particular quasi-Newton scheme — Farkas & Schlegel, PCCP
2002), it has no descent guarantee without bolting on the same line-search
and trial-step machinery we already built, and its failure profile overlaps
ours (the error vector is the gradient, so a descent-only valley defeats
it as it defeats steepest descent). Its famous habitats — SCF
convergence and constrained optimization — are out of scope here (fixed BCI
charges, no constraints). Revisit only if constrained optimization or
transition-state search ever enters scope. See Schlegel, *WIREs Comput.
Mol. Sci.* 2011, *1*, 790–809.

---

## 11. The parameter extraction pipeline

Though not part of the runtime API, the parameter extraction is a critical piece of
the build process.

```
OpenBabel .par files
  │
  ▼
scripts/extract-mmff94-par.py   ← Python, reads fixed-width columns
  │
  ▼
src/mmff94/parameters/*.ts       ← Auto-generated TypeScript files
  │
  ▼
imported by energy/ modules
  │
  ▼
TypeScript compilation
```

The extraction script handles 10 `.par` files, each with its own column
layout:

| File | Section in script | Key format |
|---|---|---|
| `mmffbond.par` | `generate_bond` | `"c-t1-t2"` |
| `mmffang.par` | `generate_angle` | `"c-t1-t2-t3"` |
| `mmffstbn.par` | `generate_stretch_bend` | `"c-t1-t2-t3"` |
| `mmfftor.par` | `generate_torsion` | `"c-t1-t2-t3-t4"` |
| `mmffvdw.par` | `generate_vdw` | atom type number |
| `mmffchg.par` + `mmffpbci.par` | `generate_bci` | `"c-t1-t2"` + per-atom defaults |
| `mmffoop.par` | `generate_oop` | `"t1-t2-t3-t4"` |
| `mmffdef.par` | `generate_atom_types` | atom type number |
| `mmffprop.par` + `mmffdef.par` (levels) | `generate_properties` | atom type number |
| `mmffdfsb.par` | `generate_default_stretch_bend` | `"row-row-row"` |

The first column in most `.par` files is the **class code** (0 = general
entry, higher = context-specific alternatives; see §6.2). Older extraction
attempts that dropped this column caused silent data loss — for example, the
C=C double bond parameter (class 0) was overwritten by the conjugated C-C
single bond parameter (class 1) because both produced the key `"2-2"`.
Including the class in the key as `"0-2-2"` / `"1-2-2"` preserves both
entries.

---

## 12. The full pipeline — data flow diagram

```
User input                      Internal state              Output
──────────                      ──────────────              ──────

SDF string (.mol/.sdf)
  │
  ▼
┌─────────────┐
│ parse_sdf() │  src/sdf.ts
└──────┬──────┘
       │
       ▼
  ┌─────────┐
  │ Molecule │  { atoms[], bonds[], name? }
  └────┬────┘
       │
       ▼
┌──────────────────┐
│assign_atom_types()│  src/mmff94/assign-atom-types.ts
└─────────┬────────┘
          │
          ▼
    ┌───────────┐
    │TypedMolecule│  { ..., atom_types[] }
    └─────┬─────┘
          │
          ├─────────────────────────────────┐
          │                                 │
          ▼                                 ▼
┌───────────────────┐           ┌──────────────────────┐
│assign_bci_charges│           │calc_energy()         │
│ (optional)        │           │  bond-stretch        │
└─────────┬─────────┘           │  angle-bend          │
          │                     │  stretch-bend        │
          ▼                     │  torsion             │
  charged molecule              │  van-der-waals       │
                                │  electrostatic       │
                                │  out-of-plane        │
                                │  ┌──────────────┐   │
                                │  │ 1-4 scaling  │   │
                                │  └──────────────┘   │
                                └──────────┬───────────┘
                                           │
                                           ▼
                                  ┌─────────────────┐
                                  │EnergyComponents │
                                  └─────────────────┘

For optimization:

  TypedMolecule
       │
       ▼
┌────────────────┐
│calc_gradient() │  →  gradient[dE/dx_i, dE/dy_i, dE/dz_i]
└───────┬────────┘
        │
        ▼
┌─────────────────┐     loop     ┌────────────────┐
│optimize_lbfgs() │  ──────────→ │ update coords  │
│ (L-BFGS)        │  ←───────── │ check gradient │
└───────┬─────────┘     until    └────────────────┘
        │             converged
        ▼
┌─────────────────────┐
│OptimizationResult   │
│  molecule (coords)  │
│  energy (final)     │
│  converged?         │
└─────────────────────┘
```

---

## 13. Module dependency graph

```
src/index.ts  (public barrel)
  │
  ├── src/types.ts              ← no deps
  │
  ├── src/sdf.ts                ← types.ts
  │
  ├── src/mmff94/index.ts       ← re-exports
  │     │
  │     ├── assign-atom-types.ts ← types.ts
  │     │
  │     ├── energy/
  │     │     ├── bond-stretch.ts    ← types, vector, parameters
  │     │     ├── angle-bend.ts      ← types, vector, parameters
  │     │     ├── stretch-bend.ts    ← types, vector, parameters
  │     │     ├── torsion.ts         ← types, vector, parameters
  │     │     ├── van-der-waals.ts   ← types, parameters
  │     │     ├── electrostatic.ts   ← types, parameters
  │     │     ├── out-of-plane.ts    ← types, parameters
  │     │     └── total.ts           ← all of the above
  │     │
  │     ├── gradient/
  │     │     ├── total.ts           ← types + the seven term gradients
  │     │     ├── derivatives.ts     ← shared geometry-derivative helpers
  │     │     └── <term>.ts (×7)     ← types, parameters, derivatives
  │     │
  │     └── parameters/
  │           ├── index.ts       ← re-exports
  │           ├── lookup.ts      ← no deps (generic helper)
  │           ├── bond.ts        ← auto-generated
  │           ├── angle.ts       ← auto-generated
  │           ├── stretch-bend.ts ← auto-generated
  │           ├── default-stretch-bend.ts ← auto-generated
  │           ├── torsion.ts     ← auto-generated
  │           ├── van-der-waals.ts ← auto-generated
  │           ├── bci.ts         ← auto-generated
  │           ├── out-of-plane.ts ← auto-generated
  │           ├── atom-types.ts  ← auto-generated
  │           ├── atom-type-properties.ts ← auto-generated
  │           ├── parameter-classes.ts ← class selection + resolution
  │           └── empirical.ts  ← part-V rules (bond/angle/torsion)
  │
  ├── src/utils/vector.ts        ← no deps
  │
  └── src/optimize/
        ├── l-bfgs.ts            ← the oracle callback is optional: the
        │                          built-in default is calc_energy +
        │                          calc_gradient on the working copy
        └── steepest-descent.ts  ← same convention as L-BFGS
```

---

## 14. Testing strategy

Every energy term is tested **in isolation** before it is tested in combination.

| Test | What it checks | How |
|---|---|---|
| **Unit** | Single energy function returns the right value for a known geometry | Compute by hand or with a reference for a 2-3 atom test case (H₂ for bond stretch, H₂O for angle bend) |
| **Regression** | Total and per-component energies match Halgren suite | OPTIMOL totals from `MMFF94.energies`, component breakdowns from `MMFF94_bmin.log`; assert \|computed − reference\| < 0.01 kcal/mol |
| **Gradient** | Analytical dE/dx matches (E(x+δ) − E(x−δ)) / (2δ) | Finite-difference on every coordinate of every atom in every fixture: δ = 10⁻⁶ Å, relative error < 10⁻⁵ |
| **Optimization** | After minimization, max gradient < threshold and energy is lower | L-BFGS and steepest descent on each `*_non-optimized.sdf` fixture (DONE — 3/3 at max\|g\| < 0.05) |

---

## 15. Current implementation status

| Component | Status | Tests |
|---|---|---|
| Data types | ✅ Complete | — |
| SDF parser | ✅ Complete | 5 tests |
| Vector math | ✅ Complete | 13 tests |
| Atom typing | ✅ Implemented | suite scoreboard: 761/761 type-exact vs OpenBabel (atom-types-suite.test.ts) |
| BCI charges | ✅ Implemented (charges.ts) | 6 tests (reference-log pins) |
| Bond stretch | ✅ Implemented | 2 tests |
| Angle bend | ✅ Implemented | 2 tests |
| Stretch-bend | ✅ Implemented | 3 tests |
| Torsion | ✅ Implemented | 4 tests |
| Van der Waals | ✅ Implemented | 5 tests |
| Electrostatic | ✅ Implemented (buffered r+0.05, 1-4 ×0.75) | reference logs + suite 759/761 ≤1e-4 (two documented exclusions) |
| Out-of-plane | ✅ Implemented | 12 tests |
| 1-4 scaling | ✅ Applied inside the electrostatic term | — |
| Total energy | ✅ Sums all seven terms | 8 tests (reference + suite comparison) |
| Gradients | ✅ Analytical (all 7 terms, shared helpers with the energy terms) | 9 tests (FD-verified, worst error 8×10⁻⁸) |
| L-BFGS | ✅ Implemented (Nocedal & Wright Alg. 7.5 + strong-Wolfe) | 16 tests (3/3 fixtures at max\|g\| < 0.05) |
| Steepest descent | ✅ Implemented (Armijo line search, wall-guarded) | same suite (3/3 fixtures at max\|g\| < 0.05) |
| MMD parser (Halgren suite) | ✅ Complete | 4 tests |
| **All tests** | **213 passing (4 skipped)** | **23 files** |
