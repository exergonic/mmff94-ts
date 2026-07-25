# mmff94-ts — End-to-End Walkthrough

**What happens when you call `calc_energy(typed_mol)` or `optimize(typed_mol)`?**

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

Partial charges (`partial_charges[]`) are filled in by `compute_bci_charges()`, which
must be called before energy terms that need electrostatics. They are stored on the
TypedMolecule rather than computed on-the-fly because the same charges are reused
by the gradient calculation. But if you only need vdW + bonds + angles, you can
skip the BCI step.

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
| `distance(a, b)` | Euclidean distance `|a - b|` | bond-stretch, vdW, electrostatic |
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

## 5. Atom typing — `src/mmff94/atom-types.ts`

This is the hardest single piece. Every atom type is an integer (1–99) that selects
the correct row from every parameter table. A wrong type means wrong bond lengths,
wrong angles, wrong torsions, wrong vdW radii — everything.

The decision tree follows Halgren 1996 (J. Comput. Chem. 17, 520–552) and considers,
in order:

1. **Element** — C, N, O, H, S, P, F, Cl, Br, I, Si, etc.
2. **Coordination number** — how many immediate neighbors (including implicit H)
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

The current implementation is a stub that only distinguishes by element (type 1 for
C, type 5 for H, type 8 for N, type 6 for O, etc.) — the full decision tree is a
work-in-progress.

Also in this file: `compute_bci_charges()`, which computes partial charges using
the Bond Charge Increment model. Each bond contributes a fixed increment to both of
its atoms; the partial charge on an atom is the sum of all its bond contributions.
This is currently a stub.

---

## 6. Parameters — `src/mmff94/parameters/`

### 6.1 Source of truth

The numeric parameters come from **OpenBabel**'s `.par` files (which are themselves
machine-readable transcriptions of Halgren's published tables). The extraction
script (`scripts/extract-params.ps1`) reads these `.par` files and writes TypeScript
literal objects. The script is a format converter — it implements zero MMFF94 logic.

### 6.2 Key format and the priority column

OpenBabel's parameter files may have multiple entries for the same type combination,
distinguished by a leading priority code (0 = default/primary, 1 = alternative,
etc.). For example, a C=C double bond and a conjugated C-C single bond both have
type pair (2,2) but different force constants:

```
0   2    2     9.505     1.333   C94    <- C=C double bond (priority 0)
1   2    2     5.310     1.430   #C94   <- conjugated single bond (priority 1)
```

The extraction preserves these as separate keys by including the priority:

```typescript
BOND_PARAMS = {
  "0-2-2": { k_b: 9.505, r0: 1.333 },   // priority 0: C=C double
  "1-2-2": { k_b: 5.310, r0: 1.430 },   // priority 1: conjugated C-C
  ...
}
```

### 6.3 The lookup system — `lookup.ts`

The `lookup_param()` function implements ordered parameter lookup:

1. **Priority 0** — try the default entry first
2. **Priority 1, 2, ...** — fall through to alternatives
3. **Wildcard at terminals** — replace type_i or type_k with 0
4. **Wildcard at both terminals** — replace both with 0

```typescript
// For a C-C bond (types [1, 1]), tries keys in this order:
//   "0-1-1" → found → return
//   "1-1-1" → skip (not found)
//   ...
//   "0-0-1" → wildcard at terminal → skip
//   "0-1-0" → wildcard at other terminal → skip
```

This mirrors OpenBabel's own lookup strategy: exact match first, then wildcards
from most-specific to least-specific.

### 6.4 Parameter tables

| File | Content | Keyed by | Entry count |
|---|---|---|---|
| `bond.ts` | Bond stretch: k_b (mdyn/Å), r₀ (Å) | `"p-t1-t2"` | ~500 |
| `angle.ts` | Angle bend: k_a (mdyn·Å/rad²), θ₀ (°) | `"p-t1-t2-t3"` | ~2350 |
| `stretch-bend.ts` | Stretch-bend: k_sb_IJK, k_sb_KJI | `"p-t1-t2-t3"` | ~290 |
| `torsion.ts` | Torsion: up to 3 Fourier terms (V_n, γ_n, n) | `"p-t1-t2-t3-t4"` | ~940 |
| `van-der-waals.ts` | VDW per-atom: R*, α_i, N_i, G_i | atom type number | ~95 |
| `bci.ts` | Bond charge increments | `"p-t1-t2"` + per-atom defaults | ~600 |
| `out-of-plane.ts` | OOP bending: k_oop | `"t1-t2-t3-t4"` | ~120 |
| `atom-types.ts` | Type definitions: symbol, element, valence | atom type number | ~95 |

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
1. Look up k_b and r₀ from `BOND_PARAMS` by the pair of atom types
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

```
E_angle = 0.043844 · k_a · (θ − θ₀)²
```

For every angle i-j-k (where j is the central atom):

The function identifies angles by walking the adjacency list: for each atom j
with N neighbors, there are N×(N−1)/2 angle pairs (i, k).

1. Look up k_a and θ₀ from `ANGLE_PARAMS` by the triplet of types (i, j, k)
2. Compute the current angle θ using `angle_in_radians()`
3. Accumulate: `E += 0.043844 * k_a * (θ_deg − θ₀)²`

The factor 0.043844 converts from mdyn·Å/rad² to kcal/mol/deg². The angle is
converted to degrees for the computation because θ₀ is stored in degrees.

### 7.3 Stretch-bend cross term — `stretch-bend.ts`

```
E_sb = 2.51210 · [k_sb_IJK · (r_IJ − r_IJ0) + k_sb_KJI · (r_KJ − r_KJ0)] · (θ − θ₀)
```

This is a CLASS II force field term — it couples bond stretching with angle bending.
Most simpler force fields (UFF, GAFF) omit this term, but MMFF94 includes it
because bond lengths and angles are physically coupled: when an H-C-H angle in
methane closes, the C-H bonds shorten slightly.

For each angle i-j-k:
1. Look up k_sb_IJK and k_sb_KJI from `STRETCH_BEND_PARAMS`
2. Look up r₀ for bonds i-j and k-j from `BOND_PARAMS`
3. Look up θ₀ from `ANGLE_PARAMS`
4. Compute current distances and angle
5. Accumulate: `E += 2.51210 * (k_IJK * dr_IJ + k_KJI * dr_KJ) * dθ`

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

The function evaluates every dihedral i-j-k-l where the central bond j-k is a
single bond. If the forward type order (i, j, k, l) doesn't match, the reverse
(l, k, j, i) is tried as fallback. Double and triple bonds have no torsion
potential — their planarity is enforced by the angle bend and out-of-plane terms.

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
- R*_ij = 0.5 × (R*_i + R*_j) (arithmetic mean)
- ε_ij = 181.16 × G_i × G_j × α_i × α_j / [α_i / √N_i + α_j / √N_j] (Slater-Kirkwood)

At r = R*, the expression simplifies to E = −ε — the well depth — because both
buffer fractions become exactly 1. This property is verified in the test suite.

### 7.6 Electrostatic — `electrostatic.ts` (stub)

```
E_elec = 332.0716 · q_i · q_j / (ε · r)
```

Coulomb's law with partial charges from the BCI model. The dielectric constant ε
depends on the environment: ε = 1 for in-vacuo, ε = r for distance-dependent
dielectric (the MMFF94 default), ε = 4.0 for protein/interior calculations.

Requires `compute_bci_charges()` to have been called first to fill in the
`partial_charges[]` array.

Current status: **stub** (returns 0).

### 7.7 Out-of-plane bending — `out-of-plane.ts` (stub)

```
E_oop = 0.043844 · k_oop · χ²
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

Applied to every sp² center (carbonyl C, olefinic C, aromatic C, trigonal planar N)
with exactly three bonded neighbors.

Current status: **stub** (returns 0).

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

  // Apply 1-4 scaling to vdW and electrostatic
  //   vdW:           multiply by 0.5
  //   electrostatic: multiply by 0.75
  // (not yet implemented)

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
| Van der Waals | 0.5 | Torsion barriers already capture some steric repulsion |
| Electrostatic | 0.75 | BCI charges are parameterized with this scaling built in |

The scaling is applied in `total.ts` rather than in the individual term functions
so that each term is simple and testable in isolation. The individual functions
compute the FULL (unscaled) energy; the orchestrator decides which pairs to scale.

**Current status**: 1-4 scaling is not yet implemented. The total is a plain sum.

---

## 9. Gradients — `src/mmff94/gradient/` (stub)

The gradient is the derivative of the total energy with respect to every atomic
coordinate:

```
∇E = [dE/dx₁, dE/dy₁, dE/dz₁, dE/dx₂, ...]
```

Each energy term contributes its own analytical derivative. The gradient files
mirror the energy/ layout: one file per term, each exporting a function that
returns the gradient contribution for that term.

The analytical gradients are derived from the same functional forms as the energy
terms and are cross-checked against finite-difference calculations (δ = 10⁻⁶ Å,
relative error < 10⁻⁵).

The gradient is the negative of the force on each atom: F_i = −∇_i E.

**Current status**: stub (returns zero vector for every atom).

---

## 10. Optimization — `src/optimize/` (stubs)

Geometry optimization finds the nearest local minimum of the energy surface by
iteratively adjusting atomic coordinates. Two algorithms are planned:

### 10.1 L-BFGS (primary)

Limited-memory Broyden–Fletcher–Goldfarb–Shanno — a quasi-Newton method that
builds an approximate Hessian from the history of gradient evaluations. Uses
cubic line search (standard algorithm from Nocedal & Wright).

### 10.2 Steepest descent (fallback)

Simple gradient descent with Armijo line search. Used when the L-BFGS line search
fails or for the first few iterations to improve the starting point.

Both optimizers stop when the maximum force component falls below a threshold
(typically 0.05 kcal/mol/Å).

**Current status**: both are stubs.

---

## 11. The parameter extraction pipeline

Though not part of the runtime API, the parameter extraction is a critical piece of
the build process.

```
OpenBabel .par files
  │
  ▼
scripts/extract-params.ps1      ← PowerShell, reads fixed-width columns
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

The extraction script handles 8 separate `.par` files, each with its own column
layout:

| File | Section in script | Key format |
|---|---|---|
| `mmffbond.par` | `Generate-Bond` | `"p-t1-t2"` |
| `mmffang.par` | `Generate-Angle` | `"p-t1-t2-t3"` |
| `mmffstbn.par` | `Generate-StretchBend` | `"p-t1-t2-t3"` |
| `mmfftor.par` | `Generate-Torsion` | `"p-t1-t2-t3-t4"` |
| `mmffvdw.par` | `Generate-Vdw` | atom type number |
| `mmffchg.par` + `mmffpbci.par` | `Generate-Bci` | `"p-t1-t2"` + per-atom defaults |
| `mmffoop.par` | `Generate-Oop` | `"t1-t2-t3-t4"` |
| `mmffdef.par` + `mmffprop.par` | `Generate-AtomTypes` | atom type number |

The first column in most `.par` files is a **priority code** (0 = default, 1+
= alternative). Older extraction attempts that dropped this column caused silent
data loss — for example, the C=C double bond parameter (priority 0) was overwritten
by the conjugated C-C single bond parameter (priority 1) because both produced the
key `"2-2"`. Including the priority in the key as `"0-2-2"` / `"1-2-2"` preserves
both entries.

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
│assign_atom_types()│  src/mmff94/atom-types.ts
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
│compute_bci_charges│           │calc_energy()         │
│ (optional)        │           │  bond-stretch        │
└─────────┬─────────┘           │  angle-bend          │
          │                     │  stretch-bend        │
          ▼                     │  torsion             │
  partial_charges[]             │  van-der-waals       │
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
┌────────────────┐     loop     ┌────────────────┐
│optimize()      │  ──────────→ │ update coords  │
│ L-BFGS or SD   │  ←───────── │ check gradient │
└───────┬────────┘     until    └────────────────┘
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
  │     ├── atom-types.ts       ← types.ts
  │     │
  │     ├── energy/
  │     │     ├── bond-stretch.ts    ← types, vector, parameters
  │     │     ├── angle-bend.ts      ← types, vector, parameters
  │     │     ├── stretch-bend.ts    ← types, vector, parameters
  │     │     ├── torsion.ts         ← types, vector, parameters (stub)
  │     │     ├── van-der-waals.ts   ← types, parameters (stub)
  │     │     ├── electrostatic.ts   ← types, parameters (stub)
  │     │     ├── out-of-plane.ts    ← types, parameters (stub)
  │     │     └── total.ts           ← all of the above
  │     │
  │     ├── gradient/
  │     │     └── total.ts       ← types (stub)
  │     │
  │     └── parameters/
  │           ├── index.ts       ← re-exports
  │           ├── lookup.ts      ← no deps (generic helper)
  │           ├── bond.ts        ← auto-generated
  │           ├── angle.ts       ← auto-generated
  │           ├── stretch-bend.ts ← auto-generated
  │           ├── torsion.ts     ← auto-generated
  │           ├── van-der-waals.ts ← auto-generated
  │           ├── bci.ts         ← auto-generated
  │           ├── out-of-plane.ts ← auto-generated
  │           └── atom-types.ts  ← auto-generated
  │
  ├── src/utils/vector.ts        ← no deps
  │
  └── src/optimize/
        ├── l-bfgs.ts            ← types, vector, gradient (stub)
        └── steepest-descent.ts  ← types, vector, gradient (stub)
```

---

## 14. Testing strategy

Every energy term is tested **in isolation** before it is tested in combination.

| Test | What it checks | How |
|---|---|---|
| **Unit** | Single energy function returns the right value for a known geometry | Compute by hand or with a reference for a 2-3 atom test case (H₂ for bond stretch, H₂O for angle bend) |
| **Regression** | Total energy matches OpenBabel/RDKit for ~20 molecules | `.sdf` fixture files with sidecar `.json` reference energies; assert `|computed − reference| < 0.01 kcal/mol` |
| **Gradient** | Analytical dE/dx matches (E(x+δ) − E(x−δ)) / (2δ) | Finite-difference on every coordinate of every atom in every fixture: δ = 10⁻⁶ Å, relative error < 10⁻⁵ |
| **Optimization** | After minimization, max gradient < threshold and energy is lower | Run L-BFGS on each fixture from the SDF geometry and from a perturbed geometry |

---

## 15. Current implementation status

| Component | Status | Tests |
|---|---|---|
| Data types | ✅ Complete | — |
| SDF parser | ✅ Complete | 5 tests |
| Vector math | ✅ Complete | 12 tests |
| Atom typing | ✅ Implemented | 9 tests |
| BCI charges | ⚠️ Stub | 0 tests |
| Bond stretch | ✅ Implemented | 2 tests |
| Angle bend | ✅ Implemented | 1 test |
| Stretch-bend | ✅ Implemented | 2 tests |
| Torsion | ✅ Implemented | 3 tests |
| Van der Waals | ✅ Implemented | 4 tests |
| Electrostatic | ❌ Stub (returns 0) | 0 tests |
| Out-of-plane | ❌ Stub (returns 0) | 0 tests |
| 1-4 scaling | ❌ Not implemented | 0 tests |
| Total energy | ⚠️ Sums all stubs | 0 tests |
| Gradients | ❌ Stub | 0 tests |
| L-BFGS | ❌ Stub | 0 tests |
| Steepest descent | ❌ Stub | 0 tests |
| MMD parser (Halgren suite) | ✅ Complete | 4 tests |
| **All tests** | **44 passing** | **9 files** |
