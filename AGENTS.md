# mmff94-ts

**Merck Molecular Force Field (MMFF94) in pure TypeScript** — energy evaluation, analytical gradients, and geometry optimization, running in the browser or Node.js without native dependencies.

---

## Reader's compact

This code is written for **chemists who code**, not for software engineers. Every decision in the architecture, naming, and commenting conventions serves that reader.

- **Names use chemistry vocabulary.** Where a concept has a well-known chemistry name, that name is the function or variable name. `calc_bond_stretch_energy`, `partial_charges_from_bci`, `buffered_14_7_vdw`. Generic programming concepts (loop counters, array indices, plain math helpers) get plain names. A forced chemistry name on a generic thing is worse than a plain one.
- **Comments explain *why*, not *what***. A comment that restates the next line of TypeScript is noise. A comment that explains *why* this term takes this functional form, *why* this atom is typed differently from its cousin, or *why* the parameter lookup falls through in a particular order — that is the target. Decision points where the chemistry, not the code, drives a branch are the places to comment.
- **Simple beats clever.** If a block takes more than a few sentences to explain why it is built the way it is, simplify it. The slightly longer but straightforward version is always preferred over the compact-but-opaque version.
- **No architecture astronautics.** Every abstraction layer must justify its existence with at least two concrete users. Abstract base classes, registries, and plugin systems are not added on day one — they are extracted once a pattern has proven itself. v1 is a single force field (MMFF94) in a single module, with plain functions.

---

## Project notes

Project notes and errata live in `NOTES.md` — read it when relevant, update it freely. This file stays the stable constitution; keep volatile context in the notes file.

---

## Why this project exists

RDKit.js (WASM-based) and OpenBabel's JavaScript bindings did not meet our needs — large download size, complex build chains, error-prone fallback paths for the browser. The only other pure-JS MMFF94 we know of is openchemlib-js (a GWT transpile of the Actelion/DataWarrior force field): it implements all seven terms with the correct functional forms, but publishes no validation against the original MMFF94 program, exposes only a total-energy/minimiser API, and refuses molecules without explicit hydrogens. No pure-JS library offers what this one does — per-term numeric validation against the original program's own suite, cross-checked three ways against openchemlib in `tests/openchemlib-comparison.test.ts`.

`mmff94-ts` exists to fill that gap: a **small, fast, correct, browser-native** molecular mechanics library that the entire web-chemistry ecosystem can depend on. It happens to power Valence's fallback geometry, but Valence bends to `mmff94-ts` — not the other way around.

---

## Stack

| Tool | Purpose |
|---|---|
| **TypeScript** (strict) | Language — pure, no runtime dependencies |
| **tsc** | Compiler — ESM output in `dist/` |
| **Vitest** | Unit tests, regression tests, benchmarks |
| **Python** (CPython 3.10+) | Parameter extraction script (build-time only) |

Node.js lives at `C:/Users/mccan/scoop/apps/nodejs/current/node` (Scoop install).
Add it to `$env:Path` or use the full path when running scripts.

Zero npm dependencies at runtime. The library compiles to ESM in `dist/`
(`npm run build` = `tsc`).

---

## Project structure

```
mmff94-ts/
├── AGENTS.md                  # This file — project guide for agents and humans
├── package.json
├── tsconfig.json
├── vitest.config.ts
├── scripts/
│   └── extract-mmff94-par.py   # Reads OpenBabel .par files (from temp_ob/data/) → TS parameter tables
├── src/
│   ├── index.ts               # Public barrel: what a consumer imports from 'mmff94-ts'
│   ├── types.ts               # Molecule, Atom, Bond, TypedMolecule — the data model
│   ├── sdf.ts                 # SDF / MOL V2000 parser (standalone, not borrowed from Valence)
│   ├── mmff94/                # Everything MMFF94-specific
│   │   ├── index.ts           # Functions: assign_atom_types, assign_bci_charges, calc_energy, calc_gradient
│   │   ├── assign-atom-types.ts # MMFF94 atom type assignment — the decision tree (the hardest piece)
│   │   ├── charges.ts         # BCI partial charges (assign_bci_charges) — the electrostatics model
│   │   ├── prepare.ts         # prepare_molecule: type + charge on demand (the simple path)
│   │   ├── energy/
│   │   │   ├── bond-stretch.ts
│   │   │   ├── angle-bend.ts
│   │   │   ├── stretch-bend.ts
│   │   │   ├── torsion.ts
│   │   │   ├── van-der-waals.ts
│   │   │   ├── electrostatic.ts   # Coulombic term (eq. 6, part III: buffered r+0.05, 1-4 ×0.75)
│   │   │   ├── out-of-plane.ts
│   │   │   └── total.ts        # Sums all terms (electrostatic applies its own 1-4 ×0.75)
│   │   ├── gradient/
│   │   │   └── total.ts        # Sums the seven per-term gradients (Phase 5 complete)
│   │   └── parameters/
│       │   ├── index.ts            # Aggregates and exports all parameter tables + the class helpers + lookup_param
│       │   ├── atom-types.ts       # MMFF94 atom type definitions
│       │   ├── atom-type-properties.ts  # Per-type flags (crd/val/pilp/mltb/arom/lin/sbmb + EqLvl3/4/5)
│       │   ├── bond.ts             # Bond stretch parameters (k_b, r₀)
│       │   ├── angle.ts            # Angle bend parameters (k_a, θ₀)
│       │   ├── stretch-bend.ts     # Stretch-bend cross term parameters (k_sb)
│       │   ├── default-stretch-bend.ts  # Element-row default k_sb values (mmffdfsb.par)
│       │   ├── torsion.ts          # Torsion parameters (V_n, γ_n)
│       │   ├── van-der-waals.ts    # Van der Waals parameters (R*, ε, G_i, α_i, N_i)
│       │   ├── bci.ts             # Bond charge increments for electrostatics
│       │   ├── out-of-plane.ts    # Out-of-plane bending parameters (k_oop)
│       │   ├── parameters/parameter-classes.ts    # BTij/ATijk/TTijkl/STijk class selection + class-scoped resolution
│       │   ├── parameters/empirical.ts            # Part V on-the-fly generation (eqs. 18-20 + the torsion rules, Tables V/VI/X: bond k_b/r₀, angle θ₀/k_a, torsion V_n — the single home for the empirical machinery)
│       │   └── lookup.ts                # lookup_param() wildcard fallbacks (class-0 only)
│   ├── optimize/
│   │   ├── l-bfgs.ts              # L-BFGS minimizer (Nocedal & Wright, Alg. 7.5 + strong-Wolfe line search)
│   │   └── steepest-descent.ts    # Steepest descent fallback (Armijo line search)
│   └── utils/
│       ├── vector.ts               # 3D vector math (add, sub, dot, cross, normalize, distance, angles, rotate)
│       └── mmd-parser.ts           # BatchMin-format .mmd parser (validation suite structures)
├── tests/
│   ├── fixtures/                   # SDF files + known MMFF94 energies from OpenBabel/RDKit
│   │   ├── sdf/                    # Molecule SDFs for unit and regression tests
│   │   │   ├── ethane.sdf
│   │   │   ├── methane.sdf
│   │   │   └── ... (16 molecules)
│   │   └── validation-suite/      # Halgren's 761-molecule MMFF94 Validation Suite (Nov 1998 revision)
│   │       ├── README
│   │       ├── ANNOUNCE
│   │       ├── index.html
│   │       ├── MMFF94.energies    # Reference total energies (OPTIMOL + BatchMin)
│   │       ├── MMFF94.titles      # Molecule names
│   │       ├── MMFF94_bmin.log    # BatchMin per-component energy breakdowns
│   │       ├── MMFF94.dative_molecules
│   │       ├── MMFF94.mmd         # BatchMin-format structure files
│   │       ├── MMFF94_dative.mol2
│   │       ├── MMFF94_hypervalent.mol2
│   │       └── mmff94-atom-types.json  # Canonical atom types from OpenBabel (761 molecules; see tests/scripts/extract_suite_types.py)
│   ├── scripts/                    # Tools for generating reference energies
│   │   ├── get_mmff94_breakdown.py # Shells out to obabel.exe for reference energies
│   │   ├── extract_suite_types.py  # OpenBabel MMFF94 atom types for the suite → mmff94-atom-types.json
│   │   ├── obenergy.sh
│   │   ├── main.py
│   │   └── pyproject.toml          # uv project (openbabel-wheel)
│   ├── references/                 # obenergy output logs for regression comparison
│   ├── sdf.test.ts
│   ├── vector.test.ts
│   ├── mmd-parser.test.ts
│   ├── atom-types.test.ts          # Hardcoded OpenBabel types for the SDF fixtures
│   ├── atom-types-suite.test.ts    # Suite-wide typing coverage vs mmff94-atom-types.json + KNOWN_GOOD regressions
│   ├── bond-stretch.test.ts
│   ├── angle-bend.test.ts
│   ├── stretch-bend.test.ts
│   ├── torsion.test.ts
│   ├── van-der-waals.test.ts
│   ├── out-of-plane.test.ts
│   ├── reference-comparison.test.ts  # Fixture energies vs obenergy logs
│   └── validate-against-suite.test.ts # Per-component energies vs BatchMin log
└── examples/
    └── quickstart.ts               # Minimal usage example
```

---

## Architecture

### The data model (src/types.ts)

```typescript
interface Atom {
  index: number;
  element: string;         // 'C', 'N', 'O', 'H', etc.
  x: number;
  y: number;
  z: number;
}

interface Bond {
  atom1: number;           // index into atoms[]
  atom2: number;
  bond_order: number;      // 1 = single, 2 = double, 3 = triple
}

interface Molecule {
  atoms: Atom[];
  bonds: Bond[];
  name?: string;
}
```

A `TypedMolecule` extends `Molecule` with a per-atom MMFF94 atom type (an integer index into the atom type parameter table). This is the output of atom typing and the input to every energy term.

### The flow — one step at a time

```
Molecule (from SDF or sketcher)
    │
    ▼
assign_atom_types(molecule)   →   TypedMolecule
    │  Decides the MMFF94 type of every atom by walking the
    │  molecular graph. This is the hardest piece: the decision
    │  tree depends on element, neighbor elements, bond orders,
    │  ring membership, and formal charge.
    │
    ▼
assign_bci_charges(typed_mol) →  TypedMolecule (partial_charges[] attached)
    │  Bond charge increment model: each bond contributes a fixed
    │  charge increment to both atoms. Sum per atom = partial charge.
    │  Pure like the rest of the pipeline: returns a copy of the
    │  molecule with the charges attached, so the value flows on to
    │  the energy terms (which also compute the charges on demand if
    │  given a bare typed molecule).
    │
    │  SIMPLE PATH: the top-level functions (calc_energy,
    │  calc_gradient, the optimizers) accept a bare Molecule and run
    │  typing + charging on demand (prepare.ts) — one call, no
    │  preparation. The results are identical to the explicit path.
    │
    ▼
calc_energy(typed_mol)  →  EnergyComponents
    │  Seven energy terms, summed. Each term is a standalone
    │  function that reads molecule geometry + parameter tables
    │  and returns a kcal/mol value.
    │  The electrostatic term applies the 1-4 ×0.75 scaling itself
    │  (it is the only scaled term; terms return totals, not pairs).
    │
    ▼
calc_gradient(typed_mol)  →  dE/dx_i, dE/dy_i, dE/dz_i  (per atom)
    │  Analytical derivative of every term. Cross-checked against
    │  finite differences in the test suite (δ = 1e-6 Å, relative
    │  error < 1e-5).
    │
    ▼
optimize_lbfgs(typed_mol, energy_gradient_fn)  →  OptimizationResult
    │  Iteratively adjusts coordinates to minimize total energy.
    │  L-BFGS is the primary algorithm (Nocedal & Wright Alg. 7.5:
    │  two-loop recursion with γ = sᵀy/yᵀy initial-Hessian scaling,
    │  strong-Wolfe line search with cubic-interpolation zoom); steepest
    │  descent (Armijo backtracking) is the fallback. Stops when max
    │  gradient < threshold.
```

Everything is a pure function. There is no global state, no singleton, no hidden context.

### Energy terms — what each one does

Each term lives in its own file under `src/mmff94/energy/`. Every file exports a single function with a consistent signature:

```typescript
function calc_term_energy(molecule: TypedMolecule): number
```

| File | Chemistry term | Functional form | Why this form |
|---|---|---|---|
| `bond-stretch.ts` | Bond stretching | `E = 143.9325 · (k_b/2) · Δr² · [1 + cs·Δr + 7/12·cs²·Δr²]` | Halgren1996 eq. (2). The factor 143.9325 converts from mdyn/Å to kcal/mol/Å². |
| `angle-bend.ts` | Angle bending | `0.043844 · (k_a / 2) · Δθ² · (1 + cb · Δθ)` | Halgren1996 eq. (3). Factor 0.043844 converts mdyn·Å/rad² to kcal/mol/deg². Eq (4) for linear centers (the `lin` type flag — e.g. sp carbon); stretch-bend is omitted for those angles. |
| `stretch-bend.ts` | Stretch-bend cross term | `E = 2.51210 · [k_sb_IJK · (r₁−r₁₀) + k_sb_KJI · (r₂−r₂₀)] · (θ−θ₀)` | Halgren1996 eq. (5). Class-II term coupling bond and angle coordinates. Two k_sb values for asymmetric environments. Skipped for linear centers (eq. 4 angles) — the `lin` flag from angle_parameters. |
| `torsion.ts` | Torsion (dihedral) | `E = Σ (V_n / 2) · [1 + cos(n·τ − γ_n)]` for n=1,2,3 | Halgren1996 eq. (7). Fourier series: γ₁=0°, γ₂=180°, γ₃=0°. All bonds evaluated — an alkene's C=C torsion (V₂ ≈ 12) is real; double/triple bonds use the class-0 entries with V₁=V₃=0. |
| `van-der-waals.ts` | Van der Waals (non-bonded) | `E = ε_ij · [(1.07·R* / (r + 0.07·R*))⁷ · (1.12·R*⁷ / (r⁷ + 0.12·R*⁷) − 2)]` | Halgren1996 eq. (8). Buffered 14-7 potential. R* from eq. (9), Waldman-Hagler combination eqs. (10-11), Slater-Kirkwood ε eq. (12). |
| `electrostatic.ts` | Electrostatic (Coulombic) | `E = 332.0716 · q_i · q_j / (r + 0.05)` | Coulomb's law with partial charges from the BCI model (charges.ts). The factor 332.0716 converts from e²/Å to kcal/mol. Eq. (6) of part III adds the electrostatic buffering constant S = 0.05 Å to every distance. Pairs 1-2 and 1-3 are excluded (the closest pairs are 1-4); 1-4 pairs are scaled by 0.75. The dielectric D = 1.0 in vacuo (the MMFF94 default; D = r is the alternative solvent model). |
| `out-of-plane.ts` | Out-of-plane bending (tri-coordinate) | `E = 0.043844 · k_oop/2 · χ²` | Harmonic in the Wilson out-of-plane angle χ at **any tri-coordinate center** — planar (carbonyl C, olefinic C) or pyramidal (amine N, amide N). The sign of k_oop encodes the chemistry: zero for amine N (pyramidalization comes from angle-bend reference values), negative for amide N (see peculiarity #7). Functionally different from an "improper torsion". |

> **Note on OpenBabel factor discrepancy.** OpenBabel's `obenergy -ff MMFF94` uses the MM2 convention for its internal conversion constants (`bondunit 71.94`, `angleunit 0.02191418` — exactly half of the MMFF94 values listed above). These are defined in `data/mm2.prm`. The Halgren-paper factors (143.9325, 0.043844, 2.51210) are the authoritative MMFF94 values and are what this library implements. OpenBabel's reported energies are therefore roughly half of what this library computes for bond stretch and angle bend. Stretch-bend and torsion energies match directly since they use different conversion paths. Cross-check against Halgren's tabulated values or RDKit (which uses the full Halgren factors), not raw `obenergy` totals.

### 1-4 scaling

Atoms exactly three bonds apart (1-4 pairs) have their **electrostatic** interaction scaled by **0.75**. Their van der Waals interaction is **not** scaled — Halgren 1996 (p. 496): "1,4-vdW interactions are not differentially scaled in MMFF94" (the ×0.5 common in MM2/GAFF is a different convention). The scaling is applied inside the electrostatic term — it is the only scaled term, and the term functions return totals (not pair lists), so total.ts cannot rescale individual pairs. The electrostatic term also excludes 1-2 and 1-3 pairs entirely (the closest pairs it evaluates are 1-4), matching the van der Waals pair list and the suite's BatchMin energies (ammonia's electrostatic energy is zero).

---

## MMFF94 peculiarities — things that differ from UFF/GAFF

These are the features that make MMFF94 more accurate and also more complex. Getting them wrong means the library is not MMFF94.

1. **Buffered 14-7 vdW.** The functional form is distinct from Lennard-Jones. See `van-der-waals.ts` for the exact formula.
2. **Bond charge increments (BCI).** Partial charges are not stored atom-wise in the parameter tables. Instead, each bond type has a charge increment. The partial charge on an atom is the sum of the BCI values of every bond it participates in, plus the formal-charge correction of part V eq. (15): each charged type carries a primary formal charge q⁰ (e.g. +1 on quaternary N, −0.5 on a carboxylate oxygen, +1/3 on a guanidinium N) and atoms with a NEGATIVE q⁰ share half of it with their bonded neighbors — the neighbor sum uses the NEIGHBOR's α (fcadj): q_i = (1 − α_i·crd_i)·q⁰_i + Σ_k α_k·q⁰_k + Σ w_ik. Transcribing the sharing with the atom's own α (as some implementations do) breaks the carboxylate: the carbon (α = 0) must still receive half of each attached oxygen's −0.5. The q⁰ values were verified against the suite's per-atom reference charges (`.mmd` pchg column): all 138 charge-comparable typing-exact molecules reproduce them to < 1e-3 (the two thiosulfinate anions are dative-adjusted in the reference). Type 32 is environment-dependent: −0.5 on a carboxylate oxygen, 0 on sulfone/nitro/nitrate oxygens (whose polarization lives in the BCI).
3. **Stretch-bend cross term.** A class-II term that couples bond and angle coordinates. Most simplified force fields omit it.
4. **Torsion barriers are Fourier series, not single terms.** Each dihedral type has up to three V_n values (n=1,2,3). The periodicity and phase of each term must be respected — many ports incorrectly reduce this to a single-cosine.
5. **Parameter classes, not priorities.** The first column of every parameter table is the parameter CLASS: 0 is the general entry, and higher classes hold context-specific alternatives — conjugated single bonds (bond class 1), angles with BT-flagged bonds (1/2), 3-ring angles (3/5/6), 4-ring angles (4/7/8), torsion classes 1/2/4/5. Which class applies is a chemical question answered by `parameter-classes.ts`: the bond-type flag BTij (part V p. 620 — sbmb pairs, with aromatic ring bonds reading 0) feeds the angle class ATijk, the torsion class TTijkl, and the stretch-bend class STijk (a remap of ATijk). Within a class, the step-down chain (part I p. 513) tries exact types, then the EqLvl3/4/5 equivalence levels of the terminal types, then the empirical rules (the angle and torsion rules in `empirical.ts`). `lookup_param()`'s wildcards (leading → trailing → both terminal zeros) are the class-0 fallback only — do NOT use it to "find" class entries; the class must be selected explicitly. Callers whose wildcard shape differs (e.g. the out-of-plane term's per-central-type `"0-j-0-0"`) fall back explicitly after the helper misses.
6. **1-4 scaling.** Electrostatic × 0.75 only — Halgren 1996 p. 496: "1,4-vdW interactions are not differentially scaled in MMFF94". Do NOT apply a 0.5 vdW factor at 1-4; that is MM2/GAFF convention.
7. **Out-of-plane bending, not improper torsion.** MMFF94 uses a dedicated oop term — harmonic in the Wilson out-of-plane angle χ (eq. 6) — rather than the improper torsion formalism used by UFF/GAFF. It applies to **every tri-coordinate center, planar or not**, and the sign of k_oop encodes real chemistry:
   - **Zero for amine N (type 8).** The term is evaluated but contributes nothing; amine pyramidalization comes from angle-bend reference values that average below 120°. This is exactly Halgren's stated purpose (part I): "For trigonal nonplanar centers, this formulation allows angle-bending reference values that average less than 120° to be used to make the center pyramidal; the out-of-plane term can then be employed to improve the fit to the inversion barrier."
   - **Negative for amide N (type 10).** MMFF94 deliberately gives *pyramidal* amide nitrogen — an intentional departure from the planar-amide idealization, fit to MP2/6-31G* reference structures that show "nonplanarity at nitrogen ... for most amides." BatchMin's own validation-suite header states it flatly: MMFF94 "(gives pyramidal delocalized trigonal nitrogens)", and the suite's component energies include large negative OOP terms (e.g. −4.4 kcal/mol) from exactly this. The equilibrium pyramidalization is modest — the negative oop curvature is balanced against angle-bend and torsion.
   - **Do not "fix" the negative constants — they are the specification.** MMFF94s is the variant that restores planarity: it was created expressly "to yield nearly planar energy-minimized geometries for delocalized trigonal nitrogen," replacing the negative amide-N constants with small positive ones (+0.015 vs −0.02 for the type-10 wildcard).

---

## Atom typing — the hardest piece

Assigning the correct MMFF94 atom type is the gate to everything else. A wrong type cascades into wrong bond/angle/torsion/vdW parameters.

The typing rules are a **decision tree** based on:

1. **Element** — C, N, O, H, S, P, F, Cl, Br, I, etc.
2. **Coordination number** — number of bonded neighbors (including implicit H)
3. **Bond orders to neighbors** — single, double, triple, aromatic
4. **Neighbor elements** — what the bonded atoms are (e.g., C=O vs. C=C)
5. **Neighbor typing** — what type the bonded atoms are (e.g., carbonyl C vs. alkene C)
6. **Ring membership** — is the atom in a ring? Small ring (3-4 membered)?
7. **Formal charge** — positively or negatively charged variants of neutral types

The implementation in `src/mmff94/assign-atom-types.ts` should be a **flat series of conditions**, not a deep inheritance hierarchy. Comments at each branch explain the chemical reason:

```typescript
// Example comment style — not actual code
// Carbon with 3 neighbors and one double bond: could be alkene (sp²)
// or carbonyl (sp²). Check if the double bond goes to oxygen.
// If yes, and the other two neighbors are either carbon or hydrogen,
// assign type 3 (carbonyl C in aldehyde/ketone).
// If the double bond goes to carbon, assign type 2 (alkene sp² C).
```

The decision tree is derived from Halgren's original 1996 paper and cross-checked against OpenBabel's output. See the extraction script and parameter source notes below.

Typing coverage is measured objectively: `tests/atom-types-suite.test.ts` compares our types against OpenBabel's canonical types for the full 761-molecule validation suite (November 1998 revision; **761/761 type-exact — 100%**, grown from a 32/550 baseline, with the pre-recovery exact set pinned as `KNOWN_GOOD` regression guards). The same assignments are cross-checked against the ORIGINAL program's own per-atom types (the suite's OPTIMOL log, `MMFF94_opti.log` from the CCL archive): byte-identical on the 753-suite check (and on all 16 revision-affected molecules of the 761 suite), with the four remaining atoms parameter-inert (the metal-hydrate cations' oxidation-state types — the original program's FE+2/CU+1 vs OpenBabel's FE+3/CU+2 — and the dative sulfone-O class), all proven parameter-identical by the three-way energy checks; Tinker's prm atom table uses the same class numbering. Energy and charge coverage: per-term residuals vs BatchMin sit at ≤10⁻⁴ for all seven terms on the whole suite — bend/oop/vdw/torsion 761/761, strbnd 760/761, stretch 759/761 (the two generated-bond rows of the ERULE fragments sit within the reference's 3-decimal print precision), electrostatics 759/761 (the two type-76-anion exclusions, cross-verified three ways against Tinker's independent transcription on 2026-08-05; FAPLUD's q⁰(72) split closed 2026-08-07); the full accounting lives in `tests/VALIDATION.md`. Partial charges match to <10⁻³ on the 757-molecule charge-comparable set. The charged variants are structural: quaternary N (34), oxide O (35), oxenium O (51), iminium N (54), amidinium N (55), guanidinium N/C (56/57), pyridinium N (58), sp3 N-oxide N (68), anionic terminal S (72), sulfinate S (73), perchlorate Cl (77), imidazolium C/N (80/81), halide anions (89/90/91).

> **Note on third-party implementations.** The Halgren papers are the sole
> authoritative specification for MMFF94. We may study third-party implementations
> (OpenBabel, RDKit, TINKER, etc.) to understand difficulties or to numerically cross-check
> our output, but the code in this repository is written from the published papers,
> not adapted from other projects. Explicit references to other projects' source
> files, line numbers, or internal logic must **not** appear in our code comments
> or documentation — only their published parameter files and observable numerical
> output (e.g., `obenergy`, `obabel`) may be cited for cross-validation.

---

## Parameters — source and extraction

### Source

The numeric parameters are extracted from **OpenBabel**'s text-format `.par` files (`mmffbond.par`, `mmffang.par`, `mmffstbn.par`, `mmfftor.par`, `mmffvdw.par`, `mmffchg.par`, `mmffpbci.par`, `mmffoop.par`, `mmffdef.par`, `mmffprop.par`). These files contain the actual parameter values (k_b, r₀, k_a, θ₀, V_n, γ_n, R*, ε, G_i, α_i, N_i, BCI values, and k_oop) as plain text tables. The parameter values originate from Halgren's publications (1996–1999); the `.par` files are a mechanical reformatting of those values.

The extraction script reads them from `temp_ob/data/` — populate that directory from an OpenBabel install (e.g. copy `tests/scripts/.venv/Lib/site-packages/openbabel/bin/data/mmff*.par`).

The extraction script is a build-time Python utility that **converts format only** — it reads OpenBabel's `.par` text format and writes TypeScript literal objects. It does not implement any MMFF94 logic itself. The generated `.ts` files are committed to the repository so that `npm install && npm run build` works without Python.

### The extraction script (`scripts/extract-mmff94-par.py`)

```
python scripts/extract-mmff94-par.py
```

What the script does:

1. Reads each `.par` file in `temp_ob/data/`
2. Parses parameter sections by filename:
   - `mmffbond.par` → `bond.ts`
   - `mmffang.par` → `angle.ts`
   - `mmffstbn.par` → `stretch-bend.ts`
   - `mmfftor.par` → `torsion.ts`
   - `mmffvdw.par` → `van-der-waals.ts`
   - `mmffchg.par` + `mmffpbci.par` → `bci.ts`
   - `mmffoop.par` → `out-of-plane.ts`
3. For each section, emits a TypeScript file exporting a typed array or record:
   ```typescript
   // Example output for bond.ts
   export interface BondParams {
     k_b: number;   // force constant in mdyn/Å
     r0: number;    // equilibrium bond length in Å
   }

   /** Keyed by "priority-type1-type2" string. Wildcard types use 0. */
   export const BOND_PARAMS: Record<string, BondParams> = {
     "0-1-1":   { k_b: 4.258, r0: 1.508 },
     "0-1-2":   { k_b: 4.539, r0: 1.482 },
     // … ~500 entries
   };
   ```
   Note: the script converts format only — there is no validation step
   (duplicate-key or range checks) yet.

### Cross-checking

We will test against **Dr. Halgren's MMFF94 Validation Suite**
(https://server.ccl.net/cca/data/MMFF94/ — the November 1998 revision,
761 structures: 698 from the Cambridge Structural Database plus 63
small molecules and ions, including the ERULE_01–08 empirical-rule
fragments), with tabulated total energies,
individual energy component breakdowns, and optimized geometries
from the original MMFF94 development. This is the definitive
reference — Halgren's own test set, not a third-party implementation.
A local copy lives at `tests/fixtures/validation-suite/`.

We may also test against OpenBabel's `obenergy -ff MMFF94` or RDKit's `MMFF94GetEnergy()` and compare against `mmff94-ts`'s `calc_energy()`. But only after confirmation that the reference implementation is a faithful implementation of the original MMFF94 algorithm. 

---

## Implementation order

Work the modules in this sequence. Each step produces testable, shippable output.

### Phase 1 — Data model and I/O (1-2 days)

1. `src/types.ts` — `Atom`, `Bond`, `Molecule`, `TypedMolecule`
2. `src/sdf.ts` — parse MOL V2000 blocks from SDF files
3. `src/utils/vector.ts` — 3D vector math (add, subtract, cross, dot, norm, normalize, angle, rotate_around_axis)
4. `tests/sdf.test.ts` — round-trip test fixtures

### Phase 2 — Parameter extraction (1-2 days)

5. `scripts/extract-mmff94-par.py` — reads OpenBabel `.par` files, emits TS tables
6. Run it, commit generated `src/mmff94/parameters/*.ts`
7. Spot-check: open a few parameter files and verify they match Halgren's published values

### Phase 3 — Atom typing (3-5 days — hardest single piece)

8. `src/mmff94/assign-atom-types.ts` — the decision tree
9. `tests/atom-types.test.ts` — for each test fixture, verify every atom's assigned type matches OpenBabel's output

### Phase 4 — Energy terms (one day each)

10. `src/mmff94/energy/bond-stretch.ts` + test
11. `src/mmff94/energy/angle-bend.ts` + test
12. `src/mmff94/energy/stretch-bend.ts` + test (requires bond + angle working)
13. `src/mmff94/energy/out-of-plane.ts` + test
14. `src/mmff94/energy/torsion.ts` + test
15. `src/mmff94/energy/van-der-waals.ts` + test
16. `src/mmff94/energy/electrostatic.ts` (BCI charges live in `charges.ts` — `assign_bci_charges`, computed on demand) + test
17. `src/mmff94/energy/total.ts` — sum all terms, apply 1-4 scaling
18. `tests/reference-comparison.test.ts` + `tests/validate-against-suite.test.ts` — compare total energy against reference for all fixtures

### Phase 5 — Gradients (one day per term)

19. `src/mmff94/gradient/*.ts` — analytical derivatives, one file per term, matching `energy/` layout (DONE — Phase 5 complete)
20. `tests/gradient.test.ts` — finite-difference check on each term (δ = 1e-6 Å, relative error < 1e-5) (DONE)

### Phase 6 — Optimization (2-3 days)

21. `src/optimize/l-bfgs.ts` — L-BFGS with cubic line search (standard algorithm, nothing custom) (DONE — Nocedal & Wright Alg. 7.5, strong-Wolfe line search with cubic zoom)
22. `src/optimize/steepest-descent.ts` — simple fallback with Armijo line search (DONE — 15/16 fixtures converge at the spec from both starts; nicotine's vdW canyon is the documented boundary — descent-only there)
23. Test: optimize each fixture to max |gradient| < 0.05 kcal/mol/Å (DONE — L-BFGS 16/16, SD 15/16 from both SDF and perturbed starts; see `tests/optimization.test.ts`)

### Phase 7 — Polish and publish (1-2 days)

24. `bench/energy.bench.ts` — benchmark on drug-size molecules (200-500 atoms)
25. `examples/quickstart.ts` — clean example showing the full pipeline
26. `README.md` — API docs, quick start, references
27. Publish `0.1.0-alpha.1` to npm

### Phase 8 — Stretch goals (post-0.1)

28. **Per-interaction energy breakdowns** (see README "Stretch goals"):
    expose every component of every energy term — each bond stretch,
    angle bend, torsion, out-of-plane bend, stretch-bend cross term,
    vdW pair and electrostatic pair, with atoms, types, parameters and
    contribution — as first-class queryable data alongside the seven
    term totals. Every per-term mismatch found during development was
    diagnosed through OpenBabel's HIGH-verbosity interaction log
    (`tests/scripts/ob_energy_breakdown.py --verbose TERM`); the goal
    is to make that data simpler to access than OpenBabel's, where the
    per-interaction listing requires the verbose-log capture hack
    (buffered `std::cout`, space-separated headers, fd redirects).

---

## Testing philosophy

Every energy term is tested **in isolation** before it is tested in combination. This means a bug in bond stretching does not hide behind a bug in torsion.

| Test | What it checks | Method |
|---|---|---|
| **Unit** | A single energy function returns the right value for a known geometry | Compute by hand or with a reference tool for a 2-3 atom test case (e.g., H₂ for bond stretch, H₂O for angle bend) |
| **Regression** | Energy terms match reference values for the fixtures and the validation suite | `obenergy` logs in `tests/references/` (`reference-comparison.test.ts`) and BatchMin per-component energies (`validate-against-suite.test.ts`); OOP also matches several suite molecules exactly |
| **Typing coverage** | `assign_atom_types` matches OpenBabel's canonical types | `atom-types-suite.test.ts` against `mmff94-atom-types.json` (761 molecules — the full Nov-1998-revision suite; **761/761 type-exact — 100%**, with the exact set pinned as `KNOWN_GOOD` regression guards) |
| **Gradient** | Analytical dE/dx matches (E(x+δ) − E(x−δ)) / (2δ) | Finite-difference on each coordinate of each atom in each fixture: δ = 10⁻⁶ Å, relative error < 10⁻⁵ (DONE — worst observed error 8×10⁻⁸) |
| **Optimization** | After minimization, max \|gradient\| < threshold and total energy is lower than starting | L-BFGS on each fixture from the SDF geometry and from a perturbed geometry (16/16 fixtures converge at max\|g\| < 0.05 — see `tests/optimization.test.ts`) |
| **Browser** | The built `dist/` loads in headless Chromium (the system Edge, driven by playwright) and gives the same energies as the Node path to 1e-7 | `tests/browser-smoke.test.ts` — rebuilds dist with tsc, serves it over HTTP, runs the full pipeline in-page on 4 fixtures, and compares against the Node-side results and the reference logs |

Test fixtures are real molecules with known MMFF94 energies, stored as `.sdf` files in `tests/fixtures/sdf/`. Reference energies live in `tests/references/*.mmff94.log` (obenergy output).

---

## Documentation rules

Two rules keep the docs from drifting (the 753→761 upgrade proved
what happens without them):

1. **Never restate a generated number.** The validation census lives
   in `docs/validation/report.md`, regenerated by `npm run docs`. The
   prose docs (README, VALIDATION.md, compliance.md) point at it —
   they never repeat its counts. A hand-written number in a prose doc
   is a future staleness bug.

2. **No PR that changes a number ships without regenerating the
   report.** The generator is the source of truth; the committed report
   is a build artifact. `npm run docs:check` fails CI if the committed
   report differs from a fresh run — it runs in the CI step. If a PR
   changes any energy term, parameter, or typing rule, regenerate the
   report and commit the updated `docs/validation/` files with it.

Closure narratives (the "why" behind a code decision) belong in the
code comment of the file that implements it, with the commit message
that introduced it. The implementer's notes file is a map to those
commitments, not a parallel narrative.

---

## Integration into Valence (secondary)

After `mmff94-ts` is alpha-quality:

1. Valence adds `mmff94-ts` as an npm dependency
2. `web-vbvis/src/services/resolve3d.ts` gains a new path: when PubChem, CIR, and RDKit.js all fail, instead of calling `place3D()` (graph-walk), it calls:
   - `place3D()` for an initial guess (existing code)
   - `mmff94-ts`'s `optimize()` to refine that guess to an MMFF94-quality geometry
3. The fallback warning popup in `setStatus()` can be removed once the geometry quality is comparable to PubChem's
4. No `mmff94-ts` file depends on Valence or knows it exists

---

## MMFF94 Compliance Statement

The compliance statement lives at `docs/mmff94-compliance.md` (modeled on
Wavefun's https://downloads.wavefun.com/FAQ/MMFF94_compliance.html). The
full census (every per-term residual, total energy, charge count) is in
the generated `docs/validation/report.md` (`npm run docs` regenerates it) —
both documents point at the report for their numbers.

It certifies:

- **Which MMFF94 variant** — MMFF94 (not MMFF94s, not the MMFF94-like approximation
  some libraries ship under the MMFF94 name).
- **Energy terms** — all seven implemented with the correct Halgren functional forms.
- **Parameter provenance** — extracted from OpenBabel's `.par` files, cross-checked against Tinker's `mmff94.prm`.
- **Cross-check results** — per-term comparison against Halgren's 761-molecule suite (see `report.md` for the full census).
- **Limitations** — known deviations noted explicitly.
- **Dependencies** — zero runtime dependencies; compiles to ESM in `dist/`.

The implementer-level forensics live in the code comments of the
implementing files; the notes file (`docs/implementer-notes.md`) maps
to them.

---

## References

1. **Halgren, T. A.** *J. Comput. Chem.* 1996, *17*, 490–519. — MMFF94 definition, energy terms, and parametrization philosophy.
2. **Halgren, T. A.** *J. Comput. Chem.* 1996, *17*, 520–552. — Atom types, bond parameters, angle parameters, and stretch-bend parameters.
3. **Halgren, T. A.** *J. Comput. Chem.* 1996, *17*, 553–586. — Van der Waals parameters and electrostatic model (BCI).
4. **Halgren, T. A.** *J. Comput. Chem.* 1996, *17*, 616–641. — MMFF94s (the s = "static" variant for conformational energies).
5. **Halgren, T. A.** *J. Comput. Chem.* 1999, *20*, 720–729. — Torsion parameters.
6. **OpenBabel** — `mmff*.par` text-format parameter files in OpenBabel's data directory. Source for the extraction script (read from `temp_ob/data/`).
7. **Nocedal, J.; Wright, S. J.** *Numerical Optimization*, 2nd ed., Springer, 2006. — L-BFGS algorithm reference.
8. **Halgren MMFF94 Validation Suite** — https://server.ccl.net/cca/data/MMFF94/ — 761 structures (698 from the Cambridge Structural Database plus 63 small molecules and ions, November 1998 revision) with tabulated energies, component breakdowns, and optimized geometries from the original MMFF94 development.
9. **Wavefun MMFF94 Compliance Statement** — https://downloads.wavefun.com/FAQ/MMFF94_compliance.html — reference format for MMFF94 implementation certification.

---

## License note

`mmff94-ts` is MIT-licensed. The parameter data is extracted from OpenBabel (GPL-licensed) `mmff*.par` data files which contain only numerical facts derived from Halgren's published papers. The extraction script is a mechanical format converter. If licensing concerns arise, the parameter tables can be re-typed by hand from the original publications and the GPL dependency eliminated entirely.
