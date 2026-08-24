# mmff94-ts

**The MMFF94 force field, in pure TypeScript.**

No WebAssembly. No native builds. No server. Every term of Halgren's
force field — bond, angle, stretch-bend, torsion, out-of-plane, van der
Waals, electrostatic — implemented from the papers, validated per-term
against the original 761-molecule validation suite, and fast enough to
hold its own against the canonical Fortran.

[![tests](https://img.shields.io/badge/tests-336%20passed%20%2F%201%20documented%20skip-brightgreen)](tests/)
[![typing](https://img.shields.io/badge/atom_typing-761%2F761%20exact-blue)](docs/validation/report.md)
[![energies](https://img.shields.io/badge/per--term_%E2%89%A41e--4-759%2B%2F761-blue)](docs/validation/report.md)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

```typescript
import { parse_sdf, calc_energy, optimize_lbfgs } from 'mmff94-ts';

const mol = parse_sdf(sdfText);        // V2000 SDF/MOL in, molecule out
const e = calc_energy(mol);            // every term + total (kcal/mol)
const opt = optimize_lbfgs(mol);       // gradient-based minimization

console.log(e.bond_stretch, e.angle_bend, e.torsion, /* … */ e.total);
console.log(opt.converged, opt.energy.total, opt.final_max_gradient);
```

---

## Why this exists

JavaScript force fields to date have made one of two trades: ship the
science as an opaque WebAssembly blob (heavy, sandbox-hostile,
unreadable), or ship a toy (total energy only, no gradients, no
validation). Neither serves real work.

mmff94-ts takes a third road — **transcribe the science faithfully and
prove it**:

- **All seven terms**, each with Halgren's published functional form:
  the buffered 14-7 vdW, the double-cubic stretch-bend, the
  three-term torsion Fourier series, the Wilson out-of-plane, the
  cubic bond and sextic angle with CB anharmonicity, and
  charge–charge electrostatics with the ×0.75 1–4 scaling.
- **Analytical gradients for every term**, verified against central
  finite differences.
- **Numbers you can check.** Not benchmarked "on typical molecules" —
  compared, term by term, against BatchMin's own output on Halgren's
  full validation suite.

## Validation — measured, not claimed

Every number below is generated from the suite files and re-checked on
every test run; nothing here is transcribed by hand.

| Check | Result |
|---|---|
| Atom typing | **761/761 exact** vs OpenBabel's MMFF94 typer |
| Per-term energies | all seven terms ≤1e-4 kcal/mol vs BatchMin (759/761 within that gate; the rest are documented ERULE print-precision artifacts) |
| Partial charges | ≤1e-3 e⁻ on **757/757** comparable molecules |
| Total energies | 758/761 within 1e-3 kcal/mol |
| Gradient correctness | finite-difference worst relative error 4e-8 over drug-like fixtures (5e-7 on the 304-atom trp-cage) |
| Independent implementations | final energies match RDKit MMFF94 to ~0.5 kcal/mol and Tinker 26.2 to four decimals on shared minima |

Three details worth knowing, because they're where hand-waved
implementations hide:

- Two suite entries (**AN11A**, **DOZNIP**) are *excluded* — Halgren
  himself flagged the type-76 anionic-nitrogen reference as unreliable.
  We document rather than average away.
- Two **ERULE** residuals sit above 1e-4 only because BatchMin printed
  the empirical-rule parameters at 3 decimal places. The gate is the
  reference's print precision, not our arithmetic.
- A third-party bug surfaced during benchmarking: Tinker 26.2 assigns a
  phantom −1 e net charge to neutral sulfones (hard-coded −0.5 base
  charge on type 107 with no sulfone compensation). We filed it
  upstream ([TinkerTools/tinker#185](https://github.com/TinkerTools/tinker/issues/185))
  rather than quietly absorbing it into our comparison tables.

The complete census — every residual, every exception, every exclusion
with its reason — lives in the
[validation report](docs/validation/report.md).

## Performance

The hot path is compiled once per molecule: every interaction parameter
resolved into flat typed arrays, then energy and gradient evaluated
allocation-free — no objects, no maps, no garbage collector.

On 29 drug-like molecules (OpenFF Industry set, 25–70 heavy+H atoms,
QM-start geometries), minimizing to convergence:

| engine | converged | median wall time |
|---|---|---|
| **mmff94-ts** — L-BFGS, pure JS | **26/29** | **242 ms** |
| Tinker 26.2 `minimize` — Fortran | 29/29 | 296 ms |

Same order as forty years of Fortran, in an interpreted language, on a
2014-era desktop core. Full methodology, caveats, and per-molecule data:
[docs/benchmark.md](docs/benchmark.md).

## Quick start

```bash
npm install mmff94-ts        # not yet published — see Status
# meanwhile: git clone && npm install && npm run build
```

```typescript
import { parse_sdf, assign_atom_types, assign_bci_charges, calc_energy, optimize_lbfgs } from 'mmff94-ts';

// Rich path: type and charge explicitly, inspect everything.
const mol = assign_bci_charges(assign_atom_types(parse_sdf(sdfText)));
const energy = calc_energy(mol);
console.log(energy.total);      // kcal/mol, with .bond_stretch, .angle_bend,
                                // .stretch_bend, .torsion, .out_of_plane,
                                // .van_der_waals, .electrostatic beside it

const result = optimize_lbfgs(mol);
console.log(result.converged);            // true
console.log(result.energy.total);         // energy at the minimum
console.log(result.final_max_gradient);   // convergence quality
console.log(result.molecule.atoms);       // the minimized geometry

// Simple path: bare SDF text straight to a minimized molecule.
const done = optimize_lbfgs(parse_sdf(sdfText));
```

Convergence defaults to max |gᵢ| < 0.05 **or** RMS gradient < 0.02
kcal/mol/Å (TINKER-style dual gate). `criterion: 'max'` restores the
strict single-coordinate rule; `'rms'` matches TINKER exactly.

A complete walkthrough — parsing, typing, charging, per-term analysis —
is in [`examples/quickstart.ts`](examples/quickstart.ts), with the
pipeline documented end-to-end in the [walkthrough](docs/walkthrough.md).

## API surface

| Function | Purpose |
|---|---|
| `parse_sdf(text)` | V2000 SDF/MOL → `Molecule` (bond-index-safe on malformed input) |
| `assign_atom_types(mol)` | MMFF94 symbolic typing — aromaticity-aware, valence-driven |
| `assign_bci_charges(typed)` | Bond-charge-increment partial charges (eq. 15 sharing) |
| `calc_energy(typed)` | Per-term + total energy |
| `calc_gradient(typed)` | Flat analytical gradient (3·N) |
| `optimize_lbfgs(mol, opts?)` | Limited-memory BFGS with strong-Wolfe line search |
| `optimize_steepest_descent(mol, opts?)` | Robust fallback minimizer |
| `parameter_gap_report(typed)` | Atoms outside MMFF94's parameter space (hypervalent centers, untyped elements) — diagnostic, never silently wrong |

Everything is plain data in and plain data out: no classes to instantiate,
no state to manage, every function pure over its inputs.

## Documentation

| Document | Contents |
|---|---|
| [Walkthrough](docs/walkthrough.md) | the full pipeline end to end |
| [Validation report](docs/validation/report.md) | the complete census — every residual, exception, exclusion |
| [Benchmark](docs/benchmark.md) | optimizer methodology + per-molecule data |
| [Numerical precision](docs/numerical-precision.md) | error-budget math |
| [Implementer's notes](docs/implementer-notes.md) | resolved forensics: the trp-cage zwitterion, ERULE provenance, reference inconsistencies |

## Status

**Implemented and gated:** all seven energy terms, analytical gradients,
BCI charges, L-BFGS and steepest-descent optimization, aromatic-ring
perception, hypervalent-center diagnostics. 336 tests, 1 intentional skip.

**Known limitations, stated plainly:**

- Explicit hydrogens required (MMFF94's own contract).
- All-pairs nonbonded evaluation — O(N²) is fine at drug scale (≤ ~500
  atoms), not protein scale. Cutoffs + neighbor lists are designed-for
  but not built.
- Single-geometry minimization; no conformer ensemble generation yet.
- Reads V2000 SDF/MOL; MOL2 and PDB readers not yet provided.

## License

MIT. The MMFF94 functional forms are those published by Thomas A.
Halgren in J. Comput. Chem. 17, 490–641 (1996); if you use this library
for scientific work, please cite those papers alongside this repository.
