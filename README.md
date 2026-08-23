# mmff94-ts

**MMFF94 force field in pure TypeScript.** Runs in the browser and Node.js.
Zero dependencies. Validated against Halgren's 761-molecule suite.

[![tests](https://img.shields.io/badge/tests-245%20passed-brightgreen)](tests/)
[![validation](https://img.shields.io/badge/validation-761%2F761-blue)](docs/validation/report.md)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

```typescript
import { parse_sdf, calc_energy, optimize_lbfgs } from 'mmff94-ts';

const mol = parse_sdf(sdfText);
console.log(calc_energy(mol));   // every term + total, in kcal/mol
console.log(optimize_lbfgs(mol)); // MMFF94-minimized geometry
```

## Why

Existing JavaScript force fields either ship WebAssembly (big downloads,
sandbox-hostile), expose only a total-energy API (no per-term
breakdown), or refuse molecules without explicit hydrogens.
`mmff94-ts` is different: **pure TypeScript, all seven terms with the
correct Halgren functional forms, per-term and per-atom numeric
validation against the original program's own suite, cross-checked
three ways.**

## Accuracy

| Check | Result |
|---|---|
| Atom typing | **761/761** type-exact vs OpenBabel |
| Per-term residuals | **≤1e-4** kcal/mol on all seven terms |
| Partial charges | **<1e-3** e per atom on 757 molecules |
| Total energies | **758/761** within 1e-3 of BatchMin |
| Cross-checks | OpenBabel + Tinker, 16 molecules, three-way |

The full census — per-term tables, worst residuals, the ERULE
exceptions, and the AN11A/DOZNIP reference-inconsistency exclusions —
is generated from the suite files: see the
[Validation report](docs/validation/report.md).

## Features

- **Bond, angle, stretch-bend, torsion, vdW, electrostatic, out-of-plane** —
  all seven MMFF94 terms with the published Halgren factors
- **Analytical gradients** — finite-difference verified (worst 8e-8)
- **L-BFGS + steepest descent** — converge from perturbed starts
- **BCI partial charges** — from the bond-charge-increment model
- **SDF/MOL parser** — standalone, no dependencies
- **Browser + Node** — identical results (verified in headless Edge)
- **Zero runtime deps** — one `import` and you're done

## Quick Start

```bash
npm install mmff94-ts
```

```typescript
import { parse_sdf, calc_energy, optimize_lbfgs } from 'mmff94-ts';

const mol = parse_sdf(sdfText);
const energy = calc_energy(mol);
console.log(energy.total);       // -4.73435 (kcal/mol)
console.log(energy.torsion);     // per-term breakdown

const optimized = optimize_lbfgs(mol);
console.log(optimized.converged); // true
console.log(optimized.energy);    // energy at the minimum
```

See [`examples/quickstart.ts`](examples/quickstart.ts) for a complete
walkthrough — parsing SDF, assigning types, computing BCI charges,
and printing per-term energies.

## What people use it for

- **Interactive molecular viewers** — energy evaluation in the browser,
  no server round-trip
- **Conformer refinement** — L-BFGS minimization from an embedder's
  starting geometry
- **Teaching** — per-term energy decomposition that students can inspect
- **Fallback geometry** — when PubChem, CIR, and RDKit.js all fail
- **Tooling** — custom analysis pipelines that need per-term MMFF94
  energies in TypeScript

## Documentation

- [Walkthrough](docs/walkthrough.md) — the full pipeline, end to end
- [Validation report](docs/validation/report.md) — the complete census
- [Numerical precision](docs/numerical-precision.md) — error-budget math
- [Implementer's notes](docs/implementer-notes.md) — forensics and scripts

## Status

All seven terms, gradients, and optimization are implemented and gated
in `npm run test`. The public API is stable for the energy/gradient/optimize
pipeline; the pre-release tag reflects the documentation polish still in
progress.

## License

MIT.
