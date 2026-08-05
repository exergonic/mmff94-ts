# mmff94-ts

**Merck Molecular Force Field (MMFF94) — pure TypeScript, zero runtime dependencies.**

> ⚠️ **Work in progress — pre-alpha.**
> The public API is not yet stable and may change before 0.1.0. See the [Status](#status)
> table for the full picture.

Energy evaluation for organic molecules, running in the browser or
Node.js without WebAssembly, native binaries, or external services.

## Why

A pure TypeScript force field means your molecular mechanics code
bundles trivially with any web framework, works in sandboxed
environments (no WASM, no native binaries), and composes naturally
with the rest of your TypeScript toolchain.

`mmff94-ts` implements the complete MMFF94 energy functional form —
the buffered 14-7 van der Waals, stretch-bend cross term,
Fourier-series torsion, bond-charge-increment electrostatics, and
the rest of the seven terms — so you can run energy evaluations
entirely on the client side.

## Status

| Term | Status |
|---|---|
| Bond stretch | ✅ |
| Angle bend | ✅ |
| Stretch-bend (class-II cross term) | ✅ |
| Torsion (Fourier series) | ✅ |
| Van der Waals (buffered 14-7) | ✅ |
| Electrostatic (BCI model) | ✅ |
| Out-of-plane bending | ✅ |
| 1-4 scaling | ✅ |
| Analytical gradients | ✅ |
| Optimization (L-BFGS + steepest descent) | ✅ |


## Validation

We checked every energy term against two independent references.
We checked every analytical gradient against finite differences.

1. **Halgren's 753-molecule MMFF94 validation suite.** This is the
   reference for the library. It gives the per-component energies and
   the per-atom partial charges for each molecule. We compared every
   term with the reference values. All comparable terms match within
   0.0001 kcal/mol:

| Term | Molecules that match | Worst \|Δ\|  (kcal/mol) |
|---|------|---|
| Bond stretch | 753 of 753 | 5.0e-5 |
| Angle bend | 753 of 753 | 3.4e-5 |
| Stretch-bend | 753 of 753 | 4.3e-5 |
| Torsion | 753 of 753 | 4.7e-5 |
| Out-of-plane | 753 of 753 | 1.6e-5 |
| Van der Waals | 751 of 753 | 4.4e-5 |
| Electrostatic | 751 of 753 | 6.8e-5 |

   The atom types match the reference types for all 753 molecules.
   The partial charges match the reference values to 0.001 e per atom
   on 749 molecules. The four remaining molecules have one term where
   the reference itself is inconsistent. The
   [Validation document](tests/VALIDATION.md) has the full details.

   The hard evidence is committed: [`docs/validation/total-energies.txt`](docs/validation/total-energies.txt)
   lists all 753 totals side by side with the suite's own OPTIMOL and
   BatchMin values (749/753 to 0.001 kcal/mol; the four exceptions are
   the documented reference anomalies), and
   [`docs/validation/per-term-and-charges.txt`](docs/validation/per-term-and-charges.txt)
   lists the per-term and per-atom charge deltas molecule by molecule.
   Regenerate both with `npm run validation:doc`.

2. **OpenBabel.** Per-term energies and per-atom partial charges for
   16 small organic molecules. All values match to five decimal
   places.


## Usage

```typescript
import { parse_sdf, calc_energy, optimize_lbfgs } from 'mmff94-ts';

const mol = parse_sdf(sdfText);
const energy = calc_energy(mol);     // every term, plus the total
const optimized = optimize_lbfgs(mol);
console.log(optimized, energy.total, energy.bond_stretch, energy.torsion);
```

See [`examples/quickstart.ts`](examples/quickstart.ts) for a complete
walkthrough — parsing SDF, assigning types, computing BCI charges,
and printing per-term energies.

## Documentation

- **[`docs/walkthrough.md`](docs/walkthrough.md)** — traces the full
  pipeline from raw SDF to energy components: data model, geometry
  primitives, atom typing, parameter lookup, every energy term's
  functional form, gradient layout, and optimization strategy.
- **[`docs/numerical-precision.md`](docs/numerical-precision.md)** —
  addresses whether a JavaScript MMFF94 can match C++ reference
  energies. Contains the error-budget math (IEEE 754 doubles,
  accumulation analysis, Kahan summation) and recommended validation
  tolerances.
- **[`docs/implementer-notes.md`](docs/implementer-notes.md)** — the
  forensics behind the validation claims: the numbering systems, the
  closure narratives, the per-anomaly numbers, and the commands that
  regenerate every number.

## License

MIT.
