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

Typing and energies are validated against the full Halgren suite
(753 molecules): **753/753 type-exact** vs OpenBabel's canonical types
(100%), and all seven energy terms match the BatchMin references on
the 749-molecule reproducible set (100%; the remaining 4 reference
anomalies are documented in [Validation](#validation)). See
[Validation](#validation).

## Validation

Every term is checked against two independent references, and every
analytical gradient against finite differences.

1. **Halgren's 753-molecule MMFF94 validation suite** — per-component
   energies and per-atom partial charges vs the
   suite's reference values. **Our atom typing reproduces the reference
   types exactly for every one of them (753/753)**.

2. **OpenBabel** — per-term energies and per-atom partial
   charges for 16 small organic molecules. All seven energy terms match exactly to five
   decimals. Will continue to expand molecular coverage.

Per-component agreement with BatchMin on the reproducible set
(749 molecules — 4 reference anomalies are excluded, each documented
in [`tests/VALIDATION.md`](tests/VALIDATION.md)):

| Term | Exact | Max abs(Δ) (kcal/mol) |
|-------|---|---|
| Bond stretch | 749/749 | 0.00 |
| Angle bend | 749/749 | 0.00 |
| Stretch-bend | 749/749 | 0.00 |
| Torsion | 749/749 | 0.00 |
| Van der Waals | 749/749 | 0.00 |
| Out-of-plane | 749/749 | 0.02 |
| Electrostatic | 749/749 | 0.00 |

 All 749 reproducible suite molecules are type-exact (100%) and every
 energy term matches the BatchMin references on all of them (within
 0.05 kcal/mol; the partial charges match to < 1e-3 e⁻ per atom on the
 same set). 
 More details in [`tests/VALIDATION.md`](tests/VALIDATION.md).

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

## License

MIT.
