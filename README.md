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
| Geometry optimization (L-BFGS) | ✅ |
| Geometry optimization (SD fallback) | ✅ (15/16; nicotine exception) |

## Validation

Every term is checked against two independent references, and every
analytical gradient against finite differences.

1. **OpenBabel** — per-term energies and per-atom partial
   charges for 16 small organic molecules. All seven energy terms match exactly to five
   decimals. Will continue to expand molecular coverage.
2. **Halgren's 753-molecule MMFF94 validation suite** — per-component
   energies and per-atom partial charges vs the
   suite's reference values. The comparison runs on all 550
   suite molecules OpenBabel can set up — **our atom typing now
   reproduces the reference types exactly for every one of them
   (550/550)**.

Per-component agreement with BatchMin on those 550 molecules:

| Term | Exact | Max abs(Δ) (kcal/mol) |
|-------|---|---|
| Bond stretch | 550/550 | 0.00 |
| Angle bend | 550/550 | 0.00 |
| Stretch-bend | 550/550 | 0.00 |
| Torsion | 550/550 | 0.00 |
| Van der Waals | 550/550 | 0.00 |
| Out-of-plane | 550/550 | 0.02 |
| Electrostatic | 550/550 | 0.00 |

 All 550 suite molecules are type-exact (100%) and every energy term
 matches the BatchMin references on all of them (within 0.05 kcal/mol).
 The "parameter-gap" workstream closed itself: the 18 remaining
 mismatches were all lookup or constant bugs in our transcription —
 the oop term's missing EqLvl3 step-down chain, the stretch-bend
 class-0 lookup scanning other classes, and the rounded cubic-bend
 constant (−0.007 vs BatchMin's precise −0.4·π/180).

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
