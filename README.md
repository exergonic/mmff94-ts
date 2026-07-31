# mmff94-ts

**Merck Molecular Force Field (MMFF94) — pure TypeScript, zero runtime dependencies.**

> ⚠️ **Work in progress — pre-alpha.**
>
> The six implemented energy terms (bond stretch, angle bend,
> stretch-bend, torsion, buffered 14-7 van der Waals, out-of-plane) are
> validated against OpenBabel `obenergy` logs and Halgren's 753-molecule
> validation suite. The remaining pieces — the electrostatic (BCI) term,
> 1-4 scaling, analytical gradients, and geometry optimization — are
> stubs that return zeros. The public API is not yet stable and may
> change before 0.1.0. See the [Status](#status) table for the full
> picture.

Energy evaluation for organic molecules, running in the browser or
Node.js without WebAssembly, native binaries, or external services.
Analytical gradients and geometry optimization are in progress.

## Why

A pure TypeScript force field means your molecular mechanics code
bundles trivially with any web framework, works in sandboxed
environments (no WASM, no native binaries), and composes naturally
with the rest of your TypeScript toolchain.

`mmff94-ts` implements the MMFF94 energy functional form — the
buffered 14-7 van der Waals, stretch-bend cross term, Fourier-series
torsion, and the rest of the six implemented terms — so you can run
energy evaluations entirely on the client side. Bond charge increment
electrostatics, 1-4 scaling, analytical gradients, and geometry
optimization are still in progress (see [Status](#status)).

## Status

| Term | Status |
|---|---|
| Bond stretch | ✅ |
| Angle bend | ✅ |
| Stretch-bend (class-II cross term) | ✅ |
| Torsion (Fourier series) | ✅ |
| Van der Waals (buffered 14-7) | ✅ |
| Electrostatic (BCI model) | ⬜ |
| Out-of-plane bending | ✅ |
| 1-4 scaling | ⬜ |
| Analytical gradients | ⬜ |
| Geometry optimization (L-BFGS / SD) | ⬜ |

## Usage

```typescript
import { parse_sdf, assign_atom_types, calc_energy } from 'mmff94-ts';

const mol = parse_sdf(sdfText);
const typed = assign_atom_types(mol);
const energy = calc_energy(typed);

console.log(energy.total);           // kcal/mol
console.log(energy.bond_stretch);    // per-component breakdown
console.log(energy.torsion);
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
