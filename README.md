# mmff94-ts

**Merck Molecular Force Field (MMFF94) — pure TypeScript, zero runtime dependencies.**

Energy evaluation, analytical gradients, and geometry optimization
for organic molecules, running in the browser or Node.js without
WebAssembly, native binaries, or external services.

## Why

A pure TypeScript force field means your molecular mechanics code
bundles trivially with any web framework, works in sandboxed
environments (no WASM, no native binaries), and composes naturally
with the rest of your TypeScript toolchain.

`mmff94-ts` implements the full MMFF94 functional form — including
the buffered 14-7 van der Waals, bond charge increment electrostatics,
stretch-bend cross term, and Fourier-series torsion — so you can
run energy evaluations and geometry optimizations entirely on the
client side.

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
