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

Every energy term is checked against Halgren's 761-molecule MMFF94
validation suite (November 1998 revision), with the per-term residuals
gated at ≤1e-4 kcal/mol in `npm run test`
(`tests/compliance-gate.test.ts`). Analytical gradients are
finite-difference checked. The full census — per-term tables, worst
residuals, coarse-precision exceptions, and the AN11A/DOZNIP
reference-inconsistency exclusions — lives in the generated
**[Validation report](docs/validation/report.md)** (`npm run docs`
regenerates it from the suite files).

Independent cross-checks: OpenBabel and Tinker on 16 small
organic molecules (per-term and per-atom charges). See the
[Implementer's notes](docs/implementer-notes.md) for the three-way
fixture comparison.


## Usage

```typescript
import { parse_sdf, calc_energy, optimize_lbfgs } from 'mmff94-ts';

const mol = parse_sdf(sdfText);
const energy = calc_energy(mol);     // every term, plus the total
const optimized = optimize_lbfgs(mol);
console.log(optimized, energy.total, energy.bond_stretch, energy.torsion);
```

MMFF94 is evaluated in vacuo: the dielectric is the default D = 1.0.

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
