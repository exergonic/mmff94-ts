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
(100%), and **all seven energy terms are machine-exact (747/747 at
≤10⁻⁴ kcal/mol per term)** vs BatchMin's per-component references —
including the suite's only empirical-rule bond (OHMW1, closed at
1.4×10⁻⁶). Partial charges match the reference to <10⁻³ e⁻ on the
749-molecule charge-comparable set. The documented reference
anomalies are in [`tests/VALIDATION.md`](tests/VALIDATION.md). See
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

Per-component agreement with BatchMin on the reproducible set —
every molecule within 10⁻⁴ kcal/mol per term (the residual census,
`tests/scripts/residual-census.ts`; the 0.05 gate in the suite tests
is looser than this):

| Term | Molecules ≤10⁻⁴ | Worst |Δ| (kcal/mol) |
|-------|---|---|
| Bond stretch | 747/747 | 5.0e-5 (DEWJEU) |
| Angle bend | 747/747 | 3.4e-5 (BEVJER10) |
| Stretch-bend | 747/747 | 0.0 |
| Torsion | 747/747 | 1.4e-5 (BEVJER10) |
| Van der Waals | 747/747 | 4.4e-5 (MG2PW3) |
| Out-of-plane | 747/747 | 1.6e-5 (ARGIND11) |
| Electrostatic | 747/747 | 6.8e-5 (DONFOB) |

Atom typing is type-exact on all 753 suite molecules (100%), and the
partial charges match the BatchMin references to <10⁻³ e⁻ per atom on
the 749-molecule charge-comparable set. The suite's only
empirical-rule bond (OHMW1's hydroxide O–H, generated from part V
eqs. 18-19) matches to 1.4×10⁻⁶. Every empirical rule (bond, angle,
torsion, BCI fallback) is additionally pinned by hand-computed unit
tests. More details in [`tests/VALIDATION.md`](tests/VALIDATION.md).

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
