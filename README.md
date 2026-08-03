# mmff94-ts

**Merck Molecular Force Field (MMFF94) — pure TypeScript, zero runtime dependencies.**

> ⚠️ **Work in progress — pre-alpha.**
>
> All seven energy terms (bond stretch, angle bend, stretch-bend,
> torsion, buffered 14-7 van der Waals, electrostatic (BCI),
> out-of-plane) are implemented and validated against OpenBabel
> `obenergy` logs and Halgren's 753-molecule validation suite, the
> analytical gradients of every term are cross-checked against
> finite differences (worst relative error 8×10⁻⁸), and L-BFGS
> geometry optimization is working (16/16 test molecules converge to
> max |gradient| < 0.05 kcal/mol/Å from both the SDF and perturbed
> geometries — see `tests/optimization.test.ts`). The steepest-descent
> fallback converges 15/16 at the spec (nicotine's vdW canyon
> is its documented boundary). The public API is not yet stable and
> may change before 0.1.0. See the [Status](#status) table for the
> full picture.

Energy evaluation for organic molecules, running in the browser or
Node.js without WebAssembly, native binaries, or external services.
Analytical gradients and geometry optimization (L-BFGS primary,
steepest-descent fallback) are complete.

## Why

A pure TypeScript force field means your molecular mechanics code
bundles trivially with any web framework, works in sandboxed
environments (no WASM, no native binaries), and composes naturally
with the rest of your TypeScript toolchain.

`mmff94-ts` implements the complete MMFF94 energy functional form —
the buffered 14-7 van der Waals, stretch-bend cross term,
Fourier-series torsion, bond-charge-increment electrostatics, and
the rest of the seven terms — so you can run energy evaluations
entirely on the client side. Geometry optimization via L-BFGS is
working (see [Status](#status)).

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
| 1-4 scaling | ✅ (electrostatic ×0.75, inside the term) |
| Analytical gradients | ✅ (all 7 terms, FD-verified < 10⁻⁵) |
| Geometry optimization (L-BFGS) | ✅ (16/16 at max\|g\| < 0.05) |
| Geometry optimization (SD fallback) | ✅ (15/16; nicotine's vdW canyon is the documented boundary) |

## Validation

Every term is checked against two independent references, and every
analytical gradient against finite differences.

1. **OpenBabel `obenergy`** — per-term energies and per-atom partial
   charges for 16 small organic molecules (`tests/reference-comparison.test.ts`,
   `tests/charges.test.ts`). All seven terms match exactly to five
   decimals.
2. **Halgren's 753-molecule MMFF94 validation suite** — per-component
   energies vs BatchMin 5.5, and per-atom partial charges vs the
   suite's reference values (`tests/validate-against-suite.test.ts`,
   `tests/charges-suite.test.ts`). The comparison runs on the 140/550
   suite molecules whose atom typing reproduces the reference types
   exactly (pinned as `KNOWN_GOOD` in `tests/atom-types-suite.test.ts`),
   so a component delta can never be blamed on a typing gap.

Per-component agreement with BatchMin on those 140 molecules:

| Term | Exact | Max Δ (kcal/mol) |
|-------|---|---|
| Bond stretch | 140/140 | 0.00 |
| Angle bend | 140/140 | 0.02 |
| Stretch-bend | 140/140 | 0.00 |
| Torsion | 140/140 | 0.00 |
| Van der Waals | 140/140 | 0.00 |
| Out-of-plane | 140/140 | 0.01 |
| Electrostatic | 140/140 | 0.00 |

Partial charges reproduce the suite's reference values to < 10⁻³ e⁻
per atom on 138/140 molecules (the two thiosulfinate anions are
excluded — BatchMin's dative representation). Gradients are checked
against central finite differences on every atom of every reference
molecule (`tests/gradient.test.ts`, δ = 10⁻⁶ Å, worst relative error
8×10⁻⁸).

The full molecule-by-molecule ledger — every delta and the known open
questions — lives in [`tests/VALIDATION.md`](tests/VALIDATION.md).

## Usage

```typescript
import { parse_sdf, calc_energy, optimize_lbfgs } from 'mmff94-ts';

// The simple path: a parsed molecule is all you need — typing and
// charges happen on demand, and the full per-term breakdown comes back.
const mol = parse_sdf(sdfText);
const energy = calc_energy(mol);     // every term, plus the total
console.log(energy.total, energy.bond_stretch, energy.torsion);

// Geometry optimization in one call — no preparation, no callback.
const result = optimize_lbfgs(mol);  // typed/charged molecule at the
                                     // minimum rides along in result.molecule
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
