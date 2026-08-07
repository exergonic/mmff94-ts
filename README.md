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

We checked every energy term against three references: the original
validation suite, OpenBabel, and Tinker.
We checked every analytical gradient against finite differences.

1. **Halgren's 761-molecule MMFF94 validation suite** (November 1998
   revision, from https://server.ccl.net/cca/data/MMFF94/). This is the
   reference for the library. It gives the per-component energies and
   the per-atom partial charges for each molecule. We compared every
   term with the reference values:

| Term | ≤1e-4 | Worst \|Δ\|  (kcal/mol) |
|---|------|---|
| Bond stretch | 759 of 761 | 2.0e-3 (ERULE_03, generated P–Si at the reference's print precision) |
| Angle bend | 761 of 761 | 4.3e-5 |
| Stretch-bend | 760 of 761 | 3.4e-4 (ERULE_03, inherited from the P–Si) |
| Torsion | 761 of 761 | 4.7e-5 |
| Out-of-plane | 761 of 761 | 1.6e-5 |
| Van der Waals | 761 of 761 | 4.4e-5 |
| Electrostatic | 759 of 761 | 7.1e-5 |

   The atom types match the reference types for all 761 molecules.
   The partial charges match the reference values to 0.001 e per atom
   on 757 molecules. The two AN11A/DOZNIP electrostatics exclusions
   are terms where the reference itself is inconsistent (FAPLUD's
   q⁰(72) split closed 2026-08-07 — see
   `docs/implementer-notes.md` §5.5). The
   [Validation document](tests/VALIDATION.md) has the full details.

   The tabulated total energies [total-energies.txt](docs/validation/total-energies.txt), 
   lists all 761 totals side by side with the suite's own OPTIMOL and
   BatchMin values (the largest residual on the matching set is
   ~8.0e-5). The exceptions are the documented
   AN11A/DOZNIP electrostatics anomalies, not shortcomings of the
   implementation.
   See the [Validation document](tests/VALIDATION.md) for details.

2. **OpenBabel.** Per-term energies and per-atom partial charges for
   16 small organic molecules. Most values match to five decimal
   places. The angle term differs by up to 0.0007 kcal/mol for three
   molecules. OpenBabel uses a rounded constant in the angle formula.
   See the [Implementer's notes](docs/implementer-notes.md).

3. **Tinker.** A second independent implementation. For the same 16
   molecules, total energies match to four decimal places. This is
   the print precision of Tinker. One angle term differs by 0.0001
   kcal/mol. See the [Implementer's notes](docs/implementer-notes.md).


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
