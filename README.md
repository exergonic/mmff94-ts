# mmff94-ts

**Merck Molecular Force Field (MMFF94) — pure TypeScript, zero runtime dependencies.**

> ⚠️ **Work in progress — pre-alpha.**
>
> All seven energy terms (bond stretch, angle bend, stretch-bend,
> torsion, buffered 14-7 van der Waals, electrostatic (BCI),
> out-of-plane) are implemented and validated against OpenBabel
> `obenergy` logs and Halgren's 753-molecule validation suite, and
> the analytical gradients of every term are cross-checked against
> finite differences (worst relative error 8×10⁻⁸). The remaining
> piece — geometry optimization — is a stub. The public API is not
> yet stable and may change before 0.1.0. See the [Status](#status)
> table for the full picture.

Energy evaluation for organic molecules, running in the browser or
Node.js without WebAssembly, native binaries, or external services.
Analytical gradients are complete; geometry optimization is in
progress.

## Why

A pure TypeScript force field means your molecular mechanics code
bundles trivially with any web framework, works in sandboxed
environments (no WASM, no native binaries), and composes naturally
with the rest of your TypeScript toolchain.

`mmff94-ts` implements the complete MMFF94 energy functional form —
the buffered 14-7 van der Waals, stretch-bend cross term,
Fourier-series torsion, bond-charge-increment electrostatics, and
the rest of the seven terms — so you can run energy evaluations
entirely on the client side. Geometry optimization is still in
progress (see [Status](#status)).

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
| Geometry optimization (L-BFGS / SD) | ⬜ |

## Validation

Every implemented term is checked against reference energies: the
SDF fixtures against OpenBabel `obenergy` logs
(`tests/reference-comparison.test.ts`), and the out-of-plane and
electrostatic terms additionally against Halgren's own BatchMin
energies from the 753-molecule validation suite
(`tests/validate-against-suite.test.ts`). Every analytical gradient
is checked against central finite differences on every atom of every
fixture (`tests/gradient.test.ts`, δ = 10⁻⁶ Å, relative error
< 10⁻⁵ — worst observed 8×10⁻⁸).

All seven terms match the obenergy references **exactly (to 5
decimals) on all 15 asserted fixtures** (formamide is a documented
typing-gap skip until amide N typing lands), and every fixture
**total** matches exactly. The BCI partial charges also match the
reference logs per atom (`tests/charges.test.ts`).

| Term | Exact on fixtures | Notes |
|---|---|---|
| Bond stretch | 15/15 | |
| Angle bend | 15/15 | |
| Stretch-bend | 15/15 | |
| Torsion | 15/15 | |
| Van der Waals | 15/15 | |
| Electrostatic | 15/15 | per-atom charges pinned in charges.test.ts |
| Out-of-plane | 15/15 | 0 on most fixtures; see BatchMin table below |

Per-component energies vs BatchMin on the 91 typing-exact suite
molecules: bond/angle/vdW/oop 91/91, strbnd 91/91, torsion 91/91 (the
three residuals were real bugs — degenerate i = l "torsions" in FUVDOP's
3-ring, FILNOD's 5-ring torsions classed by atom flags instead of ring
aromaticity, and an entry OpenBabel's strbnd transcription lost, JIYJAC —
all fixed), electrostatic 89/91 (two metal-carboxylate salts awaiting
formal-charge input).

### Out-of-plane vs BatchMin (8 suite molecules, kcal/mol)

| Molecule | Ours | BatchMin | Δ |
|---|---|---|---|
| DADDAN | 0.255548 | 0.255547 | 0.000000 |
| GIDJUY | 0.216936 | 0.216938 | −0.000002 |
| VEJWOW | 0.176902 | 0.177154 | −0.000252 |
| DIKGAF | 0.160155 | 0.158925 | +0.001230 |
| FAXVAB | 0.127921 | 0.126658 | +0.001263 |
| GEXGIZ | 0.122862 | 0.123820 | −0.000958 |
| VIRBON | 0.101801 | 0.102969 | −0.001167 |
| AMHTAR01 | 0.203026 | 0.224486 | −0.021460 |

The full molecule-by-molecule ledger — every delta, the known open
questions, and how to update it — lives in
[`tests/VALIDATION.md`](tests/VALIDATION.md).

## Usage

```typescript
import { parse_sdf, assign_atom_types, compute_bci_charges, calc_energy } from 'mmff94-ts';

const mol = parse_sdf(sdfText);
const typed = assign_atom_types(mol);
compute_bci_charges(typed);       // BCI partial charges (the electrostatic
                                  // term also computes them on demand)
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
