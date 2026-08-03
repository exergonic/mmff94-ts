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
> geometry optimization is working (14/16 fixtures converge to
> max |gradient| < 0.05 kcal/mol/Å from both the SDF and perturbed
> geometries — see `tests/optimization.test.ts`). The steepest-descent
> fallback converges 15/16 fixtures at the spec (nicotine's vdW canyon
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
| Geometry optimization (L-BFGS) | ✅ (16/16 fixtures to max\|g\| < 0.05) |
| Geometry optimization (SD fallback) | ✅ (15/16 fixtures; nicotine's vdW canyon is the documented boundary) |

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
decimals) on all 16 fixtures** — formamide included since the amide-N
typing (types 10/28) landed — and every fixture **total** matches
exactly. The BCI partial charges also match the reference logs per
atom (`tests/charges.test.ts`).

| Term | Exact on fixtures | Notes |
|---|---|---|
| Bond stretch | 16/16 | |
| Angle bend | 16/16 | |
| Stretch-bend | 16/16 | |
| Torsion | 16/16 | |
| Van der Waals | 16/16 | |
| Electrostatic | 16/16 | per-atom charges pinned in charges.test.ts |
| Out-of-plane | 16/16 | 0 on most fixtures; see BatchMin table below |

Per-component energies vs BatchMin on the 123 typing-exact suite
molecules: bond/angle/vdW/oop/strbnd/torsion all 123/123 (the three
residuals were real bugs — degenerate i = l "torsions" in FUVDOP's
3-ring, FILNOD's 5-ring torsions classed by atom flags instead of ring
aromaticity, and an entry OpenBabel's strbnd transcription lost, JIYJAC —
all fixed), electrostatic 140/140 and angle 140/140 — the formal-charge model
(part V eq. 15: primary charges + negative-charge sharing) reproduces
the reference partial charges per atom, and the sulfinate S=O angles
resolved into the 32-keyed parameters (the sulfinate S=O is typed 7 by
the reference typing rules but keyed 32 in every parameter table —
TINKER's mmff94.prm and OpenChemLib's angle.csv confirm the same five
73-center entries).

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
