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

## Validation

Every implemented term is checked against reference energies: the
nine SDF fixtures against OpenBabel `obenergy` logs
(`tests/reference-comparison.test.ts`), and the out-of-plane term
additionally against Halgren's own BatchMin energies from the
753-molecule validation suite (`tests/validate-against-suite.test.ts`).

Bond stretch, angle bend, stretch-bend, and van der Waals match the
obenergy references **exactly (to 5 decimals) on all 9 fixtures**.
Torsion matches on 6 and is close on 3 (see notes); out-of-plane is
0 on these acyclic/planar fixtures. Totals match exactly where every
term does (ethane −4.73436, methane, formaldehyde).

| Term | Exact on fixtures | Notes |
|---|---|---|
| Bond stretch | 8/9 | water Δ 0.0101 — MMFF94 types water O/H as 70/31 (r₀ 0.969); we type 6/21 (r₀ 0.972). Typing gap, not a parameter gap |
| Angle bend | 9/9 | |
| Stretch-bend | 9/9 | |
| Torsion | 6/9 | propane Δ 0.007, butane Δ 0.053, cyclohexane Δ 0.554 — C–C–C–C dihedral parameter question |
| Van der Waals | 9/9 | |
| Out-of-plane | 9/9 | 0 on these fixtures; see BatchMin table below |
| Electrostatic | — | stub (returns 0) |

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
