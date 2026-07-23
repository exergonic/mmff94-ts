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
| Van der Waals (buffered 14-7) | ⬜ |
| Electrostatic (BCI model) | ⬜ |
| Out-of-plane bending | ⬜ |
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

Every function is pure — no hidden state, no global setup, no
`new ForceField()`. Data flows through a pipeline:

```
SDF → parse_sdf() → Molecule → assign_atom_types() → TypedMolecule
                                                         ↓
                                                    calc_energy()
                                                         ↓
                                                  EnergyComponents
```

## Validation

Tested against **OpenBabel 3.1.1** (`obenergy -ff MMFF94`) and
**Dr. Halgren's MMFF94 Validation Suite** (753 structures from the
Cambridge Structural Database + small molecules/ions) which lives at
`tests/fixtures/validation-suite/`.

## Installation

```bash
npm install mmff94-ts
# or
bun add mmff94-ts
```

## Development

```bash
# Install dependencies
npm install

# Run tests
npm test

# Regenerate parameter tables from OpenBabel .par files
python scripts/extract-mmff94-par.py
```

## License

MIT. Parameter data extracted from OpenBabel `.par` files which
contain only numerical facts from Halgren's published papers.
