# MMFF94 Validation Suite — per-molecule SDF files

One SDF per suite molecule (`<CODE>.sdf`), generated from the suite's
own `MMFF94.mmd` (the BatchMin **hypervalent** representation — the
same input the library's pipeline and the BatchMin reference energies
use). Open any file in Avogadro, PyMOL, or any SDF reader.

Each file carries:

- the MOL V2000 block: atoms, bonds with orders, and the formal
  charges from `MMFF94.fc_hypervalent` (as `M  CHG` lines);
- two data fields in atom order:
  - `MMFF94_atom_types` — the MMFF94 atom types assigned by
    `assign_atom_types` (the library's own typing);
  - `MMFF94_reference_charges` — the reference partial charges (the
    `.mmd` pchg column, what BatchMin used).

Caveats:

- The hypervalent representation shows hexavalent sulfur and
  pentavalent phosphorus (two double bonds from formally neutral
  oxygens), as in the reference — the suite's `MMFF94_dative.mol2`
  is the alternative dative representation.
- Atom types are class numbers (the OpenBabel/Tinker/OPTIMOL
  numbering), not the original 1-178 MMFF94 type numbers.
- Titles come from `MMFF94.titles`; the eight ERULE fragments have
  no names there and carry just their code.

Regenerate with:

```
npx tsx tests/scripts/export-suite-sdf.ts
```

The script round-trips every file through the library's own SDF
parser (atoms/bonds must match the mmd exactly) and reports any
mismatch.
