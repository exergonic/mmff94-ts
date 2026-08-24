# MMFF94: spurious −1 e net molecular charge on neutral sulfones (`kchargem` hard-codes −0.5 on type 107)

## Environment

- Tinker 26.2 (April 2026), stock `params/mmff94.prm`, no custom parameters
- Linux x86_64
- Observed through both `analyze` (net-charge report) and `minimize` (energies)

## Summary

MMFF94 reuses atom type 107 / class 32 (`O2CM`) for the terminal oxygens of
sulfones and sulfonyl groups — that is per Halgren's published type table,
and it is what both RDKit and our independent implementation assign.
`kchargem` then gives every type-107 atom a hard-coded base partial charge
of **−0.5 e** (`source/kcharge.f`):

```fortran
if (it .eq. 107)  pchg(i) = -0.5d0
```

Type 64 (`SO2` sulfur) receives no compensating base charge, so any
molecule containing an S(VI)(=O)₂ group comes out with a **net molecular
charge of −1.0 e**, even when it is formally neutral.

## Minimal reproduction

Dimethylsulfone — 11 atoms, SMILES `CS(=O)(=O)C`, no formal charges anywhere.

```bash
$ obabel -:CS(=O)(=O)C -osdf --gen3d best > dimethylsulfone.sdf
```

Assign MMFF94 atom types with any correct typer. Relevant types: S hexavalent
→ 64 (`SO2`); sulfonyl O → 107 (`O2CM`); C sp3 → 1; H → 23.

`dimethylsulfone.txyz`:

```
11 dimethylsulfone
    1  C       0.985800      0.021000     -0.050300     1  2 6 7 8
    2  S       2.758800     -0.067800     -0.011300    64  1 3 4 5
    3  O       3.171300     -0.098000      1.378500   107  2
    4  O       3.171300     -1.105700     -0.935600   107  2
    5  C       3.258600      1.495000     -0.691900     1  2 9 10 11
    6  H       0.590500     -0.912200      0.356300    23  1
    7  H       0.660600      0.145100     -1.084800    23  1
    8  H       0.660300      0.863100      0.563000    23  1
    9  H       2.869900      2.296000     -0.060700    23  5
   10  H       2.870100      1.578300     -1.708500    23  5
   11  H       4.350100      1.526100     -0.705300    23  5
```

```bash
$ printf "parameters /path/to/mmff94.prm\nMMFF-PIBOND\n" > dimethylsulfone.key
$ echo M | analyze dimethylsulfone
```

### Expected

`Total Electric Charge : 0.00000 Electrons` — the molecule is neutral, and
per Halgren part V eq. (15) the q⁰ of a type-32/107 oxygen is −(n−k)/n with
n = number of terminal O's on the center and k = number of oxo oxygens in
the *neutral parent oxyacid*. For tetracoordinate S(VI), k = 2, so q⁰ = 0;
the strong S–O polarization lives in the BCI term instead.

### Actual

```
Total Electric Charge :                 -1.00000 Electrons
```

## Scope

Every sulfone/sulfonyl compound we tried reproduces it; controls without an
S(VI)(=O)₂ group are clean:

| molecule | true charge | Tinker net charge |
|---|---|---|
| dimethylsulfone | 0 | **−1.00** |
| methanesulfonamide | 0 | **−1.00** |
| pyridine | 0 | 0 ✓ |
| imidazole | 0 | 0 ✓ |
| indole | 0 | 0 ✓ |
| aniline | 0 | 0 ✓ |
| N-methylacetamide | 0 | 0 ✓ |
| urea | 0 | 0 ✓ |
| acetophenone | 0 | 0 ✓ |

On a 29-molecule drug-like benchmark (OpenFF Industry set, QM-start
geometries), every molecule containing an S(VI)(=O)₂ group showed the same
phantom −1 e. The practical effect is large: minimizing from an identical
starting geometry with identical atom types, the artificial ion-stabilized
surface leads `minimize` to basins up to **56 kcal/mol** below what RDKit's
independent MMFF94 (and ours) reach for the same molecules. Sulfoxides (one
terminal O) escape with −0.5; we would expect phosphine oxides and nitro
groups to be mildly affected through the same type-reuse path.

## Cross-validation

This is not a disagreement between two independent implementations about
force-field conventions:

- On the 11 benchmark molecules whose chemistry does **not** trigger the
  issue, `minimize` started *from our converged geometry* returns the final
  energy identical to ours **to four decimal places** — same force field,
  same answer.
- On the affected molecules, RDKit's independent MMFF94 derives net 0 and
  lands where we do (e.g. mol with SMILES containing a sulfonamide: ours
  −27.51, RDKit −27.56, Tinker −83.75).
- Our charge model reproduces BatchMin reference charges across Halgren's
  761-molecule validation suite to <1e-3 e per atom, including all charged
  species in the suite.

## Root cause

`source/kcharge.f`, subroutine `kchargem` (26.2 source, ~line 230): the
−0.5 base charge is applied unconditionally per type-107 atom. In a
carboxylate this is balanced by +1 on the carboxyl carbon via its own
q⁰/base contribution, but for the sulfone/sulfonyl/sulfonyl-imine uses of
type 107 there is no such compensation anywhere, and no (n, k) environment
logic distinguishes the two cases.

## Suggested fix

Compute q⁰ for type 107 (and its sibling reused types) by the eq. (15)
environment rule — −(n−k)/n with the center-specific k — instead of a flat
per-type constant. That reproduces carboxylate −½, sulfate −½, phosphate
−¾ …, and 0 for neutral sulfones/sulfoxides/nitro/phosphine-oxide groups,
matching the published spec values. Alternatively, a compensating base
charge on type 64 would patch the sulfone case specifically.

Happy to provide more failing molecules or fuller datasets if useful.
