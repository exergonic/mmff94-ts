# Tinker 26.2 MMFF94: spurious −1 e net charge on neutral sulfones

Dear Dr. Ponder and Tinker developers,

We maintain an independent TypeScript implementation of MMFF94
([mmff94-ts](https://github.com/exergonic/mmff94-ts)), validated against
Halgren's original 761-molecule validation suite (per-term energies ≤1e-4
vs BatchMin, atom typing 761/761 identical to OpenBabel, BCI charges
<1e-3 e vs the published reference charges). While benchmarking our
geometry optimizer against `minimize`, we found that Tinker 26.2's MMFF94
charge derivation appears to assign a spurious **net molecular charge of
−1.0 e** to simple, formally **neutral sulfones**.

## Minimal reproduction

`dimethylsulfone` (CS(=O)(=O)C), 10 atoms, SMILES `CS(=O)(=O)C`, no
formal charges anywhere:

```
$ obabel -:"CS(=O)(=O)C" -osdf --gen3d best > dimethylsulfone.sdf

# Assign MMFF94 atom types with any correct typer. The relevant types:
#   S hexavalent = class 18 -> Tinker type 64 ("SO2")
#   sulfonyl O   = class 32 -> Tinker type 107 ("O2CM" — Halgren's
#                 published typing for sulfone oxygens)
#   C sp3 = class 1 (type 1), H = class 3 (type 23)
```

`dimethylsulfone.txyz` (Tinker canonical numbering; S = type 64 "SO2",
the two sulfonyl O = type 107, which is MMFF94 class 32 / `O2CM` —
Halgren's published typing for sulfone oxygens):

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

 Total Electric Charge :                 -1.00000 Electrons
```

The molecule is neutral by construction; both the MMFF94 bond-charge-
increment model as we implement it (and RDKit's independent MMFF94
implementation, which lands at the same minima and energies on these
molecules) derive a total charge of 0.000 for it.

## Scope

Reproduces on every small sulfone/sulfonyl we tried:

| molecule | true charge | Tinker `analyze [M]` net charge |
|---|---|---|
| dimethylsulfone | 0 | **−1.00** |
| methanesulfonamide | 0 | **−1.00** |
| pyridine, imidazole, indole, aniline, urea, N-methylacetamide, acetophenone | 0 | 0 ✓ |

In a larger benchmark (29 drug-like molecules from the OpenFF Industry
set), every molecule containing a `S(=O)2` group showed the same −1 e
phantom; the effect is large in practice — minimizing from an identical
starting geometry with identical types, Tinker descends into an
artificial ion-stabilized basin up to **56 kcal/mol deep** relative to
our implementation and RDKit, purely because the electrostatics differ.

## Guessing at the cause (from the outside)

MMFF94 reuses atom type 32 (`O2CM`, "carboxylate anion O") for the
terminal oxygens of sulfones/sulfonyl groups — that is per Halgren's
published type table, and it is what both RDKit and this library
assign. In a carboxylate, the −1 formal-charge contribution carried by
the two type-32 oxygens is balanced by +1 on the carboxyl carbon. Our
suspicion is that Tinker's MMFF94 charge assignment applies that
carboxylate formal-charge rule to type 32 without the compensating
q⁰ for the sulfone/sulfonyl family, so every S(VI)(=O)₂ group ends up
with a net −1 e. The `mmff94.prm` BCI table itself looks fine to us;
we believe the issue sits in the formal-charge/q⁰ logic that runs
before BCI summation when the MMFF94 force field is active.

## What we'd appreciate

- A confirmation whether this reproduces on your end with the files
  above (Tinker 26.2, Linux x86-64 build).
- If confirmed, a pointer to the responsible routine so we can cite it
  precisely, and any known workaround (explicit `CHARGE` directives do
  not override the internally derived values in our testing).

Happy to send more failing molecules or a fuller dataset. Everything
above was produced with the stock `mmff94.prm` shipped with Tinker 26.2
(April 2026), no custom parameters.

— Billy Wayne McCann ([@exergonic](https://github.com/exergonic))
  mmff94-ts: https://github.com/exergonic/mmff94-ts
