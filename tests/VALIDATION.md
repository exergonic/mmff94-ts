# Validation

## Purpose

This document describes the validation of the mmff94-ts library.
It states what matches the reference and what does not.
The Implementer's Notes document has the detailed technical records.

## Reference Data

We validated the library against Halgren's MMFF94 validation suite.
The suite has 753 molecules.
It gives the structure, the total energy, the seven energy
components, and the partial charges for each molecule.
It also gives the atom types from the original program.

The molecule-by-molecule comparison is committed as generated
documentation (`npm run validation:doc` regenerates it from the suite
files via `tests/scripts/generate-validation-doc.ts`):

- `docs/validation/total-energies.txt` — all 753 totals, ours side by
  side with the suite's OPTIMOL and BatchMin columns.
- `docs/validation/per-term-and-charges.txt` — every per-term energy
  delta and per-atom charge delta, molecule by molecule.

The claims below are the summary of those files; if the two ever
disagree, the generated files and the generator are the truth.

## Tolerances

We use two tolerances:

- Energies: 0.0001 kcal/mol per term.
- Partial charges: 0.001 e per atom.

## Atom Types

The library assigns the MMFF94 atom type to every atom of every
molecule.
We compared the assignments with two independent references:

1. OpenBabel's canonical types. All 753 molecules match exactly.
2. The original program's own assignments. 749 molecules match
   exactly. The remaining four atoms have different type labels,
   but their parameters are identical. The energy checks prove
   this.

## Energies

The library computes all seven MMFF94 energy terms: bond stretch,
angle bend, stretch-bend, torsion, out-of-plane bend, van der Waals,
and electrostatic.

We compared each term with the reference components.
The reference components are single-point calculations at the stored
geometries.
Therefore, a difference on a correctly typed molecule is a bug in
that term.

The table gives the results at the 0.0001 kcal/mol level.

| Term | Molecules that match |
|---|---|
| Bond stretch | 753 of 753 |
| Angle bend | 753 of 753 |
| Stretch-bend | 753 of 753 |
| Torsion | 753 of 753 |
| Out-of-plane bend | 753 of 753 |
| Van der Waals | 753 of 753 |
| Electrostatic | 751 of 753 |

All comparable terms match within the 0.0001 kcal/mol tolerance for
all 753 molecules. The two electrostatics exclusions (AN11A, DOZNIP)
are the delocalized-anion reference anomalies below.

The suite has one bond with parameters from the MMFF94 empirical
rules (part V of the papers).
This bond is the hydroxide O–H bond of OHMW1.
The library generates those parameters with the published rules.
The bond matches to 0.0000014 kcal/mol.

We also compared the fixture molecules with OpenBabel's energy
output.
All seven terms and the total match to five decimal places.

## Partial Charges

The partial charges come from the bond charge increment (BCI) model.
They match the reference values to 0.001 e on 749 molecules.

## Gradients

We verified every analytical gradient with a finite-difference
calculation.
We perturbed each coordinate of each atom of each reference molecule
by 0.000001 Å.
The worst relative error is 8.5e-8.

## Outliers

Two molecules have one energy term that we cannot reproduce: for
each, the reference itself is inconsistent for that term. We verified
all other terms of these molecules against two independent
implementations: Tinker and OpenBabel — both agree with our values.
Two further molecules (JALSOE, SO18A) have all seven energy terms
reproduced; only their reference partial charges are not comparable
(the reference adjusts them to the dative representation).

| Molecules | Term | Reason |
|---|---|---|
| AN11A, DOZNIP | Electrostatic | The anionic five-ring nitrogen has no uniform primary charge (Halgren states this). Each implementation gives a different value. Tinker does not give the term for AN11A. |
| JALSOE, SO18A | Partial charges | The reference adjusts the sulfur-sulfur bonds to the dative representation. The reference charges are therefore not comparable. All seven energy terms match. |

The former FE2PW3/CU1PW1 van der Waals split is closed: the +2/+1
metal-hydrate cations carry their own vdW rows (which differ from the
+3/+2 rows only in the polarizability). OpenBabel's canonical typing
collapses them onto the +3/+2 classes; the vdW term now bridges to
the +2/+1 rows by formal charge (see implementer-notes §5.1), and
both molecules rejoin the census on all terms.
