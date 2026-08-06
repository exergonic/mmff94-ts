# MMFF94 Compliance Statement

**mmff94-ts 0.1.0-alpha.1** — a statement of what this library
implements, how it was validated, and where it differs from the
original MMFF94 program. The format follows Wavefun's MMFF94
compliance statement
(https://downloads.wavefun.com/FAQ/MMFF94_compliance.html).

## Energy terms

MMFF94 has seven energy terms. Each term uses the functional form
published by Halgren:

| Term | Functional form | Reference |
|---|---|---|
| Bond stretch | 143.9325·(k_b/2)·Δr²·[1 + cs·Δr + 7/12·cs²·Δr²] | eq. (2), part I |
| Angle bend | 0.043844·(k_a/2)·Δθ²·(1 + cb·Δθ) (eq. 4 for linear centers) | eqs. (3)-(4), part I |
| Stretch-bend | 2.51210·[k_sb_IJK·(r₁−r₁₀) + k_sb_KJI·(r₂−r₂₀)]·(θ−θ₀) | eq. (5), part I |
| Torsion | Σ (V_n/2)·[1 + cos(n·τ − γ_n)], n = 1, 2, 3 | eq. (7), part I |
| Van der Waals | ε_ij·[(1.07·R*/(r+0.07·R*))⁷·(1.12·R*⁷/(r⁷+0.12·R*⁷) − 2)] | eq. (8), part I |
| Electrostatic | 332.0716·q_i·q_j/(r + 0.05), 1-4 pairs ×0.75 | eq. (6), part III |
| Out-of-plane | 0.043844·(k_oop/2)·χ² (Wilson angle, any tri-coordinate center) | eq. (6), part I |

The conversion factors are the published Halgren values (143.9325,
0.043844, 2.51210, 332.0716). 

Pairs separated by three bonds (1-4 pairs) have their electrostatic
interaction scaled by 0.75. Their van der Waals interaction is not
scaled. Halgren states this explicitly (part III, p. 496).

## Parameters

The numeric parameters come from Halgren's published tables.

The papers also define rules to generate parameters when a lookup
misses (part V): the bond rules (eqs. 18-19, Table V), the angle
rules (eq. 20, Table VI), the torsion rules (pp. 631-632, Table X),
and the charge fallback (eq. 17). This library implements these
rules as published. Each rule is pinned by hand-computed unit tests
and cross-checked against Tinker's independent implementation.

## Validation

The library was validated against three references: the original
program's own 753-molecule validation suite (BatchMin and OPTIMOL
output), OpenBabel, and Tinker.

**Atom types.** The types match OpenBabel's canonical types on all
753 molecules. They match the original program's own types on 749
of 753 atoms. The four remaining atoms are described under
"Differences".

**Energy terms.** Six of the seven terms match the BatchMin
reference on all 753 molecules to 0.0001 kcal/mol or better. The
electrostatics match on 751 molecules. The worst residual is
6.8e-5 kcal/mol. The two excluded molecules are described under
"Differences".

**Total energies.** The totals match the BatchMin reference to
0.001 kcal/mol on 751 of 753 molecules.

**Partial charges.** The charges match the reference values to
0.001 e per atom on 749 molecules.

**Gradients.** Every analytical gradient is checked against a
central finite difference (step 1e-6 Å). The worst relative error
is 8e-8.

The molecule-by-molecule evidence is committed with the library:
the 753 totals side by side with the reference values, and the
per-term and per-atom-charge deltas for every molecule
(`docs/validation/`).

## Differences from the original

The library differs from the original program in the following
places. Each difference is deliberate and documented.

**The two delocalized-anion electrostatics (AN11A, DOZNIP).** The
anionic five-ring nitrogen (type 76) has no uniform primary charge.
Halgren states this in the papers. The three implementations each
produce a different value. The reference itself uses different
charges in different environments: on the symmetric anion (JILWUW)
our charge matches the reference exactly; on the two asymmetric
anions it does not. No single value can match all three
implementations. The electrostatics of the two molecules are
therefore not compared. Every other term of these molecules
matches.

**The metal-hydrate cations (FE2PW3, CU1PW1).** The original
program types these cations with the +2/+1 oxidation states and
their own van der Waals rows. OpenBabel's canonical typing assigns
them the +3/+2 classes. The two row sets differ only in the
polarizability. The van der Waals term selects the +2/+1 rows by
the formal charge. With this bridge, both molecules match the
reference to 1e-6 kcal/mol on all seven terms.

**The angle cubic constant.** The paper gives the cubic constant as
−0.007 per degree, "or, more precisely, −0.4 per radian". This
library, BatchMin, and Tinker use the precise value. OpenBabel uses
the rounded value. The difference grows with the cube of the angle
deviation. It shows only on strained angles: three small molecules
differ from OpenBabel by up to 0.0007 kcal/mol in the angle term.

**The empirical bond rule.** The bond-length rule of part V (eq.
18) is implemented in its plain published form, without the delta
correction. The reference and Tinker do the same. OpenBabel's
transcription includes the correction for that one bond. The
reference bond matches to 1.4e-6 kcal/mol.

**The JALSOE/SO18A reference charges.** The reference adjusts these
molecules to the dative representation. Their reference charges are
therefore not comparable to any BCI model. Their energies match on
all seven terms.

## Limitations

- The part V Badger's-rule fallback for bonds outside Table V is
  implemented: the length from eq. (18) and the force constant from
  the Herschbach-Laurie parameterization of Badger's rule
  (k = 1.86/(r − d)³), with the d values derived from the E94 rows of
  Table V. The validation suite never needs it, so it is pinned by
  unit tests rather than reference energies.
- The dielectric is the in-vacuo value D = 1.0. The alternative
  solvent model (D = r) is not exposed.
