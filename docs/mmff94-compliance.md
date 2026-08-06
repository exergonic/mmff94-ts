# MMFF94 Compliance Statement

**mmff94-ts 0.1.0-alpha.1** — a statement of what this library implements,
how it was validated, and what it deliberately does not do. Modeled on
Wavefun's MMFF94 compliance statement
(https://downloads.wavefun.com/FAQ/MMFF94_compliance.html).

## Which MMFF94 variant

**MMFF94 proper** — the original force field as defined in Halgren's
series (J. Comput. Chem. 1996, 17, 490-641; 1999, 20, 720-748). The
MMFF94s "steric" variant (part VI) is **not** implemented: MMFF94s
replaces the negative amide-N out-of-plane constants with small
positive ones and otherwise alters the balance of terms, and this
library deliberately implements the standard variant only.

## Energy terms — all seven, with the published functional forms

| Term | Functional form | Reference |
|---|---|---|
| Bond stretch | 143.9325·(k_b/2)·Δr²·[1 + cs·Δr + 7/12·cs²·Δr²] | eq. (2), part I |
| Angle bend | 0.043844·(k_a/2)·Δθ²·(1 + cb·Δθ) (eq. 4 for linear centers) | eqs. (3)-(4), part I |
| Stretch-bend | 2.51210·[k_sb_IJK·(r₁−r₁₀) + k_sb_KJI·(r₂−r₂₀)]·(θ−θ₀) | eq. (5), part I |
| Torsion | Σ (V_n/2)·[1 + cos(n·τ − γ_n)], n = 1, 2, 3 | eq. (7), part I |
| Van der Waals | ε_ij·[(1.07·R*/(r+0.07·R*))⁷·(1.12·R*⁷/(r⁷+0.12·R*⁷) − 2)] | eq. (8), part I |
| Electrostatic | 332.0716·q_i·q_j/(r + 0.05), 1-4 pairs ×0.75 | eq. (6), part III |
| Out-of-plane | 0.043844·(k_oop/2)·χ² (Wilson angle, any tri-coordinate center) | eq. (6), part I |

The conversion factors are the authoritative Halgren values
(143.9325, 0.043844, 2.51210, 332.0716). OpenBabel's `obenergy`
output uses MM2-style half-factors internally; this library does not.

1-4 electrostatics are scaled ×0.75 and 1-4 vdW is **not** scaled,
per Halgren's explicit statement (part III p. 496).

## Parameter provenance

All numeric parameters come from OpenBabel's text-format `.par` files
(`mmffbond.par`, `mmffang.par`, `mmffstbn.par`, `mmfftor.par`,
`mmffvdw.par`, `mmffchg.par`, `mmffpbci.par`, `mmffoop.par`,
`mmffdef.par`, `mmffprop.par`), which are mechanical transcriptions of
Halgren's published tables. The build-time extraction script
(`scripts/extract-mmff94-par.py`) converts format only; the generated
tables are committed so the library builds without Python.

The part V empirical-rule generation (the designed fallback when the
class-scoped lookup misses) is implemented from the published rules:
bond eqs. (18)-(19) + Table V, the θ₀ protocol + eq. (20) + Table VI,
the torsion rules (pp. 631-632) + Table X, and the eq. (17) BCI
fallback. All four families are pinned by hand-computed unit tests
(41 tests) and cross-checked against Tinker's independent
transcription (kbond.f, kangle.f, ktors.f — local build).

## Cross-check results

Validated against Halgren's own 753-molecule MMFF94 Validation Suite
(the definitive reference — the original program's output):

- **Atom typing: 753/753 type-exact** vs OpenBabel's canonical types,
  cross-checked against the original program's own per-atom
  assignments (the suite's OPTIMOL log): 749/753 identical; the four
  remaining atoms — the two metal-hydrate cations' oxidation-state
  types (the original program's FE+2/CU+1 vs OpenBabel's FE+3/CU+2),
  bridged at the vdW lookup by formal charge (the parameter rows
  differ only in the polarizability; see below), and the two dative
  sulfone-O atoms, parameter-identical — as the three-way energy
  checks prove. Tinker's prm atom table shares the same class
  numbering, so the parameter-level cross-checks exercise the same
  types.
- **Energies: every BatchMin per-component reference reproduced to
  ≤10⁻⁴ kcal/mol — six of seven terms on all 753 molecules,
  electrostatics on 751** (worst residual 6.8×10⁻⁵; the two per-term
  reference anomalies are itemized below and were cross-verified
  three ways against Tinker's independent transcription on
  2026-08-05). The former FE2PW3/CU1PW1 van der Waals split is
  closed: those cations carry their own +2/+1 vdW rows (differing
  from the +3/+2 rows only in the polarizability), and the vdW term
  now selects them by the formal charge — both molecules rejoin the
  census on all seven terms. JALSOE/SO18A's full seven terms are in
  the census. The suite's only empirical-rule bond (OHMW1's
  hydroxide O–H) matches to 1.4×10⁻⁶; Tinker independently
  reproduces it (0.7654 vs 0.765397246118).
- **Partial charges: 749/749 to <10⁻³ e⁻** per atom (BCI model,
  eqs. (14)-(17)).
- The molecule-by-molecule evidence is committed in
  [`docs/validation/`](validation/): `total-energies.txt` lists all 753
  totals side by side with the suite's own OPTIMOL and BatchMin values
  (751/753 to 10⁻³ kcal/mol; the two exceptions are the reference
  anomalies below), and `per-term-and-charges.txt` lists every
  per-term energy delta and per-atom charge delta. Regenerate with
  `npm run validation:doc`.
- Gradients are checked against finite differences (δ = 10⁻⁶ Å,
  relative error < 10⁻⁵; worst observed 8×10⁻⁸).
- A second independent check against OpenChemLib (devDependency)
  reproduces total energies on the fixtures.

The documented reference anomalies (the two type-76-anion
electrostatics, cross-checked three ways against Tinker — see
`tests/VALIDATION.md`) are per-term: every other component of those
two molecules is reproduced. The JALSOE/SO18A dative molecules'
energies are fully verified (their reference partial charges are
dative-adjusted and not comparable).

## Browser and build

- Zero runtime dependencies; the library compiles to plain ESM in
  `dist/` (NodeNext emit — every import carries its `.js` extension).
- Smoke-tested in headless Chromium: parse → type → energy →
  gradient run against the built dist with results identical to Node
  (E = −4.7344 kcal/mol on ethane).

## Limitations (known and deliberate)

- MMFF94s (part VI) is not implemented; MMFF94 only.
- The part V eq. (18) bond rule is implemented in its **plain measured
  form** (no δ = 0.008, no mltb/BOij corrections): the reference's own
  output matches the plain form to 1.4×10⁻⁶, and Tinker's kbond.f
  implements the same form. OpenBabel's transcription (δ-inclusive)
  is the documented outlier for that one bond.
- The Badger's-rule fallback beyond Table V bond pairs is not
  implemented (the suite never needs it).
- The optimizer (L-BFGS) converges 16/16 fixtures from both starts at
  max|g| < 0.05 but has documented stalls on some geometries
  (see `tests/VALIDATION.md` — a flat-lands line-search pathology);
  steepest descent is the robust fallback.
- The dielectric is the in-vacuo default D = 1.0 (the D = r solvent
  alternative is not exposed).
- Per-interaction energy breakdowns (each bond/angle/torsion/vdW pair
  as queryable data) are a documented stretch goal, not yet shipped.
