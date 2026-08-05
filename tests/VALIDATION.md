# Validation

Facts: what `mmff94-ts` has been validated against, and what remains
undone. All numbers below are produced by `npm run test`
(regenerate the suite scoreboard with
`npx tsx tests/scripts/energy-scoreboard.ts`).

Reference data: Halgren's MMFF94 validation suite (753 structures,
per-component energies in `MMFF94_bmin.log`, per-atom reference
partial charges in `MMFF94.mmd`, canonical OpenBabel types in
`mmff94-atom-types.json`) lives in `tests/fixtures/validation-suite/`.
Tolerances: energies within 0.05 kcal/mol per term, charges within
10⁻³ e⁻ per atom.

## Validated

**Atom typing — 753/753 vs OpenBabel (100%).**
Every suite molecule reproduces OpenBabel's canonical MMFF94 atom
types exactly (`atom-types-suite.test.ts`). The 550 pre-recovery
exact codes are pinned as `KNOWN_GOOD` regression guards
(`tests/known-good.ts`); zero regressions since the pin.

**Energies — 749/749, all seven terms vs BatchMin.**
Every typing-exact molecule with a BatchMin per-component reference
matches all seven terms within 0.05 kcal/mol
(`energy-scoreboard.ts`, asserted per-term in
`validate-against-suite.test.ts`):

| term | exact | max \|Δ\| | worst |
|---|---|---|---|
| bond stretch | 749/749 | 0.0005 | OHMW1 |
| angle bend | 749/749 | 0.0007 | GAKTAN |
| stretch-bend | 749/749 | 0.0000 | — |
| torsion | 749/749 | 0.0001 | TAJSUS |
| out-of-plane | 749/749 | 0.0163 | KESNEB |
| van der Waals | 749/749 | 0.0000 | — |
| electrostatic | 749/749 | 0.0035 | SOHXOC |

BatchMin's log is a single-point calculation at the `.mmd` geometry
(per the suite README), so a delta on a typing-exact molecule is a
term or lookup bug by construction.

**Partial charges — 749/749 vs the suite's reference values.**
Per-atom |Δ| < 10⁻³ e⁻ on the same reproducible set
(`charges-suite.test.ts`).

**obenergy cross-check — 16 molecules, exact.**
All seven terms, the total, and the per-atom partial charges match
the obenergy logs to 5 decimals (4 for charges)
(`reference-comparison.test.ts`).

**Out-of-plane pins vs BatchMin (|Δ| < 0.05).**

| code | ours | BatchMin | Δ |
|---|---|---|---|
| DADDAN | 0.255548 | 0.255547 | +0.000000 |
| GIDJUY | 0.216936 | 0.216938 | −0.000002 |
| VEJWOW | 0.177150 | 0.177154 | −0.000004 |
| DIKGAF | 0.160155 | 0.158925 | +0.001230 |
| FAXVAB | 0.127921 | 0.126658 | +0.001263 |
| GEXGIZ | 0.123819 | 0.123820 | −0.000001 |
| VIRBON | 0.101801 | 0.102969 | −0.001167 |
| AMHTAR01 | 0.224483 | 0.224486 | −0.000003 |

**Gradients.** Every analytical gradient is finite-difference verified
on every atom of every reference molecule (worst relative error
8.5×10⁻⁸).

**Tests: 187 passing, 0 skipped.**

## Undone

| codes | terms | why the reference is not reproduced |
|---|---|---|
| AN11A, DOZNIP | electrostatic, charges | Anionic 5-ring N⁻ (type 76): no uniform primary charge q⁰(76) reproduces the reference charges from part V eq. (15) |
| FE2PW3, CU1PW1 | van der Waals | BatchMin's hydrated-metal vdW predates the X94 metal parameters; OpenBabel matches this transcription exactly (55.84481 vs 55.8448) and Tinker agrees, so the reference is the outlier |
| JALSOE, SO18A | charges only | BatchMin adjusts their S–S bonds to the "MMFF dative representation" (its log states this); all seven energy terms match |
