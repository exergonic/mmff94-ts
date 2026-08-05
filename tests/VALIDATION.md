# Validation

What `mmff94-ts` has been validated against, and what remains
undone. All numbers below are produced by `npm run test`
(regenerate the suite scoreboard with `npx tsx tests/scripts/energy-scoreboard.ts`).

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
| bond stretch | 747/747 | 0.0001 | DEWJEU |
| angle bend | 747/747 | 0.0000 | — |
| stretch-bend | 747/747 | 0.0000 | — |
| torsion | 747/747 | 0.0000 | — |
| out-of-plane | 747/747 | 0.0000 | — |
| van der Waals | 747/747 | 0.0000 | — |
| electrostatic | 747/747 | 0.0001 | DONFOB |

(The "exact" columns count molecules within 10⁻⁴ kcal/mol of the
BatchMin component — the per-term residual census
(`tests/scripts/residual-census.ts`) — not the 0.05 gate. At that
bar ALL SEVEN TERMS are machine-exact on every molecule (747/747).
The last stretch residual — OHMW1's hydroxide O–H, the only bond in
the suite resolved by MMFF94's empirical-rule parameter generation
(part V; the CCL `MMFF94.empirical_rule_parameters` file lists
"Empirical rule bond parameters: 0 4 5" — one of only three
original-suite molecules exercising any empirical rule, with
CEWYIM30/KEPKIZ angles) — closed 2026-08-05 at 1.4×10⁻⁶. No stored
row exists in any numbering; BatchMin generated k_b/r₀ at runtime
from eqs. (18)-(19) of part V. mmff94-ts implements that generation
(`parameters/empirical.ts`: eq. 18 Schomaker-Stevenson/Blom-Haaland
r₀, eq. 19 inverse-sixth force-constant rule with Table V). The
implemented form is MEASURED against the reference, not the paper's
literal eq. (18): the paper adds a δ = 0.008 Å shrinkage and the
mltb/BOij radius reductions, but the reference's O–H matches the
plain form — r₀ = 0.72 + 0.33 − 0.050·1.30^1.4 = 0.9778 — to 1.4e-6,
and Tinker's independent transcription (kbond.f: bl = rad0a + rad0b
− cst·|χa−χb|^1.4, no δ, no corrections) implements exactly that
form and reproduces the reference's Bond Stretching (0.7654 vs
0.765397246118) on its own build (2026-08-05). The published tables
are used as-is (r(O) = 0.72 — the CCL errata and Tinker's
mmffcovrad agree). BatchMin's log flags the bond: "NO CHARGE
INCREMENT CAN BE FOUND FOR THE BOND BETWEEN ATOMS O #4 AND H #5" and
"low quality stretch parameters = 1"; the BCI side needs no
empirical code — charges.ts's unparametrized fallback
w = pbci(I) − pbci(K) IS the paper's eq. (17) (the refined rule).
TAJSUS's torsion
(1.2×10⁻⁴) closed 2026-08-04: the class-2 torsion branch must not
fire when the central bond is an aromatic ring bond — type 80 (CIM+)
lacks the arom par flag, so is_aromatic_bond now treats it as
aromatic (mirroring bond_type_flag); the triazole ring bonds then
read class 0 and the reference's V₂ = 4.0 rows resolve.

BatchMin's log is a single-point calculation at the `.mmd` geometry
(per the suite README), so a delta on a typing-exact molecule is a
term or lookup bug by construction. The unit conversions use the
exact factors (143.9325·(π/180)² = 0.0438443467… for angle/oop, not
the rounded 0.043844): the rounded form left all 168 bend components
~1e-5 relative low.

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
