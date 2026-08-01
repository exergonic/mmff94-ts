# Validation Ledger

Living record of what `mmff94-ts` has been validated against —
molecule by molecule, term by term. The numbers here are produced by
`npm run test`; update this file whenever the reference tests change.
The README's [Validation](../README.md#validation) section is the
condensed public version of this ledger.

## Method

- **Fixtures**: 9 molecules in `tests/fixtures/sdf/`, typed with our
  own `assign_atom_types()`.
- **Fixture references**: `obenergy -ff MMFF94` logs in
  `tests/references/` (OpenBabel 3.1.0, generated via
  `tests/scripts/get_mmff94_breakdown.py`).
- **Suite references**: Halgren's 753-molecule MMFF94 validation
  suite in `tests/fixtures/validation-suite/` — structures from
  `MMFF94.mmd`, BatchMin component energies from `MMFF94_bmin.log`.
- **Asserted in CI** (`tests/reference-comparison.test.ts`,
  `tests/validate-against-suite.test.ts`):
  - bond, angle, stretch-bend, torsion, vdW vs obenergy: all 9
    fixtures, |Δ| < 0.02
  - out-of-plane vs BatchMin: 8 suite molecules, |Δ| < 0.05

## Fixtures — per-term deltas vs obenergy (kcal/mol, |ours − ref|)

From the printed comparisons in `reference-comparison.test.ts`
(2026-08-01). "stub" in the elec column = we return 0; the reference
value is shown in parentheses when nonzero.

| molecule | bond | angle | strbnd | torsion | vdw | oop | elec |
|---|---|---|---|---|---|---|---|
| benzene | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | stub (ref 3.07810) |
| butane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | stub |
| cyclohexane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | stub |
| ethane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | stub |
| ethene | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | stub (ref 8.05300) |
| formaldehyde | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | stub |
| methane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | stub |
| propane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | stub |
| water | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | stub |

Total energies match exactly where every term matches: ethane
(−4.73436), butane (−5.07596), cyclohexane (−3.56091), propane
(−4.89729), methane (0.02638), formaldehyde (0.05416), water
(0.00000). The only remaining total deltas are the electrostatic
stub: benzene (3.078) and ethene (8.053).

## Suite — out-of-plane vs BatchMin (8 molecules, kcal/mol)

Asserted in `validate-against-suite.test.ts` (|Δ| < 0.05). These
molecules are chosen because our typing reproduces the reference
types exactly, so the comparison isolates the oop term.

| code | ours | BatchMin | Δ |
|---|---|---|---|
| DADDAN | 0.255548 | 0.255547 | +0.000000 |
| GIDJUY | 0.216936 | 0.216938 | −0.000002 |
| VEJWOW | 0.176902 | 0.177154 | −0.000252 |
| DIKGAF | 0.160155 | 0.158925 | +0.001230 |
| FAXVAB | 0.127921 | 0.126658 | +0.001263 |
| GEXGIZ | 0.122862 | 0.123820 | −0.000958 |
| VIRBON | 0.101801 | 0.102969 | −0.001167 |
| AMHTAR01 | 0.203026 | 0.224486 | −0.021460 |

## Known open questions

1. ~~Water typing (bond Δ 0.01008)~~ — **RESOLVED 2026-08-01.**
   MMFF94 has dedicated water types — O = 70 ("OXYGEN IN WATER"),
   H = 31 ("H-OH") — with bond `'0-31-70': r₀ 0.969` (already in our
   parameter table). We typed water as generic alcohol (6/21,
   r₀ 0.972). Fix: `assign_atom_types` now pre-scans for water (O with
   exactly two H neighbors, each H bonded only to that O → O=70,
   H=31). Water's bond term is now exact and its total matches
   obenergy exactly (0.00000). Note: the validation suite has no
   bare-water BatchMin reference (only hydrates), so obenergy is the
   only cross-check. The typing reference (`atom-types.test.ts`) pins
   water as [70, 31, 31]; the suite scoreboard is unaffected (no bare
   water in the 550 typed molecules).

2. ~~Cyclohexane torsion (Δ 0.554, 5% low)~~ — **RESOLVED 2026-08-01.**
   Root cause: `lookup_param`'s generic wildcard fallback ran before the
   reversed-direction exact match. An H–C–C–C dihedral (5,1,1,1) has no
   forward exact entry, so it matched the generic `'0-0-1-1-0'`
   (`*-1-1-*`, V3 = 0.300) instead of the exact reversed entry
   `'0-1-1-1-5'` (0.639 / −0.630 / 0.264) that OpenBabel uses (its
   direction canonicalization + exact lookup never reaches the
   wildcard). Fix: torsion.ts now tries exact keys in both directions
   before any wildcard (Halgren part I, p. 513 step-down protocol).
   All 9 fixture torsions now match obenergy exactly; the torsion
   assertion is extended to every fixture. The RING = AL log column is
   informational — the torsion class (TTijkl, part IV p. 609) is 0 for
   six-membered rings; class 5 applies only to five-membered rings.

3. **AMHTAR01 oop (Δ 0.021).** Largest suite oop delta, within
   tolerance. Ester/carboxyl pyramidalization; reference types the
   CO₂M carbon 41 and oxygens 32/32 where we type 3/7/6 — typing gap
   (see the atom-types suite scoreboard).

4. **AGLYSL01 components.** Stretch-bend 0.00000 vs 0.24423, torsion
   −3.00231 vs −4.71331, vdW 3.51572 vs 2.78652 — all traceable to
   typing (sulfonate S, ammonium N). The OOP row is exact (0/0).

5. **Electrostatic stub.** obenergy's benzene (+3.078) and ethene
   (+8.053) show what the stub leaves out; no fixture total with
   nonzero BCI charges can match until Phase 4.

6. **Typing scoreboard.** 65/550 suite molecules type-exact vs
   OpenBabel (11.8%; 165 ≥90% atoms; see `atom-types-suite.test.ts`).
   The torsion assertions on non-ethane fixtures stay skipped until typing
   isolates the term. Biggest remaining classes (from the worst-10):
   thioamide/thiazole S and ring C/N (44, 63–66, 78–81), sulfonyl/
   sulfate S (18), phosphate P (25) with H–S (71), amide/cyanamide N
   (10, 28, 42/43).

## Suite — per-component report (typing-exact molecules)

Regenerated by `npm run test` (validate-against-suite.test.ts "reports
per-component energies for typing-exact molecules"): per-term deltas vs
BatchMin on the currently typing-exact molecules (65 of 550, computed
fresh from `mmff94-atom-types.json`). Snapshot 2026-08-01 (after the
parameter-class work, commits `9997d0d`–`e75805c`):

| term | \|Δ\|≤0.05 | mean\|Δ\| | max\|Δ\| | worst |
|---|---|---|---|---|
| bond | 65/65 | 0.000 | 0.00 | — |
| angle | 65/65 | 0.001 | 0.02 | FAGVEO −0.02 |
| strbnd | 64/65 | 0.001 | 0.05 | JIYJAC +0.05 |
| torsion | 63/65 | 0.021 | 1.12 | FUVDOP +1.12, FILNOD +0.22 |
| vdw | 65/65 | 0.000 | 0.00 | — |
| oop | 65/65 | 0.000 | 0.00 | — |
| elec (stub) | — | 17.90 | 218.72 | magnitude the stub omits |

The comparison is valid: BatchMin's log is a single-point calculation at
the .mmd geometry (suite README), so deltas on typing-exact molecules
are term bugs by construction. What closed each thread (all verified
against OpenBabel ≡ BatchMin on the affected molecules):

- **FAGVEO +88.94 angle / JIYJAC +25.80 bond** — the parameter-class
  system was missing: BTij (part V p. 620, sbmb flag — case (a) applies
  to ANY sbmb pair, the "not aromatic" clause is not enforced) and
  ATijk select class 1/2 (BT sum), 3/5/6 (3-ring), 4/7/8 (4-ring).
  The old "priority sweep" found class entries only by accident — e.g.
  squarate C–C (3,3) hit class-1 `1-3-3` while butadiene's C2–C2 used
  the C=C double-bond params (`0-2-2`) for a conjugated single bond.
- **CUDNEU strbnd −11.35** — stretch-bend missed angles fall to the
  default-fsb table (mmffdfsb.par, element-row keyed); the Si–C bonds
  are stretched 0.27 Å so the Si angles dominate. The strbnd class is a
  REMAP of the angle class (GetStrBndType: 1↔2 by BT side, 2→3, 4→4,
  5→6/7, 8→11); remapped classes miss and use the fsb defaults.
- **KAGBOJ/METBZC10/FORGOI angles** — the angle step-down chain is the
  EqLvl3/4/5 equivalence levels (mmffdef.par; type 37's level 3 is 2),
  not generic wildcards; the empirical θ₀/ka rules (part II) are the
  final fallback; the lin flag selects the eq-4 cosine form.
- **DUYNOA/DIKWID/GEKXEZ vdW** — any pair with a donor atom uses the
  arithmetic-mean combination (ε halved and R*×0.8 only for
  donor-acceptor pairs); Waldman–Hagler only when neither is a donor.
- **FUVDOP −4.92 → +1.12** — the TTijkl torsion classes (part IV p.
  609: 1 = central BT, 2 = terminal BT, 4 = 4-ring, 5 = non-aromatic
  5-ring with an sp3 C) + the asymmetric EqLvl3/EqLvl5 chain in the
  order-canonical direction + double-bond-centered torsions now
  evaluated. The residual +1.12: the class-5 set matches OpenBabel's
  SSSR basis exactly (15/15 cycles verified) and OpenBabel ≡ BatchMin
  (28.5478 vs 28.5477), so the remainder is a BatchMin-internal detail
  of the N-cage torsion handling, not reachable from the par data.
  FILNOD +0.22 (benzothiazole S-oxide) is a separate small residual.

## Updating this ledger

1. Run `npm run test`; the reference-comparison tests print every
   fixture's per-term comparison, and the suite test prints AGLYSL01.
2. Update the fixture delta table and the oop table with the new
   numbers (re-run the oop table script if the values drift — the
   numbers live in the test's printed output).
3. When an open question is resolved, move it to a short "Resolved"
   note at the bottom with the commit that fixed it, and tighten the
   corresponding assertion in the test files.
