# Validation Ledger

Living record of what `mmff94-ts` has been validated against —
molecule by molecule, term by term. The numbers here are produced by
`npm run test`; update this file whenever the reference tests change.
The README's [Validation](../README.md#validation) section is the
condensed public version of this ledger.

## Current status (2026-08-03, night)

- **Atom typing**: **753/753 suite molecules type-exact vs OpenBabel
  (100%)** — the full Halgren suite. The 550 exact codes from the
  pre-recovery plateau are pinned as `KNOWN_GOOD` in
  `atom-types-suite.test.ts`. Three jumps today: the phosphate/sulfonate
  oxygen cluster (140 →
  241, morning), the amidine/sp-N cluster (241 → 512, evening), and
  the aromatic-perception cluster (512 → 550, night): the π-count
  Kekulé rule with the fused-aromatic exocyclic-double allowance
  (fixpoint over the candidate rings — benzimidazole/furan/isoxazole/
  thiophene fusions found, sulfolenes and COGDEH's triazine-fused
  ring excluded), sulfone/sulfoxide S's dropped from the lone-pair
  donor role, nitroso N 46, sulfonyl-imine N 48, triazene N 10,
  P=C carbon 3, allene C 4, the amidinium 57/55 refined by the =N
  coordination with the pyridinium exclusion (COJFIQ's adenine), the
  imidazolium C 80 via the N-type gate (CUDREY's 2-aminoimidazolium),
  and the hydroxide O 35 (OHMW1).
- **Energy terms** vs BatchMin on all 753 typing-exact molecules:
  **all 7 terms 745/745 exact** (within 0.05 kcal/mol; worst |Δ| 0.016,
  KESNEB oop) on the reproducible set — 8 reference anomalies excluded
  (the six part-V delocalized anions AN11A/DAKBAS/AN06A/AN08A/TAJVUV/
  DOZNIP, whose reference charges eq. (15) cannot reproduce — Halgren's
  own vinyl-oxide/vinyl-sulfide caveat — and FE2PW3/CU1PW1, whose
  BatchMin vdW predates the X94 metal parameters; OB and Tinker agree
  with our transcription there, 55.84481 vs 55.8448). The 18 "parameter-gap" molecules from the evening
  session were diagnosed one by one against the BatchMin per-component
  log, OpenBabel's per-interaction HIGH-verbosity log
  (`tests/scripts/probe-ob-log.py`) and Tinker's independent mmff94.prm
  transcription (`temp_tinker/`) — **every one was a lookup or
  constant bug in our code, not a missing table entry**:
  - oop (8 molecules): the term lacked the part-I p. 513 step-down
    chain — the EqLvl3 equivalence of the substituent types, one at a
    time, re-sorted (COYVIV's delocalized N(40) with [28,28,63]
    resolves to the 2-40-28-28 entry via EqLvl3(63) = 2, k = −0.007,
    not the −0.005 wildcard). Level 4/5 must NOT be tried (FUDPOJ's
    cyclopropenone centers would wrongly resolve 1-3-2-7 / 2-2-3-5).
  - strbnd (9 molecules): the class-0 lookup scanned priorities 0–5
    and grabbed entries of other classes — VIYPAU's O–C(20)–C(3)
    angle resolved the class-4 ring entry 4-6-3-20 (k = 1.179)
    instead of the element-row default. The lookup is now locked to
    the computed class (`max_priority = 0`).
  - bend (GESNIB): the cubic-bend constant was the rounded −0.007
    deg⁻¹; BatchMin uses the paper's "more precise" −0.4 rad⁻¹
    (i.e. −0.4·π/180). GESNIB's near-linear C(37)–C(37)–C(22) angles
    (Δ ≈ 47°) differ by 0.035 kcal/mol each between the forms.
- **Partial charges** vs the suite's reference values: **745/745** to
  < 10⁻³ e⁻ per atom on the reproducible set (the 8 anomalies above). The default-BCI sign convention is part V
  eq. 14 (the unparametrized pair's increment flows TO the smaller
  type — the hydroxide's 21–35 pair, OHMW1); the parametrized pairs
  keep the par's own direction. The imidazolium N q⁰ is
  environment-dependent (1/(#N3 on the central C80): 1/2 plain, 1/3
  2-aminoimidazolium), and CIM+ (80) is treated as aromatic for the
  BTij/Bci class (its par entry lacks the arom flag).
- **The 16 obenergy reference molecules**: all seven terms and the
  total exact to 5 decimals.
- **Gradients**: every analytical gradient FD-verified on every atom
  of every reference molecule (worst relative error 8.5×10⁻⁸).
- **Tests**: 187 passing, 0 skipped.

## Reference molecules vs obenergy (16 molecules, kcal/mol)

Asserted in `reference-comparison.test.ts` (exact to 5 decimals).
Deltas are `ours − obenergy`; every row is all zeros.

| molecule | bond | angle | strbnd | torsion | vdw | elec | oop | total |
|---|---|---|---|---|---|---|---|---|
| ammonia | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| benzene | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 16.22697 |
| butane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | −5.07596 |
| cyclohexane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | −3.56091 |
| ethane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | −4.73436 |
| ethene | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 8.20010 |
| formaldehyde | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.05416 |
| formamide | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | −30.90735 |
| methane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.02638 |
| nicotine | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 30.25415 |
| piperidine | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | −5.99401 |
| propane | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | −4.89729 |
| pyridine | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 15.52345 |
| pyrrole | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 3.28684 |
| trimethylamine | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 7.36863 |
| water | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |

The `total` column is our value (== reference). Per-atom partial
charges also match the obenergy logs (`charges.test.ts`, |Δ| < 1e-4).

## Suite — out-of-plane pins vs BatchMin (8 molecules, kcal/mol)

Asserted in `validate-against-suite.test.ts` (|Δ| < 0.05). These
molecules are pinned beyond the typing-exact set because they exercise
specific oop chemistry. AMHTAR01's −0.021 gap closed when the
formally-charged typing landed (2026-08-02, commit `4b13f34`): its
CO₂M carbon/oxygen pair now types 41/32 like the reference.

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

## Suite — per-component report (typing-exact molecules)

Regenerated by `npm run test` (`validate-against-suite.test.ts`
"reports per-component energies for typing-exact molecules"): per-term
deltas vs BatchMin on the currently typing-exact molecules (140 of
550, computed fresh from `mmff94-atom-types.json`). Snapshot
2026-08-02, after the formally-charged typing (`4b13f34`) and the
sulfinate S=O bridges (`86cc3c8`, `5834bce`):

| term | \|Δ\|≤0.05 | mean\|Δ\| | max\|Δ\| | worst |
|---|---|---|---|---|
| bond | 140/140 | 0.000 | 0.00 | — |
| angle | 140/140 | 0.001 | 0.02 | FAGVEO −0.02, DUBNET +0.01, DAWYUV −0.01 |
| strbnd | 140/140 | 0.000 | 0.00 | — |
| torsion | 140/140 | 0.000 | 0.00 | — |
| vdw | 140/140 | 0.000 | 0.00 | — |
| oop | 140/140 | 0.000 | 0.01 | COYNAF +0.01, CITNOI10 +0.01 |
| elec | 140/140 | 0.000 | 0.00 | — |

The comparison is valid: BatchMin's log is a single-point calculation
at the .mmd geometry (suite README), so deltas on typing-exact
molecules are term bugs by construction. The sulfinate family
(JALSOE/SO18A) is the instructive recent case: its S=O oxygen is typed
7 by the reference typing rules but keyed 32 in every parameter table
(OB, TINKER, OpenChemLib agree), so the angle, vdW, and bond lookups
carry documented 7→32 bridges; with them, all three terms match
BatchMin exactly on both molecules.

### How each thread closed (history)

The per-term threads below were closed across the parameter-class work
(commits `9997d0d`–`e75805c`); kept here as the design record.

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
  Both closed later by the 3-ring i = l "torsion" skip (FUVDOP) and
  ring-aromaticity classing (FILNOD).

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
   All fixture torsions now match obenergy exactly; the torsion
   assertion is extended to every fixture. The RING = AL log column is
   informational — the torsion class (TTijkl, part IV p. 609) is 0 for
   six-membered rings; class 5 applies only to five-membered rings.

3. **AMHTAR01 oop (Δ 0.021).** — **RESOLVED 2026-08-02** (commit
   `4b13f34`). Ester/carboxyl pyramidalization; the reference types
   the CO₂M carbon 41 and oxygens 32/32. The formally-charged typing
   now reproduces those types, and AMHTAR01's oop matches BatchMin to
   −0.000003 (see the pin table above).

4. **AGLYSL01 components.** Stretch-bend 0.00000 vs 0.24423, torsion
   −3.00231 vs −4.71331, vdW 3.51572 vs 2.78652 — all traceable to
   typing (sulfonate S, ammonium N); not typing-exact, so not in the
   per-component report. The OOP row is exact (0/0).

5. **Electrostatic — RESOLVED (2026-08-01, term live; fully closed
   2026-08-02, commit `4b13f34`).** The BCI charges (charges.ts) and
   the buffered Coulomb term (eq. 6, part III: E = 332.0716·q_i·q_j/
   (r + 0.05), the electrostatic buffering constant S = 0.05 Å) match
   the reference logs per-atom AND per-energy for every typed fixture
   — benzene +3.07810, ethene +8.0530, pyridine +2.0939, pyrrole
   +3.0720, nicotine −2.2135 — and the fixture totals match exactly
   (they were the last gap). The suite's per-component report covers
   electrostatics: **512/512** typing-exact molecules match BatchMin's
   electrostatic component (the earlier 89/91 with the two
   metal-carboxylate misses closed when the formal-charge model of
   part V eq. 15 — primary charges + negative-charge sharing, with the
   sharing flowing from each *neighbor's* α — landed; the imidazolium
   q⁰ rule and the CIM+ arom-flag fix closed the rest). Pairs 1-2/1-3
   are excluded (ammonia's electrostatic is zero — the BatchMin log
   confirms), 1-4 pairs ×0.75 inside the term.

6. **Typing scoreboard.** 512/550 suite molecules type-exact vs
   OpenBabel (93.1%; see `atom-types-suite.test.ts`). The climb:
   aromatic-ring perception (Kekulé pattern + lone-pair heteroatoms in
   5-rings) took 65 → 91 (2026-08-01); amide-N typing (types 10/28,
   with the sulfonyl exclusion) and the formally-charged types
   (34/35/51/54/55/56/57/58/68/72/73/77/80/81/89/90/91, structural
   valence-based branches) took it to 140 (2026-08-02); the
   phosphate/sulfonate oxygen cluster (2026-08-03 — type-32 rules for
   P=O/S=O/S–O⁻ keyed on the terminal-oxygen count, H-on-P 71,
   acid/enol/phenol H 24/33/29, H-on-cationic-N 36, C=S → 3,
   sulfine 74, P 25/75/26) took it to 241; the amidine/sp-N cluster
   (same day — C=N → 3, amidinium 57 by =N coordination, N3-on-sp2-C
   40, nitrile 42, sulfonamide/cyanamide 43, H rules 27/28/36 via an
   N-type pass, 4-ring alkene 30, cyclopropane BFS fix) took it to
   512; the aromatic-perception cluster (same night — the π-count
   Kekulé rule with the fused-aromatic exo-double fixpoint, sulfolene
   exclusion, nitroso 46, sulfonyl-imine 48, triazene 10, P=C 3,
   allene 4, amidinium refinements, imidazolium-C 80, hydroxide 35)
   took it to **550 — 100%**. Energy terms match BatchMin on all but
   the 18 parameter-gap molecules listed in Current status.

## Updating this ledger

1. Run `npm run test`; the reference-comparison tests print every
   molecule's per-term comparison, and the suite test prints the
   per-component sweep.
2. Update the fixture delta table and the oop pin table with the new
   numbers. The oop values can be regenerated with a throwaway test
   (parse `MMFF94.mmd`, `calc_oop_energy`, read `Out-of-Plane =` from
   `MMFF94_bmin.log`, remembering the Fortran D exponents).
3. When an open question is resolved, strike it with the commit that
   fixed it, and tighten the corresponding assertion in the test
   files.
