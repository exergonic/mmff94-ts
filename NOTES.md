# mmff94-ts — Project Notes & Errata

Working journal for context that changes faster than AGENTS.md. Read when
relevant; update freely. Nothing here steers agent behavior — it is context,
not instruction.

## Current status

- Suite status: GREEN (261 passed, 1 intentional skip), typecheck clean —
  2026-08-23 adversarial-review session (see Errata).
- Geometry pipeline details live in the valence-orbital-visualization skill;
  mmff94-ts carries the rigor.

## Errata

_(Corrections and caveats about prior work — dated entries.)_

- **2026-08-23 review session closed the remaining open threads**: B11
  (validation artifacts confirmed 761-row/current), B12 (silent fixture skip
  → loud MISSING-reference failure + `INTENTIONALLY_UNREFERENCED` registry),
  B13 (CUVJOS added to the FD sweep — first charged chemistry, all seven
  terms ≤1e-9), C15 (zero-length-bond guard documented in derivatives.ts),
  C16 (SDF bond indices remapped through declared→compacted slot map +
  discriminating test), D22 (prepare_molecule length-checks stale annotation
  arrays), doc-truth sweep (walkthrough 753→761 counts, implicit-H claim
  replaced with the explicit-H contract, Chromium→Edge, gradient angle-bend
  cb, validation-report generator key bug fixed — real |Δ| now printed).
  New references: dimethyl-ether, phosphine, butane/ethane/water
  non-optimized .mmff94.log (all pass exactly; water at 200.8 kcal/mol
  distortion exercises the cubic terms off-equilibrium).
- **Q3 CLOSED same session via TINKER arbitration**: trp-cage single-point
  on lenovo matches this library ≤0.001 on every term (elec −449.4932 vs
  −449.49394); OB is the outlier at neutral-pH zwitterion termini. The SDF's
  lone digit on the N-term N is a V2000 VALENCE override (not a charge) —
  type 34's q⁰=+1 supplies the +1. Log: tests/references/tinker/trpcage.log.
- **Optimizer perf round (same session)**: cached strong-Wolfe line search +
  topology-keyed nonbonded/class contexts → nicotine optimization 4.28 s →
  1.94 s, bench energy 17→9.1 ms / gradient 24→15.8 ms, bitwise-stable.
  Commits 5edde4d, c18ea93, 61afd81.

## Decisions worth remembering

- WS2+WS3 closed 2026-08-07: empirical rules arbitrated by ERULE rows.
- FAPLUD q⁰(72) split over terminal chalcogens (P(=O)(S⁻) −0.5/−0.5).
- Optimization-test fixtures are user-authored (`*_non-optimized.sdf` from
  Avogadro, source OneDrive/Desktop/structures/).

## Open threads

_(Things to pick up next session.)_

- **Dative-drawing normalization in our typer** (found 2026-08-23 during
  the Tinker sulfone investigation): a SMILES-drawn dative dimethylsulfone
  `C[S+2]([O-])([O-])C` types as S class 15 (THIOL/SULFIDE — wrong) with
  net −1, while the hypervalent drawing `CS(=O)(=O)C` types correctly as
  class 18 + two class-32 with net 0. The suite's own dative cases
  (JALSOE, SO18A) pass because the `.mmd` files use BatchMin's "MMFF
  dative" convention, which our typer reads. Question to settle: should a
  raw formal-charge SMILES normalize identically (S⁺² with 4 single bonds
  → class 18, O⁻ terminal on S(VI) → class 32)? If yes, the fix lives in
  the sulfur/oxygen typing branches of assign-atom-types.ts and needs a
  discriminating test pair like `tests/scripts/dms-two-drawings.ts` (the
  throwaway probe from this session — rewrite it). Related context:
  docs/tinker-sulfone-charge-bug.md documents Tinker's *separate* bug
  (kcharge.f line 230 hard-codes −0.5 on type 107 with no sulfone
  compensation), which is upstream's, not ours.
- **Tinker sulfone bug FILED upstream**: TinkerTools/tinker#185
  (2026-08-24), draft at docs/tinker-sulfone-charge-bug.md. If no response
  in ~2 weeks (mirroring #184's silence), plan was to email Dr. Ponder
  directly with the same text. If confirmed/fixed upstream, revisit the
  ff-bench energy comparison — the spurious-charge molecules should then
  agree with us too.

## Future direction — what industrial-strength would add

_(Mile-high gap list, written 2026-08-24 after the optimizer/benchmark
round. Ordered by tier; each entry notes what implementation likely
entails. The engine itself is complete — everything below is reach.)_

### Tier 1 — claims not yet true (cheap, high-signal)

- **CI (GitHub Actions).** No `.github/workflows/` exists today. A
  workflow running `npm run typecheck && npm test` on push/PR converts
  "tests pass on Billy's machine" into a durable property — during the
  ERULE implementation a couple of broken pushes would have been caught
  by this. Trivial to write (~20 lines); do first.
- **npm publish.** README previously claimed `npm install mmff94-ts`
  while the package was never published. Decision (2026-08-24): publish
  AFTER Tier 2 lands so the first public release carries the real
  feature set. When publishing: set final version, add
  `files`/`exports`/`engines` fields, test the packed tarball.

### Tier 2 — functional gaps real workflows hit

- **Conformer generation.** THE missing capability — MMFF94's main
  industrial use is conformational ensembles. Minimal respectable
  version: rotatable-bond enumeration from the existing connectivity
  graph → torsion grid/sampling → minimize each with the existing fast
  oracle → RMSD-cluster → rank by energy. Everything it needs already
  exists in src/ (fast evaluator, optimizer, dihedral utilities).
  Estimated: one solid session.
- **Constraints & restraints.** Freeze atoms, tether-to-reference,
  distance restraints — required for ligand-in-pocket refinement.
  Implementation: constraint projection after each optimizer step
  (SHAKE-style for frozen atoms is simplest) plus penalty terms in the
  fast evaluator for restraints.
- **Geometry writers** (SDF at minimum). Pipelines can't round-trip
  minimized structures today. parse_sdf exists; the writer is its
  mirror. Small.
- **Multi-record SDF read.** Verify whether parse_sdf reads record 1
  only; libraries process ensembles. If single-record, split-and-loop
  helper or parser extension. Small once confirmed.
- **Cutoffs + neighbor lists (Phase D).** All-pairs O(N²) fine ≤~500
  atoms, dead beyond. The nonbonded-context.ts pair cache from the perf
  round is exactly the substrate: cell-list build per N steps, skin
  margin, per-pair cutoff mask. Opt-in flag; exact mode stays default.

### Tier 3 — input breadth

- **MOL2 reader.** We already understand its dative/hypervalent pitfalls
  (see Open threads above — the two conventions must normalize).
- **PDB reader.** Proteins; pairs naturally with Tier 2 cutoffs if
  protein-scale work is wanted.

### Tier 4 — polish

- Typed error classes (`ParseError`, `TypingError`) instead of generic
  throws; parser fuzzing (malformed-SDF corpus); documented
  worker-thread pattern for browser use.

### Deliberately NOT planned

Periodic boundaries, free-energy methods, polarizability — outside
MMFF94's mission and each is project-sized. The library stays an
excellent molecular-scale MMFF94 engine rather than a mediocre
everything-tool.

---

## REVIEW PREP — adversarial review brief (written 2026-08-13)

**UPDATE 2026-08-13 (same session): the doc/comment mismatches are FIXED and
the torsion regression is CLOSED.** What happened:

1. Regenerating the census exposed a live regression: torsion residuals
   0.29-0.41 kcal/mol on ERULE_01/02/04/08 — introduced 2026-08-11 by commit
   863a70c (the both-pilp rule-(g) case-(1) skip moved outside the mltb gate),
   which contradicted BatchMin itself. Fixed: the both-pilp pair now
   suppresses only rule (g)'s V2 and falls through to rule (h)'s
   V3 = √(V_b·V_c)/N_bc — the reference totals match to 5 decimals
   (probe-erule-skipped.ts; pins rewritten in torsion-empirical.test.ts +
   phosphine-imide.test.ts). Tinker's both-pilp zeroing is a Tinker deviation.
   Census restored: torsion 761/761 ≤1e-4 (worst 4.69e-5 SO18A), vdW 761/761
   (dropped the stale FE2PW3/CU1PW1 exclusions from the scripts), totals
   worst 2.52e-3 (ERULE_03), all-seven ≤5e-5: 753/761, ≤1e-4: 757/761.
2. Doc/comment fixes: van-der-waals.ts ("applied in total.ts" → electrostatic.ts;
   no 0.5 vdW scaling), angle-bend.ts header cb (precise −0.4·π/180), charges.ts
   550→761, mmd-parser.ts (203-recovery codes now documented), VALIDATION.md
   (753→761, 749→757), compliance.md (vdW 761/761, elec-only exclusions,
   gate counts 753/757, totals 759), numerical-precision.md (746/761 reference
   self-agreement), validate-against-suite header, implementer-notes (761),
   AGENTS.md (tree 761, compliance section rewritten), README tree, tool
   comments. Stale-pattern sweep clean (550/753/749/"applied in total.ts" gone).
   Full suite GREEN: 27 files, 240 passed, 6 skipped (the +1 is the new (8,15)
   both-pilp pin).
3. **Still open for the review session (the enforcement gap, now PROVEN)**:
   the per-term ≤1e-4 census is NOT asserted in vitest — the 863a70c regression
   passed `npm test` for two days. Decide whether residual-distribution /
   energy-scoreboard gates belong in the test run (or at least a pinned ERULE
   torsion regression row — the four ERULE fragments are the natural pins).
   Remaining review targets from the original brief below: the implicit-H
   contract, the zero-length-bond gradient guard, the silent it.skip for
   log-less fixtures (dimethyl-ether.sdf is still in that state), FD coverage
   of charged chemistry, and the doc-truth sweep of the remaining files
   (walkthrough.md was not re-audited this round).

### Original brief (goal, standing, architecture, adversarial list)

See the sections below for the full original write-up — the architecture map,
validation numbers, and the remaining adversarial targets (implicit-H contract,
zero-length bond gradient, silent fixture skips, FD charged-chemistry hole,
perf probes) are unchanged and still open.

The user asked the agent to learn the codebase and play the ADVERSARY in a
thorough review. This section is the handoff. Baseline measured 2026-08-13:
`npm run test` = 27 files, 239 passed, 6 skipped, 73 s, vitest 4.1.10.
Last commit af2b844 (parameter_gap_report, 2026-08-12). Working tree: NOTES.md
and tests/fixtures/sdf/dimethyl-ether.sdf are UNTRACKED (dimethyl-ether has no
reference log → silently skipped by reference-comparison.test.ts).

### Project goal and standing (one paragraph)

Pure-TypeScript MMFF94 (zero runtime deps, browser+Node) — energy, analytical
gradients, L-BFGS + steepest-descent optimization. Validated against Halgren's
761-molecule Nov-1998 suite (BatchMin per-term energies, .mmd pchg charges,
OPTIMOL types), OpenBabel, and Tinker. All seven terms ≤1e-4 on 759-761/761
(the ERULE_03/06 generated-bond rows at the reference's 3-decimal print
precision; AN11A/DOZNIP elec excluded — reference-inconsistent, three-way
documented), typing 761/761 vs OB, charges <1e-3 on 757, gradients FD-verified
(worst 8.5e-8), optimizers converge the user-authored *_non-optimized series.
Phase 7 (publish) reached: 0.1.0-alpha.1, `npm pack`-tested; Phase 8 stretch =
per-interaction energy breakdowns.

### Architecture map

- src/types.ts — Atom/Bond/Molecule/TypedMolecule/EnergyComponents/OptimizationResult.
- src/sdf.ts — V2000 parser (explicit-H only; no implicit-H handling).
- src/mmff94/assign-atom-types.ts (1364 ln) — the decision tree; prescans for
  water/carboxylate/amide; separate N pass; aromatic perception
  (find_aromatic_rings: cycle enumeration ≤6-ring, chord check, π=6 eval,
  exo-fixpoint for fused); estimate_ring_size (shortest-cycle BFS heuristic).
- src/mmff94/charges.ts — BCI + part V eq. 15 q⁰/α sharing; environment rules
  for q⁰ of types 32/61/72/76/81; `round: true` default (BatchMin 5-dp).
- src/mmff94/parameter-gaps.ts — NEW (af2b844): flags atoms whose coordination
  EXCEEDS the type's crd + non-C atoms on the type-1 fallback. Diagnostic only.
- src/mmff94/parameters/ — data tables (bond/angle/torsion/bci/vdw/oop/...
  auto-generated, committed) + lookup.ts (priority-wildcard helper) +
  parameter-classes.ts (BTij/ATijk/TTijkl/STijk, class-scoped resolution,
  torsion canonical-direction base-136 trick) + empirical.ts (part V: eqs.
  18-19 bonds, θ₀ protocol + eq. 20, torsion rules; tables radii/χ/Z/C/U/V/W;
  deliberate deviations: χ(P)=2.04/N=3.05, V(S)=0.48).
- src/mmff94/energy/ — 7 term files + total.ts (all pure, each rebuilds its own
  ClassContext via make_class_context).
- src/mmff94/gradient/ — mirrors energy/; shared resolution helpers
  (bond_length_derivatives, stretch_bend_angle_terms, torsion_terms,
  vdw_pair_parameters, oop_force_constant, is_1_4_pair); degenerate-geometry
  guards fall back to forward differences at cusps.
- src/optimize/l-bfgs.ts — N&W 7.5 + strong-Wolfe/cubic zoom; γ-compensated
  first trial, noise-pair discard, −g fallback + history reset. steepest-descent.ts.
- src/utils/ — vector.ts (incl. dihedral_angle with the fixed v3×v2
  handedness; wilson_oop_angle), mmd-parser.ts (suite structures, MMD_ELEMENT
  incl. the 203-recovery codes 201/204/206/65-67/70/207-212).
- tests/ — 27 files; key: validate-against-suite (suite per-term vs BatchMin),
  reference-comparison (fixtures vs obenergy logs, 0.02 gate, vinylphosphine
  exceptions pinned), gradient (FD δ=1e-6), optimization (non-optimized
  fixtures), atom-types-suite (761/761 + KNOWN_GOOD), charges-suite,
  openchemlib-comparison, browser-smoke (dist in headless Edge), wittig-ylide,
  parameter-gaps. Docs: AGENTS.md (authority, gitignored), README (STE),
  tests/VALIDATION.md (STE, facts-only), docs/implementer-notes.md (gore),
  docs/mmff94-compliance.md (Wavefun-model statement), docs/walkthrough.md,
  docs/numerical-precision.md.

### Adversarial target list (found during familiarization — VERIFY each)

**A. Docs vs code (the repo's own standard: docs must match code exactly)**

1. van-der-waals.ts L39 + L152-153: comments claim the 1-4 elec ×0.75 scaling
   is "applied in total.ts" and mention a "0.5 scaling applied externally" —
   both false: no vdW 1-4 scaling exists (MMFF94 spec), and the 0.75 lives in
   electrostatic.ts. Same file's own header (L36-39) is correct. Stale from an
   earlier design.
2. angle-bend.ts header eq. block L12: "cb = −0.007 deg⁻¹" — the code
   deliberately uses the precise −0.4·π/180 (L41, documented below). The header
   equation documents a constant the implementation does not use.
3. tests/VALIDATION.md L12/L21-22: "The suite has 753 molecules", "all 753
   totals" — stale since the 761 upgrade (2026-08-07). L100: charges "749
   molecules" vs compliance's 757.
4. docs/mmff94-compliance.md §4.1 vdW row: "759 of 761" with worst 4.35e-5 —
   internally contradictory (worst ≤1e-4 ⇒ all 761 green) and stale (the
   FE2PW3/CU1PW1 vdW exclusions closed 2026-08-05; README says 761/761). The
   §4.1 "751 of 761 ≤5e-5 all seven" and "755 of 761 ≤1e-4 every comparable
   term" counts need regeneration (residual-distribution.ts) to verify.
5. tests/validate-against-suite.test.ts header: "241/550 in the suite
   scoreboard" — stale (761/761 now).
6. README structure tree: "mmff94-atom-types.json … (550 molecules …)" —
   stale; README's own validation section says 761.
7. charges.ts header: "all 550 reference molecules" — stale.
8. mmd-parser.ts L24-26 comment: "Types not listed here (31, 51, 201+) are
   ones OpenBabel cannot translate either" — contradicted by the table below
   (L40-42), which lists exactly those codes (the 203-recovery amendment).
9. AGENTS.md §MMFF94 Compliance Statement still cites the ver.98.05.22 753
   URL and a 0.01 kcal/mol target; its structure tree nests
   parameters/parameters/ (cosmetic). Grep the whole doc set for stale
   patterns: 550, 753, 749, "applied in total.ts", "0.5 scaling", "(planned)".

**B. Claims vs enforcement (the big one)**

10. The headline per-term ≤1e-4 census is NOT asserted by any vitest test.
    validate-against-suite.test.ts only PRINTS the per-term stats and asserts
    nChecked ≥ 50; the real gates are manual scripts (energy-scoreboard.ts,
    residual-distribution.ts, validation:doc). Hard gates that DO exist:
    761/761 typing, 8 oop pins, FUVDOP/FILNOD/JIYJAC, fixture comparisons
    (0.02), charges pins, FD gradients. Verdict to reach: a term regression on
    an unpinned molecule passes `npm test` while the compliance docs stay
    stale — should the census be a test?
11. docs/validation/*.txt: committed generated artifacts — check whether they
    currently hold 753 or 761 rows (the generator is the truth per VALIDATION.md).
12. reference-comparison.test.ts L62-65: a fixture with no reference log is
    skipped SILENTLY (it.skip at collection). dimethyl-ether.sdf is in exactly
    that state right now — untracked, never actually compared. New fixtures
    can quietly never be tested.
13. Gradient FD coverage: fixtures + exactly 3 suite molecules
    (FUVDOP/FILNOD/JIYJAC) — the documented residual hole: charged/ionic
    chemistry is never FD-tested (skill's own note; the natural patch is a
    charged suite molecule in SUITE_MOLECULES). Verify no charged fixture is
    in the FD sweep.

**C. Consumer-input robustness (what a real caller can feed it)**

14. Explicit-H contract: typing counts EXPLICIT neighbors only. An implicit-H
    SDF (common from sketchers/converters) types silently wrong (bare C →
    fallback type 1, no H types) and yields nonsense energies. parameter_gap_report
    does NOT catch it (fewer-neighbors-than-crd is deliberately unflagged).
    README never states the explicit-H requirement. Recommended review
    verdict: document the contract + consider a validation flag.
15. Duplicate/overlapping atoms (r ≈ 0): bond gradient returns zero silently
    (derivatives.ts L47 guard) — energy is huge-but-finite, gradient vanishes,
    optimizer stalls at high E without NaN. The angle/oop cusps got the
    forward-difference treatment; the zero-length bond is the remaining silent
    hole. (r=0 is a genuine cusp — the one-sided limit is direction-dependent;
    the guard value is defensible but should be documented.)
16. Malformed SDF: parse_sdf skips bad atom lines and compacts indices while
    bond block offsets use the declared count — atom/bond index misalignment
    on malformed input; no fuzz tests. Minor, but the parser is the front door.

**D. Chemistry/algorithmic spots to re-derive**

17. The empirical-rules arbitrations rest on a handful of 3-decimal reference
    rows: χ(P)=2.04/N=3.05 (P–Si, F–N: Δχ 0.30/1.05), V(S)=0.48 ((S,C) and
    (N,S) rows), rule-(c) order-2 gate (P–Si V3=0.285), the both-pilp
    case-(1) skip (Tinker-arbitrated, suite-blind "green both ways" — the
    skill itself flags this class). For each: how many rows pin it, and would
    the alternative value pass the suite too?
18. q⁰ environment rules (types 32/61/72/76/81) are empirical fits to the
    suite's cases (e.g. the ester-enolate C(=O)(O⁻)(OR) q⁰(32) path is
    unvalidated). The landscape probes exist (probe-s72-landscape.ts); the
    review should check coverage claims, not just values.
19. Perf: every energy term (and every gradient term) rebuilds ClassContext
    via make_class_context — 7+7 perceptions per energy+gradient call, each
    running the aromatic cycle enumeration (DFS ≤6-ring from every ring atom,
    dedup by sorted key). Bench covers only trp-cage (304 atoms, ~17 ms
    energy). Probe a cage-rich molecule (C60, porphyrin) for blow-up.
20. electrostatic/vdw pair loops are O(n²) with a per-pair BFS for 1-4
    classification (is_1_4_pair) — fine at 304 atoms; no neighbor lists.
21. lookup_torsion's base-136 mixed-radix direction trick: verify the
    canonical-direction choice matches the par's storage direction for
    asymmetric entries beyond the pinned tests (the chain is
    exact → lvl3/lvl5 asymmetric — a wrong canonical direction hits wildcard
    defaults silently).
22. prepare_molecule: `if (!prepared.atom_types)` — a Molecule carrying an
    empty/stale atom_types array skips typing silently. Weak contract.
23. assign_bci_charges round:true default — energies depend on 5-dp charge
    rounding (BatchMin artifact); documented + opt-out; the review should
    confirm the gradient uses the same rounded values (it does — same field).

**E. Review protocol for the next session**

1. Baseline: `npm run test` (73 s) + `npm run typecheck` (vitest doesn't
   typecheck). Enumerate the 6 skips and justify each.
2. Regenerate the census: `npx tsx tests/scripts/residual-distribution.ts` and
   `tests/scripts/energy-scoreboard.ts`, compare against README /
   VALIDATION.md / compliance tables — this adjudicates findings A3/A4/B10-11
   empirically.
3. Grep sweep for stale numbers/claims: 550, 753, 749, "applied in total.ts",
   "0.5 scaling", "(planned)", "NOT YET IMPLEMENTED", "STUB".
4. Read the review-rounds reference (references/review-rounds.md in the skill)
   for round-1 precedent: verify every claim against code, write the
   discriminating test, one commit per concern.
5. Deliverable: a triage table (finding → verdict → action) — the user decides
   what gets fixed. User report format: one actionable with its definite
   payoff at the end; hedge in docs, not in the reply.

## Optimizer Phase A (fast path) — 2026-08-23 (commit 3696b71)

- New src/optimize/fast-system.ts: compiled kernel (typed arrays, zero
  per-call allocation) used by both optimizers on their default path.
- Differential guard tests/fast-system.test.ts: energies bitwise-equal
  to the readable terms; per-term gradients <=1e-8 absolute.
- Measured speedups: nicotine 1925 -> 124 ms, butane 118 -> 9 ms,
  ethane 41 -> 7 ms, water 13 -> 4 ms; trpcage oracle ~21.4 -> ~1.2 ms.
- Iteration counts can shift on flat surfaces (ULP chaos) — same basin,
  energies agree; the optimization series still converges every fixture.
- Phase B (arithmetic micro-opts: pow->mults, hypot->sqrt) and Phase C
  (RMS/max-step convergence criteria, warm-started α) remain.
