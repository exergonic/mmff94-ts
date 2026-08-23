# Implementer's Notes: Validation Reference

A concise map of where validation knowledge lives. The narratives that
once filled this file now live in code commits (the "why") and the
generated report (the numbers). This file points at them.

## 1. The three type numberings

Conflating these has cost real debugging hours:

| numbering | range | where |
|---|---|---|
| **original Halgren** | 1–178 | Tinker `mmff94.prm` first column; Merck native. Water O = 177. |
| **class / OpenBabel** | 1–110 | What this library assigns; OB canonical types; OPTIMOL `# ty`. |
| **bmin / OPTIMOL-internal** | — | The `.mmd` first column. **Not public** — do not feed to anything. |

## 2. The three representations

| file | representation | used by |
|---|---|---|
| `MMFF94.mmd` | hypervalent (S=O doubles) | BatchMin; our energy census and typing |
| `MMFF94_dative.mol2` | dative (four S singles) | OPTIMOL; Tinker runs |
| `MMFF94_opti.log` | dative | original program's per-atom types and symbols |

Cross-representation comparisons must align by **coordinates**, not index.

## 3. Reference resources

| resource | location |
|---|---|
| Validation suite (761 molecules) | `tests/fixtures/validation-suite/` (mmd, mol2s, bmin log, energies, titles, atom-types.json) |
| `MMFF94_opti.log` (24.7 MB) | CCL archive only; local copy at `<TINKER_ROOT>/MMFF94_opti.log` |
| Tinker build + `mmff94.prm` | `<TINKER_ROOT>/` (local) |
| OpenBabel `.par` files | `temp_ob/data/` (OB install) |

Scripts:

| script | purpose |
|---|---|
| `tests/scripts/energy-scoreboard.ts` | per-term coverage gate |
| `tests/scripts/residual-distribution.ts` | per-term ≤10⁻⁴ census + worst residuals |
| `tests/scripts/generate-validation-doc.ts` | writes `docs/validation/` (`npm run docs`) |
| `tests/scripts/bmin-log.ts` | shared BatchMin log parser |
| `tests/scripts/gen-tinker-fixtures.ts` | Tinker xyz/prm/key for 16 fixtures |
| `tests/scripts/tinker-fixture-comparison.ts` | three-way fixture table |
| `tests/scripts/ob_energy_breakdown.py` | OB per-term + per-interaction diagnostic |

---

## Where the narratives live

The empirical-rule closures, reference anomalies, and transcription
divergences that used to fill this file are now documented where they
cannot drift:

| topic | lives in |
|---|---|
| Empirical-rule closures (OHMW1 bond, θ₀ protocol, torsion rules, BCI fallback) | `empirical.ts` and `charges.ts` code comments + their commit messages |
| Reference anomalies (FE2PW3 vdW, AN11A/DOZNIP elec, JALSOE/SO18A charges, FAPLUD q⁰) | `docs/validation/report.md` (generated) + code comments that closed them |
| P=N phosphine imide typing + both-pilp torsion fix | `empirical.ts` comments + `tests/phosphine-imide.test.ts` + `tests/torsion-empirical.test.ts` |
| Parameter-resolution subtleties (TAJSUS, step-down chains, OOP pins) | `parameter-classes.ts` and `out-of-plane.ts` code comments |

When investigating a discrepancy, read the code comment in the
implementing file first — the narrative is there, with the commit hash
that introduced it.

## Known transcription divergences

| topic | status |
|---|---|
| eq. (18) δ-transcription | closed — plain form is right (empirical.ts comment) |
| eq. (20) degrees² | closed upstream (OB PR#2741669); we and Tinker square correctly |
| torsion rules (c)/(g)/(h) | paper-arbitrated (empirical.ts comment); ERULE rows confirm |
| Table X V(S) = 0.48 | suite-arbitrated (empirical.ts comment) |
| χ(P) = 2.04, χ(N) = 3.05 | suite-arbitrated (empirical.ts comment) |
| metal vdW parameters | closed — typing collapse, formal-charge bridge (van-der-waals.ts comment) |
| type-76 charges | open — all three implementations differ (report.md) |
| OB 3.2.1 API (per-term methods removed) | worked around in `ob_energy_breakdown.py` |
| angle cubic constant | closed — OB uses rounded −0.007; we use precise −0.4·π/180 (angle-bend.ts comment) |

## Regeneration recipes

Every number in the validation docs, regenerated:

```bash
npm run test                      # the suite + compliance gate
npm run docs                      # docs/validation/report.md, total-energies.txt, per-term-and-charges.txt
npm run docs:check                # CI guard: fails if committed report differs from fresh run
npx tsx tests/scripts/energy-scoreboard.ts   # per-term coverage
npx tsx tests/scripts/residual-distribution.ts  # ≤10⁻⁴ bins + worst residuals
```

The Tinker runs need `<TINKER_ROOT>` with the built `analyze` binary
and the opti log present.

## Open questions

1. **The metal vdW**: which parameter set does Halgren's part III table specify — the X94 revision or the Merck original?
   (Resolved in practice 2026-08-23: the formal-charge bridge — type 88/+2 →
   row 87, type 98/+1 → row 97, `van-der-waals.ts` — reproduces BatchMin's
   hydrate energies to ~1e-6; the historical question is which table Halgren
   *printed*, not which one the reference program uses.)
2. **The type-76 q⁰**: Halgren's own caveat — no uniform primary charge reproduces the reference charges. Three implementations, three answers.
3. **The zwitterion electrostatics gap vs OpenBabel** (found 2026-08-23,
   trp-cage): OB evaluates the N-terminal NH₃⁺ (type 34) with q = +0.333
   per H and the COOH carbon at +0.6; this library computes −0.844 on the
   NH₃⁺ N and +1.2 on a carboxylic-acid C(57) — a ~147 kcal/mol elec gap on
   trp-cage that the BatchMin suite never exposes (all 761 suite molecules
   close to ≤1e-4). The BCI q⁰/α reading of types 34/57 for neutral-pH
   terminal groups is the open question; see `tests/fixtures/sdf/trpcage.sdf`
   (kept unreferenced in reference-comparison until resolved).
