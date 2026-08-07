# MMFF94 Compliance Statement

**mmff94-ts** implements the Merck Molecular Force Field (MMFF94) in
pure TypeScript, with zero runtime dependencies. It runs in browsers
and in Node.js, and computes energies, analytical gradients, and
optimized geometries.

The authoritative specification is Halgren's series of papers:
*J. Comput. Chem.* 1996, *17*, 490-641, and the 1999 torsion paper
(*J. Comput. Chem.* 1999, *20*, 720-729). Where the published text and
the original MMFF94 program disagree, this library follows the
program's observed behavior — the validation suite below is the
arbiter, exactly as it is for the rest of the field.

## 1. Which variant

This library implements **MMFF94**, the original force field.

- It is **not** MMFF94s (the "static" variant for conformational
  energies). Most notably, the amide nitrogen keeps MMFF94's
  deliberate pyramidalization (the negative out-of-plane constants) —
  the "s" variant restores planarity, and we do not.
- The in-vacuo dielectric (D = 1.0) is the default and the only
  dielectric option. The D = r solvent model is not implemented.

## 2. Energy terms

All seven MMFF94 terms are implemented, with the published functional
forms and conversion factors:

| Term | Functional form | Reference |
|---|---|---|
| Bond stretch | E = 143.9325·(k_b/2)·Δr²·(1 + cs·Δr + 7/12·cs²·Δr²) | eq. (2) |
| Angle bend | E = 0.043844·(k_a/2)·Δθ²·(1 + cb·Δθ) | eq. (3), eq. (4) for linear centers |
| Stretch-bend | E = 2.51210·[k_sb_IJK·(r₁−r₁₀) + k_sb_KJI·(r₂−r₂₀)]·(θ−θ₀) | eq. (5) |
| Torsion | E = Σ (V_n/2)·[1 + cos(n·τ − γ_n)], n = 1, 2, 3 | eq. (7) |
| Van der Waals | Buffered 14-7: ε_ij·[(1.07·R*/(r+0.07·R*))⁷·(1.12·R*⁷/(r⁷+0.12·R*⁷) − 2)] | eq. (8) |
| Electrostatic | E = 332.0716·q_i·q_j/(r + 0.05) | eq. (6), part III |
| Out-of-plane | E = 0.043844·(k_oop/2)·χ² | eq. (6), part I |

The 1-4 rule matches the spec: electrostatic 1-4 pairs are scaled by
0.75; van der Waals 1-4 pairs are **not** scaled (Halgren 1996,
p. 496). Partial charges come from the bond-charge-increment (BCI)
model (part V eq. 15, including the formal-charge sharing rules).

## 3. Parameter provenance

The parameter tables are extracted from OpenBabel's `.par` text files
(mmffbond, mmffang, mmffstbn, mmfftor, mmffvdw, mmffchg, mmffpbci,
mmffoop), which are mechanical transcriptions of Halgren's published
tables. The extraction script converts format only — it implements no
force-field logic. The values are committed as TypeScript tables, so
the build needs no Python and no external data.

The parameter tables are cross-checked against Tinker's
`mmff94.prm`, an independent transcription of the same originals.

## 4. Validation against Halgren's own suite

The library is validated against the MMFF94 Validation Suite
(753 molecules, the CCL archive version), using the BatchMin 5.5
per-component energies as the reference. Every molecule's atom types
are assigned independently and compared:

- **Atom typing**: 753/753 molecules match OpenBabel's canonical
  MMFF94 types; 749/753 are byte-identical to the original program's
  own types (the four remaining atoms are parameter-inert
  oxidation-state and dative-bonding cases, proven parameter-identical
  by the energy checks).

### 4.1 Per-component energy residuals

For each of the seven terms, the table shows the absolute residual
against BatchMin across the suite (typing-exact molecules; the four
documented anomaly exclusions in section 5 applied):

| Term | Molecules | Worst | Mean | RMS |
|---|---|---|---|---|
| Bond stretch | 753 | 4.97e-5 | 7.0e-6 | 1.0e-5 |
| Angle bend | 753 | 3.45e-5 | 6.0e-6 | 8.4e-6 |
| Stretch-bend | 753 | 4.34e-5 | 4.7e-6 | 7.1e-6 |
| Torsion | 753 | 4.69e-5 | 1.2e-6 | 2.6e-6 |
| Out-of-plane | 753 | 1.64e-5 | 6.0e-7 | 1.6e-6 |
| Van der Waals | 751 | 4.37e-5 | 6.7e-6 | 9.4e-6 |
| Electrostatic | 751 | 6.76e-5 | 7.1e-6 | 1.1e-5 |

At the ±5e-5 gate — the accuracy Wavefunction claims for Spartan —
**745 of 753 molecules match all seven components**. Every component
of every molecule is within 1e-4.

### 4.2 Total energies

The total-energy residuals (the 749 molecules without excluded terms):

| Statistic | Value |
|---|---|
| Worst | 5.10e-4 kcal/mol |
| Mean | 2.48e-4 kcal/mol |
| RMS | 2.86e-4 kcal/mol |

We do **not** claim the Wavefun-level ±5e-5 on totals; the per-term
residuals accumulate to the ~2.5e-4 mean shown.

### 4.3 Partial charges

BCI partial charges match the suite's reference charges to better
than 1e-3 on all 749 charge-comparable molecules.

### 4.4 Gradients

Analytical gradients exist for all seven terms and are
finite-difference checked on every fixture and the pinned suite
molecules (δ = 1e-6 Å; relative error < 1e-5; worst observed 8e-8).

## 5. Known deviations and limitations

1. **The type-76 anionic nitrogen** (AN11A, DOZNIP): the electrostatic
   component diverges across *all* implementations (BatchMin, Tinker,
   OpenBabel, and this library disagree; Tinker drops the term for
   AN11A entirely). Halgren's own papers caution that strongly
   delocalized anions have no uniform charge assignment. These two
   molecules' electrostatics are excluded from the census; their other
   six terms are verified.
2. **The hydrated-metal van der Waals** (FE2PW3, CU1PW1): a parameter
   split that proved to be a typing collapse, closed 2026-08-05 — both
   molecules now match BatchMin to ~1e-6 and rejoin the census.
3. **The empirical bond rule (eq. 18)**: the paper prints a δ = 0.008 Å
   shrinkage plus hybridization corrections. Halgren's own reference
   implementation does not apply them, and the validation suite
   reflects that. This library matches the reference implementation,
   not the literal printed equation. Tinker's transcription agrees.
4. **The empirical torsion rules**: the corroborated universal rule
   (c) reading — every non-aromatic central bond gets eq. (21) with
   π = 1.0 for mltb-2/2 pairs, else 0.4 — matching both OpenBabel and
   Tinker on the measured cases. The paper's rules (d)-(h) are
   unreachable for non-aromatic bonds and are kept as the spec. The
   suite never exercises the empirical torsion rules.
5. **OpenBabel divergences on phosphorus bonds**: on the vinyl
   phosphine C–P bond, OpenBabel's empirical bond length and
   stretch-bend constants differ from the reference behavior; Tinker
   and this library agree to four decimals. The finding is documented
   in `docs/implementer-notes.md` §4.1.
6. **Small-ring and linear-center edge cases**: the near-linear angle
   and out-of-plane guards return the true limit values (the 2026-08
   fixes), and the 3- and 4-ring class parameters follow the published
   class scheme.

## 6. Dependencies

Zero runtime dependencies. The library compiles to ESM in `dist/` and
runs in browsers and Node.js with identical results (verified in
headless Chromium to 1e-7).

## 7. Reproduction

The numbers in section 4 are regenerable:

- `npx tsx tests/scripts/residual-distribution.ts` — the residual
  table and the ±5e-5 gate counts.
- `npx tsx tests/scripts/energy-scoreboard.ts` — the per-molecule
  census.
- The validation suite tests: `npm run test` (the suite comparisons
  are part of the regular test run).
