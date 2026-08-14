# MMFF94 Validation Report

> **Generated:** 2026-08-14 from the 761-molecule MMFF94 validation suite (November 1998 revision).
> **Regenerate:** `npm run docs` — this file is the single source of truth; all prose docs point at it.

---

## At a glance

| Claim | Result |
|---|---|
| Typing-exact molecules | 761/761 (100.0%) vs OpenBabel |
| Molecules in suite | 761 |

---

## Per-term energy residuals

Every typing-exact molecule's per-term residual vs BatchMin 5.5 (kcal/mol),
computed at the .mmd geometries. The 1e-4 gate is the hard regression
threshold (`tests/compliance-gate.test.ts`, run in `npm run test`).

### Gate summary

| Term | ≤1e-5 | ≤5e-5 | ≤1e-4 | Worst | Worst molecule |
|---|---|---|---|---|---|
| bond | 577/761 (75.8%) | 759/761 (99.7%) | 759/761 (99.7%) | 1.97e-3 | ERULE_03 |
| angle | 614/761 (80.7%) | 761/761 (100.0%) | 761/761 (100.0%) | 4.30e-5 | BEVJER10 |
| strbnd | 677/761 (89.0%) | 760/761 (99.9%) | 760/761 (99.9%) | 3.36e-4 | ERULE_03 |
| torsion | 755/761 (99.2%) | 761/761 (100.0%) | 761/761 (100.0%) | 4.69e-5 | SO18A |
| oop | 755/761 (99.2%) | 761/761 (100.0%) | 761/761 (100.0%) | 1.61e-5 | ARGIND11 |
| vdw | 575/761 (75.6%) | 761/761 (100.0%) | 761/761 (100.0%) | 4.35e-5 | MG2PW3 |
| elec | 595/759 (78.4%) | 755/759 (99.5%) | 759/759 (100.0%) | 7.08e-5 | DONFOB |

### Coarse-precision exceptions

These rows pin at measured tolerances — the reference prints the
generated parameter to 3 decimals, so the residual is bounded by
the reference's own print precision:

- **ERULE_03** bond_stretch: |Δ| = 0.00e+0 (≤ 2.00e-3 — generated P–Si bond at reference 3-dp print precision)
- **ERULE_03** stretch_bend: |Δ| = 0.00e+0 (≤ 4.00e-4 — inherited from the P–Si generated bond)
- **ERULE_06** bond_stretch: |Δ| = 0.00e+0 (≤ 2.00e-3 — generated F–N bond at reference 3-dp print precision)

### Documented exclusions

- **AN11A / DOZNIP** electrostatics: reference itself is inconsistent
  (the type-76 anionic nitrogen — Halgren's own caveat). All other
  terms verified.

---

## Total energies

| Gate | Count |
|---|---|
| ≤1e-4 | 757/761 (99.5%) |
| ≤1e-3 | 758/761 (99.6%) |
| Worst | 5.81e+0 (AN11A) |

\*-marked rows: BatchMin diverges from OPTIMOL (single-precision
charge sharing — up to 0.0035 kcal/mol).

---

## Partial charges

| Gate | Count |
|---|---|
| max|Δq| ≤ 1e-3 e⁻ | 757/757 (100.0%) |
| Worst | 0.00e+0 () |

Gated on typing-exactness; JALSOE/SO18A/AN11A/DOZNIP excluded
(dative-adjusted or delocalized-anion references).

---

## Gradients

Analytical gradients for all seven terms are finite-difference checked
on every fixture and the pinned suite molecules (δ = 1e-6 Å; relative
error < 1e-5; worst observed 8e-8).

---

## Outliers

Two molecules have one energy term that cannot be reproduced: for each,
the reference itself is inconsistent for that term. All other terms of
these molecules are verified against two independent implementations
(Tinker and OpenBabel — both agree).

| Molecule | Term | Reason |
|---|---|---|
| AN11A | Electrostatic | The anionic five-ring nitrogen has no uniform primary charge (Halgren states this). Each implementation gives a different value. |
| DOZNIP | Electrostatic | Same as AN11A. Tinker drops the term entirely. |

---

## Method

- Structures: `tests/fixtures/validation-suite/MMFF94.mmd`
- Reference per-term energies: `tests/fixtures/validation-suite/MMFF94_bmin.log`
- Reference total energies: `tests/fixtures/validation-suite/MMFF94.energies`
- Reference atom types: `tests/fixtures/validation-suite/mmff94-atom-types.json`
- Raw data: `total-energies.txt`, `per-term-and-charges.txt`
- Hard gate: `tests/compliance-gate.test.ts` (runs in `npm run test`)

