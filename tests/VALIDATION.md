# Validation

This document describes **where** the validation data lives and **how**
the tolerances are defined. The actual numbers live in the generated
**[Validation report](docs/validation/report.md)** (`npm run docs`
regenerates it from the suite files) — this document never restates
its counts, so the two cannot disagree.

## Reference data

We validated the library against Halgren's MMFF94 validation suite
(761 molecules, the November 1998 revision). It gives the structure,
total energy, the seven energy components, the partial charges, and
the atom types from the original program.

The molecule-by-molecule comparison is committed as generated
documentation:

- `docs/validation/report.md` — the full census (single source of truth).
- `docs/validation/total-energies.txt` — all 761 totals side by side.
- `docs/validation/per-term-and-charges.txt` — per-term + charge deltas.

Regenerate them with `npm run docs`. CI checks they're current
(`npm run docs:check`).

## Tolerances

A three-layer strategy (full detail in the report):

1. **Hard suite gate** (`tests/compliance-gate.test.ts`, runs in
   `npm run test`): every typing-exact molecule's per-term residual
   ≤ 1e-4 kcal/mol. AN11A/DOZNIP electrostatics excluded (reference
   itself inconsistent). ERULE_03/06 generated-bond rows pin at
   measured tolerances (reference 3-dp print precision).
2. **Pinned-molecule regression rows**: the ERULE fragments' torsion
   totals assert against BatchMin's printed value (catches the
   structurally-blind class of rule change).
3. **Empirical fixtures** (`tests/wittig-ylide.test.ts`,
   `tests/phosphine-imide.test.ts`, three-way fixture comparison):
   chemistry the suite was never stressed on.

Partial charges: 0.001 e per atom. Gradients: finite-difference
checked (δ = 1e-6 Å; relative error < 1e-5; worst 8.5e-8).

## Outliers

The two AN11A/DOZNIP electrostatics exclusions and the JALSOE/SO18A
dative-adjusted charges are the reference's own inconsistencies —
documented in the report. The former FE2PW3/CU1PW1 vdW split is
closed (the +2/+1 cation rows carry their own parameters; the bridge
is by formal charge — `docs/implementer-notes.md` §5.1).

## Method

- Atom typing: `assign_atom_types`, cross-checked against OpenBabel's
  canonical types and the original program's own assignments
  (`docs/implementer-notes.md` §5.4).
- The compliance gate is the load-bearing check for term-level
  regressions; the report is the full evidence. If a number changes,
  regenerate the report — don't edit it by hand.
