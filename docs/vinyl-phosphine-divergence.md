# Vinyl phosphine: one torsional resolution accounts for the whole energy gap

Status: characterized 2026-09-18. No action — the resolution is deliberate,
and this note shows it is the entire difference.

## Summary

Vinyl phosphine minimizes to **6.5716 kcal/mol** in this library and to
**10.4515** in Tinker, both from the same starting geometry. The two C–P
dihedrals resolve through rule (g) of the empirical torsion rules here —
π = 0.15, V2 = 1.423 kcal/mol. The references resolve that same bond to
V2 = 3.795 in their per-interaction output. The paper's rule (g) is what we
implement; 3.795 is their deviation (compliance §5). The 3.9 kcal/mol gap
is that one resolution, and nothing else.

## The single-geometry evidence

At the fixture geometry (`tests/fixtures/sdf/vinylphosphine.sdf`) the two
sides disagree in the torsion term and essentially nowhere else:

| term | this library | references |
|---|---|---|
| torsion | 2.29399 | 6.1102 |

The difference, 3.82 kcal/mol, is the size of the minimized gap. The
remaining terms sit within hundredths of the references at this geometry
(bond 0.019, stretch-bend 0.003 against Tinker), and the companion angle
item from the same investigation has since closed in the build (our angle
term 0.64645 against their 0.64625/0.6464).

## The experiment

Overriding only the two rule-(g) sites in a copy of the build — the
second-row π = 0.15 → 0.4, i.e. V2 = 1.423 → 3.795 — and re-minimizing
from the identical start:

| build | torsion at the fixture | minimized total | H–P–H |
|---|---|---|---|
| this library | 2.29399 | 6.5716 | 96.7° |
| this library, V2 = 3.795 override | 6.11020 | 10.4571 | 101.1° |
| Tinker | 6.1102 | 10.4515 | ≈101° |

The override lands on the reference torsion at the fixture geometry to
five decimals, and on the reference minimum to 0.0056 kcal/mol. It also
lands on the reference geometry: the P center tightens from 96.7° to
101.1° of H–P–H, which is where Tinker's minimum sits.

## What the residual is

0.0056 kcal/mol is not a second disagreement. It is the transposition
drift in the remaining terms (each ≤0.02 against Tinker at the fixture
geometry) plus the convergence tolerance — the two optimizers stop at
different RMS gradients. The torsion was the only term large enough to
matter, and the override consumed all of it.

## Reproduction

Minimize vinyl phosphine from a fetched geometry with this library, then
repeat with a build whose rule-(g) case (2)/(3) resolution uses π = 0.4.
The two numbers should be 6.5716 and 10.4571. The measurements in this note
were taken through the WebMO engine integration (`integrations/webmo/`).
