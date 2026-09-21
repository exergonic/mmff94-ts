# Conformer generation for mmff94-ts — design

**Status:** proposal (rev 2), awaiting approval · **Date:** 2026-09-20
**Goal:** add `generate_conformers()` — the last feature before 1.0 and npm publish.

> **Rev 2 changes** (external review, 2026-09-20). The search becomes
> systematic instead of random, and the PRNG disappears with it (§3.2); the
> ring test moves from atoms to bonds (§3.1); the input's own minimum is
> always kept (§3.2); the butane numbers are measured and cross-checked
> instead of asserted (§5); `rmsd.ts` gains a reflection gate (§5); the API
> takes a plain `Molecule` and returns the full energy breakdown (§2);
> "conformer" and "duplicate" are defined (§1); the timing claim is grounded
> in `docs/benchmark.md` (§6).

## 1. Scope

**In scope:** from one 3D molecule, produce a diverse set of minimized
conformers sorted by energy — a systematic torsional search with energy and
RMSD pruning (Confab's genre), the minimizer relaxing each start. Note
"systematic": the original draft called this Confab-style but sampled
randomly (§3.2).

**Out of scope (deliberate):**

- *No de novo 3D embedding.* mmff94-ts takes coordinates as input (SDF); it
  has no embedder and this feature doesn't add one. The graph-walk embedder
  lives in Valence. Conformer generation here is **perturb-and-minimize**,
  not ETKDG.
- *Ring conformations are not sampled.* Ring bonds are never driven, so
  chair/twist-boat interconversion will not be found. Documented, not
  silent. Cheapest honest workaround today: hand the generator a structure
  that already carries the pucker you want — one call per ring family —
  which `docs/conformers.md` will say explicitly. A `seed_geometries`
  option only if users ask for one.
- *Dedup is not symmetry-aware.* v1 compares heavy-atom Kabsch RMSD without
  graph automorphism. Read the vocabulary below before treating that as a
  defect: it is narrower than it sounds.

**Vocabulary** (decided here, because the symmetry work depends on it). A
*conformer* is a distinct local minimum on the MMFF94 surface as the
coordinates see it. Mirror-related wells — butane's gauche+ and gauche− —
are **two conformers, not a duplicate**: they are two wells, related by
reflection, and no rigid superposition maps one onto the other. A
*duplicate* is the **same** geometry reached by relabeling equivalent
atoms, and v1 does not detect those. That is the real limitation, and any
future symmetry-aware RMSD must preserve the mirror pair rather than
collapse it.

The returned set is a set of **minima, not a Boltzmann population** — no
entropy, no populations, gas-phase MMFF94 energies. `docs/conformers.md`
must say so, or someone will compute populations from it.

## 2. API

```ts
export interface ConformerOptions {
  count?: number;               // target conformers; may return fewer (default 50)
  rmsd_threshold?: number;      // Å, heavy-atom Kabsch RMSD dedup cutoff (default 0.5)
  energy_window?: number;       // kcal/mol above the best minimum to keep (default 10)
  max_rotatable_bonds?: number; // sanity cap; warn and proceed (default 20)
  optimize_tolerance?: number;  // forwarded to the optimizer (default: its own criterion)
}

export interface Conformer {
  molecule: TypedMolecule;      // full molecule, coordinates replaced
  energy: EnergyComponents;     // the full per-term breakdown, not just the total
  rotor_dihedrals: number[];    // degrees, one per rotatable bond
}

export function generate_conformers(
  mol: Molecule | TypedMolecule,
  options?: ConformerOptions,
): Conformer[];
```

- **The input is prepared the way `calc_energy()` and `optimize_lbfgs()`
  prepare it** (types, then BCI charges), so an untyped `Molecule` works
  here as it does everywhere else in this library. The original draft
  required a charged `TypedMolecule` to avoid a silent fallback; that
  concern is answered now by `diagnose_molecule()`, which reports what the
  field could not represent through the warning handler instead of hiding
  it. Same honesty, no friction, one consistent entry contract.
- **The diagnostics run once on the input, not per conformer.** Atom types
  and BCI charges are graph properties — conformation cannot change them —
  so the search itself is pure geometry.
- **No `seed` option.** The enumeration is deterministic by construction
  (§3.2), so reproducible runs need no knob.
- **`energy` is the full breakdown.** Every other surface in this library
  exposes the seven terms; returning only a total would be the one place a
  caller cannot see where an energy came from.
- **`rotor_dihedrals`** is the vector the cheap dedupe walks on (§3.2); it
  makes a returned conformer inspectable without re-deriving anything.
  `count` is a target, not a promise: the returned set may be smaller, and
  is documented as such.

## 3. Pipeline

Ordered for cost: everything cheap runs before the first force-field call.

1. **Find rotatable bonds** (`src/conformers/rotatable.ts`). A bond is
   rotatable iff:
   - bond order is 1;
   - it is not a ring bond — a **bond-level** question, not an atom-level
     one. `find_ring_atoms()` answers "is this atom in a ring", and two ring
     atoms can share a bond that is in no ring: biphenyl's central bond is
     the example, and that rotor is exactly the one worth driving. The test
     is a bridge pass over the bond graph (a bond is a ring bond iff it is
     not a bridge); 3- and 4-membered ring bonds fall out of the same rule;
   - both atoms carry at least one *other* heavy-atom neighbor (rotating a
     methyl or a terminal halogen moves nothing the heavy-atom RMSD sees);
   - it is not in the partial-double exclusion list (amide, thioamide,
     conjugated ester C–O — a named constant, easy to extend).
2. **Enumerate starts — deterministic and systematic.** The staggered grid
   is {60°, 180°, 300°}^k for k rotatable bonds. If the grid fits the
   budget, every point is used; if not, a **Halton-spread subset** of budget
   points is (low-discrepancy over the rotor torus: deterministic, well
   spread, and no PRNG — which also drops the single-user module the first
   draft needed). Budget = `max(4 * count, 200)`, so the cheap rejects in
   step 3 cannot starve the target. The **unperturbed input, minimized, is
   always start #1**: the user's own geometry is a meaningful sample and the
   natural reference for `energy_window`.
3. **Cheap rejects, before any force-field call.** Drop a start whose
   rounded rotor vector (15° bins) duplicates one already accepted —
   torsion-space dedupe, no RMSD involved — and drop any start with a
   heavy-atom pair closer than a fraction of the vdW minimum (a clash
   screen). This is what keeps `count` affordable (§6), and it must not
   collapse mirror-related wells: gauche+ and gauche− differ in the sign of
   the rotor, so a bin-distance test keeps both (§5).
4. **Set the torsions.** Measure each rotor with `dihedral_angle()`, collect
   the atoms on one side by BFS, rotate them about the j→k bond *line* with
   `rotate_around_axis()`. Two details the first draft glossed over:
   - `Vec3` is a `[x, y, z]` tuple, not an atom object, so every geometry
     call needs the conversion the energy terms already do (`torsion.ts`).
     `tsc` catches the mismatch; a tsx probe does not — it silently produced
     NaN coordinates while I measured the butane gap for this review.
   - `rotate_around_axis()` rotates about an axis through the *origin*, so
     the usable form is translate–rotate–translate about the pivot atom.
   Assert every coordinate is finite after perturbation.
5. **Minimize** the survivors with `optimize_lbfgs()`; record the full
   energy components.
6. **Filter + dedupe.** Find the best minimum; drop everything more than
   `energy_window` above it; sort by energy; greedily keep a conformer only
   if its heavy-atom Kabsch RMSD to every kept conformer exceeds
   `rmsd_threshold`. The returned set may be smaller than `count`.

## 4. New code

| File | Contents |
|---|---|
| `src/conformers/rotatable.ts` | `find_rotatable_bonds()` + the bond-level ring (bridge) test |
| `src/conformers/rmsd.ts` | Kabsch alignment + heavy-atom RMSD, reflection-guarded |
| `src/conformers/generate.ts` | `generate_conformers()`: enumeration, cheap rejects, orchestration |
| `src/index.ts` | export the new API |
| `tests/conformers.test.ts` | §5 |
| `tests/references/conformers/` | recorded cross-check ensembles (§5) |
| `tests/scripts/conformer-cross-check.py` | the RDKit script that writes them |

There is no `prng.ts` in rev 2: the enumeration is deterministic, so the
seeded generator had exactly one user and no longer exists.

## 5. Validation

**Butane — measured, not asserted.** Driving the C1–C2–C3–C4 torsion with
this library's own primitives and minimizing (2026-09-20):

| minimum | τ | mmff94-ts | RDKit MMFF94 (60 ETKDG conformers) |
|---|---|---|---|
| anti | 180° | −5.07596 | −5.0760 |
| gauche | ∓64.9° | −4.29265 | −4.2938 (τ ∓65.3°) |
| gap | | **0.783** | **0.782** |

Three wells by the §1 vocabulary (anti, gauche+, gauche−); two up to mirror
symmetry, which is what a symmetry-aware grouping reports. The original
draft asserted "the gauche pair ≈ 0.9 kcal/mol (the MMFF94 value)" — wrong
by ~15% and unsourced, which AGENTS.md rule 3 forbids. The measured pair
above is what the test asserts (gap matched within 0.02). Note also that
the gauche well sits at ~65°, not 60°: starting a rotor at 60° is fine
because the minimizer relaxes it, but nobody should read 60° as the gauche
minimum.

The RDKit column comes from `tests/scripts/conformer-cross-check.py` (run on
lenovo, where RDKit lives), whose output is recorded under
`tests/references/conformers/`. Tests assert against the recording, so CI
needs no RDKit — the pattern the OpenBabel and Tinker reference logs already
use. This is the numerical cross-check AGENTS.md permits, and the first
draft had no equivalent.

Other gates:

- **`rmsd.ts`, directly:** zero for identical coordinates, invariant under
  rigid motion, and **non-zero for a mirror image** — the Kabsch reflection
  guard (reverse the sign of the smallest singular value when det < 0). That
  is the most likely silent bug in the feature and the first draft never
  tested the primitive.
- **Every returned conformer is a stationary point:** max|grad| within the
  optimizer's own criterion, and all coordinates finite.
- **The cheap dedupe keeps mirror pairs** (gauche+ and gauche− both survive
  §3.3) — the guard against a dedupe that is too eager.
- **Determinism:** two identical calls give identical output (tolerance
  bands per AGENTS.md if platform drift ever appears).
- **Energy window:** everything returned is within the window of the best.
- **Degenerate input** (benzene, ethane): one minimized conformer.
- **Over the cap:** `max_rotatable_bonds` exceeded → the warning handler
  fires and the call proceeds.
- **Pentane** (sourced per AGENTS.md — never hand-built): its staggered
  minima, recorded when measured, not predicted here.

## 6. Performance

Grounded in the repo's own benchmark (`docs/benchmark.md`: 29 drug-like
molecules, 25–70 atoms): **242 ms median, 446 ms mean per optimization**.
Fifty conformers is therefore *tens of seconds*, not seconds — "seconds in
Node" was optimistic by roughly an order of magnitude. What keeps the real
cost down is §3.3: the force field only ever sees starts that survived the
torsion-space dedupe and the clash screen. Measure the actual number when
the feature lands and record it; do not restate an estimate.

## 7. Docs (ship with the feature)

- README: "Conformers" section with a quickstart.
- New `docs/conformers.md`: algorithm, defaults, and the §1 limitations —
  including, explicitly, that the result is a set of **minima, not
  populations**, that rings are not sampled (with the workaround), and what
  v1 does not dedupe.
- CHANGELOG.md (doesn't exist yet — create it) with the feature entry.

## 8. 1.0 checklist (after this lands)

- `package.json` version 0.1.0-alpha.1 → 1.0.0
- CHANGELOG.md 1.0.0 entry
- README quickstart covers conformers
- `npm publish` — **needs Billy's npm credentials; his step**

## 9. Open questions

- ~~Default `count` scaling with rotor count~~ — settled by the budget rule
  (§3.2): the grid size decides, `count` is the target the search fills.
- ~~Charges in, or assigned internally~~ — settled (§2): prepare internally,
  as every other entry point does.
- **Symmetry-aware RMSD (v2):** a canonical atom ranking catches the real
  duplicates, and it must be validated to *keep* mirror pairs (§1).
- **`seed_geometries` (ring puckers):** only if users ask; document the
  workaround first.
- **Boltzmann weighting / populations:** out of scope by design — state it
  in the docs rather than leaving it implied.
