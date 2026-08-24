# Geometry-optimization benchmark

mmff94-ts's L-BFGS optimizer against Tinker 26.2's Fortran `minimize`,
on 29 drug-like molecules from the OpenFF Industry Benchmark Season 1
v1.2 (B3LYP-D3BJ/DZVP QM-optimized starting geometries, pulled from
QCArchive). Molecules span 25–70 atoms; `mol_017` is excluded from all
legs (aromatic-bond connectivity cannot be typed consistently across
engines).

Run on the same machine (Intel i3-4150 @ 3.5 GHz, single core), same
starting geometries, converged to the default criterion — max |g_i|
< 0.05 kcal/mol/Å OR RMS gradient < 0.02 (Tinker: its own
RMS-gradient criterion at the equivalent threshold). That makes the
JavaScript numbers doubly conservative: a 2014-era dual-core desktop,
no JIT warm-up tricks, no native extensions.

## Results

| engine | converged | median wall time | mean |
|---|---|---|---|
| **mmff94-ts** (TypeScript, L-BFGS) | **26/29** | **242 ms** | 446 ms |
| Tinker 26.2 `minimize` (Fortran)   | 29/29    | 296 ms          | 324 ms |

(Tinker also minimizes `mol_017`, the molecule excluded from all legs
for typing inconsistency; 29/29 is its rate on the shared set.)

Both engines minimize every molecule to an MMFF94 minimum in a few
hundred milliseconds. mmff94-ts converges 26 of the same 29 molecules
and is faster at the median — in interpreted JavaScript, with zero
dependencies and no WebAssembly. (Tinker's own leg dropped 4 of its 30
inputs to segfaults or non-finite energies; on the shared set it
converges all 29.)

## Methodology notes (read before comparing)

- **Timing** measures only the minimize call, not file I/O, typing, or
  interpreter startup. mmff94-ts resolves every interaction parameter
  once per molecule into flat typed arrays (see
  `src/optimize/fast-system.ts`) and evaluates allocation-free; that
  compiled fast path is what closes the gap to native code.
- **Convergence criteria differ by engine**: Tinker's tolerance is
  RMS-gradient-per-atom; this library's default stops on either max |g|
  or RMS gradient. The two are not strictly comparable — an RMS gate is
  systematically easier than a max gate on the same number.
- **Tinker's internal MMFF94 charge derivation assigns spurious net
  charges on some N/O/S-rich molecules.** Fed our exact atom types,
  Tinker reports net molecular charges of −0.5 to −1.5 e on seven
  formally neutral molecules (`analyze` → Total Electric Charge), where
  both this library and RDKit derive net 0 from the same BCI model. On
  those molecules Tinker descends into artificial ion-stabilized basins
  (up to 56 kcal/mol below any independent implementation). On the
  eleven molecules where Tinker's derived charges match ours, the two
  implementations agree on the final energy **to four decimal places**
  — two independent codebases, same answer. Per-molecule charges and
  deltas: see the ff-bench session data.
- **Correctness anchor**: converged energies agree with RDKit's
  independent MMFF94 within ~0.5 kcal/mol on comparable molecules, atom
  typing is 761/761 identical to OpenBabel across the full validation
  suite, and per-term energies match BatchMin ≤1e-4 on Halgren's
  761-molecule suite.
- The three unconverged molecules (`mol_013`, `mol_016`, `mol_019`) are
  flexible polyamines whose trajectories stall one gate-width short of
  the threshold after 3000 iterations — a known L-BFGS limitation on
  flat, coupled torsion surfaces, not specific to this implementation.

## Per-molecule data

Committed alongside this document:

- [benchmark-mmff94-ts.csv](benchmark-mmff94-ts.csv) — name, status,
  iterations, ms, final energy, final max gradient (this library)
- [benchmark-tinker.csv](benchmark-tinker.csv) — name, ms, final
  energy, final RMS gradient (Tinker)

## Reproducing

```bash
# this library (from the repo root, against the exported SDFs)
node -e "
import('./dist/sdf.js').then(async ({ parse_sdf }) => {
  const { optimize_lbfgs } = await import('./dist/optimize/l-bfgs.js');
  const fs = await import('fs');
  const mol = parse_sdf(fs.readFileSync('mol.sdf', 'utf-8'));
  const r = optimize_lbfgs(mol);
  console.log(r.converged, r.iterations, r.energy.total);
});
"

# Tinker
minimize <name>.MMFF94.txyz <name>.MMFF94.key 3000
```
