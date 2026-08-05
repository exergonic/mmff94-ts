// Drug-size benchmark: energy and gradient evaluation on trp-cage
// (1L2Y model 1, 304 atoms) — a real 20-residue peptide with an
// experimental NMR geometry (PDB, explicit hydrogens). The molecule
// is prepared ONCE (typing + charges); the timed calls are the
// per-step evaluations a consumer pays per optimization iteration.
//
// Run: npm run bench
import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import {
  parse_sdf,
  calc_energy,
  calc_gradient,
  assign_atom_types,
  assign_bci_charges,
} from '../src';

const molecule = parse_sdf(readFileSync(fileURLToPath(new URL('../tests/fixtures/sdf/trpcage.sdf', import.meta.url)), 'utf-8'));

const t_setup = performance.now();
const typed = assign_atom_types(molecule);
const charged = assign_bci_charges(typed);
const setup_ms = performance.now() - t_setup;

// One sanity pass before timing: the energy and gradient must be
// finite on the (experimental) geometry, and a sane force field
// energy is small in magnitude per atom.
const energy = calc_energy(charged);
const gradient = calc_gradient(charged);
if (!Number.isFinite(energy.total)) throw new Error(`non-finite energy: ${energy.total}`);
if (gradient.some((g) => !Number.isFinite(g[0]) || !Number.isFinite(g[1]) || !Number.isFinite(g[2]))) {
  throw new Error('non-finite gradient');
}
if (Math.abs(energy.total) > 1e5) throw new Error(`suspicious energy magnitude: ${energy.total}`);

/** Mean per-call time in ms over n iterations (3 warmup calls first). */
function time_it(fn: () => void, n: number): number {
  fn(); fn(); fn(); // warmup (JIT)
  const t0 = performance.now();
  for (let i = 0; i < n; i++) fn();
  return (performance.now() - t0) / n;
}

const N = 10;
const energy_ms = time_it(() => calc_energy(charged), N);
const gradient_ms = time_it(() => calc_gradient(charged), N);
const both_ms = time_it(() => {
  calc_energy(charged);
  calc_gradient(charged);
}, N);

console.log(`trp-cage 1L2Y: ${molecule.atoms.length} atoms, ${molecule.bonds.length} bonds`);
console.log(`setup (typing + charges): ${setup_ms.toFixed(1)} ms, E = ${energy.total.toFixed(2)} kcal/mol`);
console.log('');
console.log(`calc_energy             ${energy_ms.toFixed(1)} ms/call   ${(1000 / energy_ms).toFixed(1)} calls/s`);
console.log(`calc_gradient           ${gradient_ms.toFixed(1)} ms/call   ${(1000 / gradient_ms).toFixed(1)} calls/s`);
console.log(`energy + gradient       ${both_ms.toFixed(1)} ms/call   ${(1000 / both_ms).toFixed(1)} calls/s`);
