// TAJSUS per-dihedral torsion diff: the library's resolution vs the
// reference's HIGH-verbosity torsional log (uv run ob_energy_breakdown.py
// TAJSUS --verbose torsion). The reference rows are transcribed below;
// their energy sum must reproduce the log's TOTAL (7.48219).
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';
import { calc_energy } from '../../src/mmff94/energy/total';
import { torsion_terms } from '../../src/mmff94/energy/torsion';
import { make_class_context } from '../../src/mmff94/parameters/parameter-classes';
import { dihedral_angle } from '../../src/utils/vector';

const name = process.argv[2] ?? 'TAJSUS';
const text = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');
const mol = parse_mmd(text).find(m => m.name === name)!;
const typed = assign_atom_types(mol);
const comp = calc_energy(typed);

const adj: number[][] = Array.from({ length: mol.atoms.length }, () => []);
for (const bond of mol.bonds) {
  adj[bond.atom1].push(bond.atom2);
  adj[bond.atom2].push(bond.atom1);
}
const pos = (i: number) => [mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z] as [number, number, number];
const ctx = make_class_context(typed, adj);

const DEG = 180 / Math.PI;
// Halgren convention: gamma1 = 0, gamma2 = 180, gamma3 = 0 (degrees).
const torsion_e = (tau: number, v1: number, v2: number, v3: number) =>
  (v1 / 2) * (1 + Math.cos(tau - 0)) +
  (v2 / 2) * (1 + Math.cos(2 * tau - Math.PI)) +
  (v3 / 2) * (1 + Math.cos(3 * tau - 0));

// Reference rows: [tI, tJ, tK, tL, class, angleDeg, v1, v2, v3, energy]
const REF: [number, number, number, number, number, number, number, number, number, number][] = [
  [32, 41, 1, 5, 0, -39.484, 0.000, 0.000, -0.106, -0.028],
  [32, 41, 1, 5, 0, 77.633, 0.000, 0.000, -0.106, -0.021],
  [32, 41, 1, 81, 0, -162.029, 0.000, 0.600, 0.000, 0.057],
  [32, 41, 1, 5, 0, 141.070, 0.000, 0.000, -0.106, -0.077],
  [32, 41, 1, 5, 0, -101.813, 0.000, 0.000, -0.106, -0.084],
  [32, 41, 1, 81, 0, 18.525, 0.000, 0.600, 0.000, 0.061],
  [41, 1, 81, 80, 0, -46.989, 0.000, 0.000, 0.000, 0.000],
  [41, 1, 81, 79, 0, 131.211, 0.000, 0.000, 0.000, 0.000],
  [5, 1, 81, 80, 0, -169.656, 0.000, 0.000, 0.000, 0.000],
  [5, 1, 81, 79, 0, 8.544, 0.000, 0.000, 0.000, 0.000],
  [5, 1, 81, 80, 0, 73.976, 0.000, 0.000, 0.000, 0.000],
  [5, 1, 81, 79, 0, -107.824, 0.000, 0.000, 0.000, 0.000],
  [5, 78, 81, 80, 0, -178.261, 0.000, 4.000, 0.000, 0.004],
  [5, 78, 81, 9, 0, -2.051, 0.000, 4.000, 0.000, 0.005],
  [79, 78, 81, 80, 0, 3.238, 0.000, 4.000, 0.000, 0.013],
  [79, 78, 81, 9, 0, 179.448, 0.000, 4.000, 0.000, 0.000],
  [5, 78, 79, 81, 0, 179.929, 0.000, 6.000, 0.000, 0.000],
  [81, 78, 79, 81, 0, -1.754, 0.000, 6.000, 0.000, 0.006],
  [5, 80, 81, 78, 0, 176.511, 0.000, 4.000, 0.000, 0.015],
  [5, 80, 81, 9, 0, 0.307, 0.000, 4.000, 0.000, 0.000],
  [81, 80, 81, 78, 0, -3.155, 0.000, 4.000, 0.000, 0.012],
  [81, 80, 81, 9, 0, -179.359, 0.000, 4.000, 0.000, 0.001],
  [5, 80, 81, 1, 0, 0.955, 0.000, 4.000, 0.000, 0.001],
  [5, 80, 81, 79, 0, -177.368, 0.000, 4.000, 0.000, 0.008],
  [81, 80, 81, 1, 0, -179.378, 0.000, 4.000, 0.000, 0.000],
  [81, 80, 81, 79, 0, 2.299, 0.000, 4.000, 0.000, 0.006],
  [35, 3, 9, 81, 0, -4.720, 0.000, 16.000, 0.000, 0.108],
  [1, 3, 9, 81, 0, 177.269, 0.000, 16.000, 0.000, 0.036],
  [35, 3, 1, 5, 0, -6.928, 0.000, 0.400, 0.300, 0.296],
  [35, 3, 1, 5, 0, -122.805, 0.000, 0.400, 0.300, 0.581],
  [35, 3, 1, 37, 0, 115.489, 0.000, 0.400, 0.300, 0.622],
  [9, 3, 1, 5, 0, 171.400, 0.000, 0.400, 0.300, 0.024],
  [9, 3, 1, 5, 0, 55.523, 0.000, 0.400, 0.300, 0.276],
  [9, 3, 1, 37, 0, -66.183, 0.000, 0.400, 0.300, 0.343],
  [3, 1, 37, 37, 0, -58.384, 0.000, 0.000, 0.200, 0.000],
  [3, 1, 37, 37, 0, 122.900, 0.000, 0.000, 0.200, 0.199],
  [5, 1, 37, 37, 0, 65.476, 0.000, -0.420, 0.391, -0.340],
  [5, 1, 37, 37, 0, -113.240, 0.000, -0.420, 0.391, 0.024],
  [5, 1, 37, 37, 0, -178.485, 0.000, -0.420, 0.391, 0.000],
  [5, 1, 37, 37, 0, 2.799, 0.000, -0.420, 0.391, 0.388],
  [1, 37, 37, 37, 0, -179.417, 0.000, 7.000, 0.000, 0.001],
  [1, 37, 37, 5, 0, 1.467, 0.000, 7.000, 0.000, 0.005],
  [37, 37, 37, 37, 0, -0.677, 0.000, 7.000, 0.000, 0.001],
  [37, 37, 37, 5, 0, -179.794, 0.000, 7.000, 0.000, 0.000],
  [1, 37, 37, 5, 0, -1.247, 0.000, 7.000, 0.000, 0.003],
  [1, 37, 37, 37, 0, 179.270, 0.000, 7.000, 0.000, 0.001],
  [37, 37, 37, 5, 0, -179.982, 0.000, 7.000, 0.000, 0.000],
  [37, 37, 37, 37, 0, 0.535, 0.000, 7.000, 0.000, 0.001],
  [37, 37, 37, 5, 0, 179.470, 0.000, 7.000, 0.000, 0.001],
  [37, 37, 37, 37, 0, -0.150, 0.000, 7.000, 0.000, 0.000],
  [5, 37, 37, 5, 0, -0.016, 0.000, 7.000, 0.000, 0.000],
  [5, 37, 37, 37, 0, -179.636, 0.000, 7.000, 0.000, 0.000],
  [37, 37, 37, 5, 0, 179.469, 0.000, 7.000, 0.000, 0.001],
  [37, 37, 37, 37, 0, -0.103, 0.000, 7.000, 0.000, 0.000],
  [5, 37, 37, 5, 0, -0.150, 0.000, 7.000, 0.000, 0.000],
  [5, 37, 37, 37, 0, -179.722, 0.000, 7.000, 0.000, 0.000],
  [37, 37, 37, 5, 0, 179.404, 0.000, 7.000, 0.000, 0.001],
  [37, 37, 37, 37, 0, -0.038, 0.000, 7.000, 0.000, 0.000],
  [5, 37, 37, 5, 0, -0.169, 0.000, 7.000, 0.000, 0.000],
  [5, 37, 37, 37, 0, -179.611, 0.000, 7.000, 0.000, 0.000],
  [37, 37, 37, 37, 0, 0.434, 0.000, 7.000, 0.000, 0.000],
  [37, 37, 37, 5, 0, 179.551, 0.000, 7.000, 0.000, 0.000],
  [5, 37, 37, 37, 0, -179.012, 0.000, 7.000, 0.000, 0.002],
  [5, 37, 37, 5, 0, 0.104, 0.000, 7.000, 0.000, 0.000],
  [1, 81, 79, 78, 0, -178.859, 0.000, 6.000, 0.000, 0.002],
  [80, 81, 79, 78, 0, -0.387, 0.000, 6.000, 0.000, 0.000],
  [78, 81, 9, 3, 1, 136.536, 0.000, 4.800, 0.000, 2.271],
  [80, 81, 9, 3, 1, -48.027, 0.000, 4.800, 0.000, 2.653],
];
const refSum = REF.reduce((s, r) => s + r[9], 0);
console.log(`REF rows: ${REF.length}  sum of energies: ${refSum.toFixed(5)}  (log TOTAL: 7.48219)`);

// Model side: same enumeration as calc_torsion_energy.
type Row = { t: number[]; ang: number; v: number[]; e: number };
const rows: Row[] = [];
for (const bond of mol.bonds) {
  const j = bond.atom1;
  const k = bond.atom2;
  const i_nbrs = adj[j].filter(n => n !== k);
  const l_nbrs = adj[k].filter(n => n !== j);
  if (i_nbrs.length === 0 || l_nbrs.length === 0) continue;
  for (const i of i_nbrs) {
    for (const l of l_nbrs) {
      if (i === l) continue; // 3-ring closure
      const tau = dihedral_angle(pos(i), pos(j), pos(k), pos(l));
      const terms = torsion_terms(ctx, typed, i, j, k, l);
      if (!terms) continue;
      const e = torsion_e(tau, terms.v1, terms.v2, terms.v3);
      rows.push({
        t: [typed.atom_types[i], typed.atom_types[j], typed.atom_types[k], typed.atom_types[l]],
        ang: tau * DEG,
        v: [terms.v1, terms.v2, terms.v3],
        e,
      });
    }
  }
}
const modelSum = rows.reduce((s, r) => s + r.e, 0);
console.log(`MODEL rows: ${rows.length}  sum: ${modelSum.toFixed(5)}  (calc_energy torsion: ${comp.torsion.toFixed(5)})`);
console.log(`delta vs reference: ${(modelSum - refSum).toFixed(6)}`);

// Match each reference row to a model row by (type-quad + angle).
let maxD = 0;
const unmatched: Row[] = [...rows];
for (const r of REF) {
  const key = r.slice(0, 4).join('-');
  const idx = unmatched.findIndex(
    m => m.t.join('-') === key && Math.abs(m.ang - r[5]) < 0.01,
  );
  if (idx < 0) {
    console.log(`NO MATCH: ${key} ang=${r[5]}`);
    continue;
  }
  const m = unmatched.splice(idx, 1)[0];
  const d = m.e - r[9];
  maxD = Math.max(maxD, Math.abs(d));
  const vDiff =
    m.v[0] !== r[6] || m.v[1] !== r[7] || m.v[2] !== r[8] ? ' V-DIFF' : '';
  if (Math.abs(d) > 1e-5 || vDiff) {
    console.log(
      `DIFF ${key} ang=${m.ang.toFixed(3)} (ref ${r[5].toFixed(3)}) ` +
        `V=[${m.v.join(',')}] (ref [${r[6]},${r[7]},${r[8]}]) ` +
        `E=${m.e.toFixed(6)} (ref ${r[9].toFixed(6)}) d=${d.toExponential(2)}${vDiff}`,
    );
  }
}
if (unmatched.length) {
  console.log('MODEL rows with no reference row:');
  for (const m of unmatched) console.log(`  ${m.t.join('-')} ang=${m.ang.toFixed(3)} V=[${m.v.join(',')}] E=${m.e.toFixed(6)}`);
}
console.log(`max |row delta|: ${maxD.toExponential(2)}`);
