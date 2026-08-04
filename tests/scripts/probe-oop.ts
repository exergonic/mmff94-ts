// Per-center oop diagnosis mirroring the energy loop exactly.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';
import { calc_energy } from '../../src/mmff94/energy/total';
import { wilson_oop_angle } from '../../src/utils/vector';
import { OOP_PARAMS, lookup_param, ATOM_TYPE_PROPERTIES } from '../../src/mmff94/parameters';

const name = process.argv[2] ?? 'COYVIV';
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

// replicate the energy's oop_k (with the step-down chain)
function oop_k(sorted: number[], tj: number): number | undefined {
  let params = lookup_param(OOP_PARAMS, [sorted[0], tj, sorted[1], sorted[2]]);
  for (const lvl of ['lvl3', 'lvl4', 'lvl5'] as const) {
    if (params) break;
    for (let p = 0; p < 3; p++) {
      const t = [...sorted];
      t[p] = ATOM_TYPE_PROPERTIES[t[p]]?.[lvl] ?? t[p];
      t.sort((x, y) => x - y);
      params = lookup_param(OOP_PARAMS, [t[0], tj, t[1], t[2]]);
      if (params) break;
    }
  }
  if (!params) params = lookup_param(OOP_PARAMS, [0, tj, 0, 0]);
  return params?.k_oop;
}

let total = 0;
for (let j = 0; j < mol.atoms.length; j++) {
  if (adj[j].length !== 3) continue;
  const [a, c, d] = adj[j];
  const tj = typed.atom_types[j];
  const sorted = [typed.atom_types[a], typed.atom_types[c], typed.atom_types[d]].sort((x, y) => x - y);
  const k = oop_k(sorted, tj) ?? 0;
  const chi_a = wilson_oop_angle(pos(d), pos(j), pos(c), pos(a));
  const chi_c = wilson_oop_angle(pos(a), pos(j), pos(d), pos(c));
  const chi_d = wilson_oop_angle(pos(a), pos(j), pos(c), pos(d));
  const e = 0.043844 * (k / 2) * (chi_a * chi_a + chi_c * chi_c + chi_d * chi_d);
  total += e;
  console.log(
    `center ${j} t${tj} nbrs[${sorted}] k=${k.toFixed(4)} chis=[${chi_a.toFixed(1)},${chi_c.toFixed(1)},${chi_d.toFixed(1)}] E=${e.toFixed(4)}`,
  );
}
console.log(`sum: ${total.toFixed(4)}  calc_energy oop: ${comp.out_of_plane.toFixed(4)}`);
