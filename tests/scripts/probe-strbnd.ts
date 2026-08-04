// Per-angle stretch-bend diagnosis: key, table entry, k used.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';
import { STRETCH_BEND_PARAMS, DEFAULT_STRETCH_BEND, ELEMENT_ROW, lookup_param, bond_parameters } from '../../src/mmff94/parameters';
import { make_class_context, strbnd_type, angle_parameters } from '../../src/mmff94/parameters/parameter-classes';
import { calc_energy } from '../../src/mmff94/energy/total';
import { distance, angle_in_radians } from '../../src/utils/vector';

const name = process.argv[2] ?? 'VIYPAU';
const text = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');
const mol = parse_mmd(text).find(m => m.name === name)!;
const typed = assign_atom_types(mol);
const comp = calc_energy(typed);

const adj: number[][] = Array.from({ length: mol.atoms.length }, () => []);
for (const bond of mol.bonds) {
  adj[bond.atom1].push(bond.atom2);
  adj[bond.atom2].push(bond.atom1);
}
const ctx = make_class_context(typed, adj.map(n => [...n]));

for (let j = 0; j < mol.atoms.length; j++) {
  if (adj[j].length < 2) continue;
  const tj = typed.atom_types[j];
  for (let a = 0; a < adj[j].length; a++) {
    for (let b = a + 1; b < adj[j].length; b++) {
      const i = adj[j][a];
      const k = adj[j][b];
      const ti = typed.atom_types[i];
      const tk = typed.atom_types[k];
      const t_min = Math.min(ti, tk);
      const t_max = Math.max(ti, tk);
      const cls = strbnd_type(ctx, i, j, k);
      let entry;
      if (cls !== 0) {
        for (const key of [`${cls}-${t_min}-${tj}-${t_max}`, `${cls}-0-${tj}-${t_max}`, `${cls}-${t_min}-${tj}-0`, `${cls}-0-${tj}-0`]) {
          if (STRETCH_BEND_PARAMS[key]) { entry = { key, ...STRETCH_BEND_PARAMS[key] }; break; }
        }
      } else {
        entry = lookup_param(STRETCH_BEND_PARAMS, [t_min, tj, t_max]);
      }
      let used;
      let F: [number, number] | undefined;
      const ra = ELEMENT_ROW[mol.atoms[i].element] ?? '?';
      const rb = ELEMENT_ROW[mol.atoms[j].element] ?? '?';
      const rc = ELEMENT_ROW[mol.atoms[k].element] ?? '?';
      if (entry) {
        used = `entry k=${entry.k_sb_IJK ?? entry.k_sb}`;
      } else {
        F = DEFAULT_STRETCH_BEND[`${ra}-${rb}-${rc}`] ?? DEFAULT_STRETCH_BEND[`${rc}-${rb}-${ra}`];
        used = `DEFAULT ${ra}-${rb}-${rc} F=${F ? F.join('/') : 'NONE'}`;
      }
      // the energy contribution
      const ang = angle_parameters(ctx, i, j, k);
      if (!ang) continue;
      const b_ij = bond_parameters(ctx, i, j);
      const b_kj = bond_parameters(ctx, k, j);
      if (!b_ij || !b_kj) continue;
      const r_ij = distance(
        [mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z],
        [mol.atoms[j].x, mol.atoms[j].y, mol.atoms[j].z],
      );
      const r_kj = distance(
        [mol.atoms[k].x, mol.atoms[k].y, mol.atoms[k].z],
        [mol.atoms[j].x, mol.atoms[j].y, mol.atoms[j].z],
      );
      const theta = angle_in_radians(
        [mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z],
        [mol.atoms[j].x, mol.atoms[j].y, mol.atoms[j].z],
        [mol.atoms[k].x, mol.atoms[k].y, mol.atoms[k].z],
      );
      let k_ij = 0;
      let k_kj = 0;
      if (entry) {
        k_ij = ti <= tk ? entry.k_sb_IJK : entry.k_sb_KJI;
        k_kj = ti <= tk ? entry.k_sb_KJI : entry.k_sb_IJK;
      } else if (F) {
        const direct = DEFAULT_STRETCH_BEND[`${ra}-${rb}-${rc}`] ? true : false;
        k_ij = direct ? F[0] : F[1];
        k_kj = direct ? F[1] : F[0];
      }
      if (k_ij === 0 && k_kj === 0) continue;
      const e = 2.51210 * (k_ij * (r_ij - b_ij.r0) + k_kj * (r_kj - b_kj.r0)) * (theta - ang.theta0 * Math.PI / 180);
      console.log(
        `angle ${i}-${j}-${k} (${mol.atoms[i].element}${ti}-${mol.atoms[j].element}${tj}-${mol.atoms[k].element}${tk}) cls=${cls} ${used} E=${e.toFixed(4)}`,
      );
    }
  }
}
console.log(`our strbnd total: ${comp.stretch_bend.toFixed(4)}`);

