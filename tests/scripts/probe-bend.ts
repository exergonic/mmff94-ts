// Per-angle bend diagnosis: class, key, k, energy.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';
import { make_class_context, angle_class, angle_parameters } from '../../src/mmff94/parameters/parameter-classes';
import { calc_energy } from '../../src/mmff94/energy/total';
import { angle_in_radians } from '../../src/utils/vector';

const name = process.argv[2] ?? 'GESNIB';
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
      const cls = angle_class(ctx, i, j, k);
      const p = angle_parameters(ctx, i, j, k);
      const theta = angle_in_radians(
        [mol.atoms[j].x, mol.atoms[j].y, mol.atoms[j].z],
        [mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z],
        [mol.atoms[k].x, mol.atoms[k].y, mol.atoms[k].z],
      );
      if (!p) {
        console.log(`angle ${i}-${j}-${k} (${ti}-${tj}-${tk}) cls=${cls} NO PARAMS`);
        continue;
      }
      const theta0 = p.theta0 ?? 0;
      const delta = theta - (theta0 * Math.PI) / 180;
      const e = 0.043844 * (p.k_a / 2) * delta * delta;
      if (Math.abs(e) > 0.05 || Math.abs(delta) > 0.3) {
        console.log(
          `angle ${i}-${j}-${k} (${mol.atoms[i].element}${ti}-${mol.atoms[j].element}${tj}-${mol.atoms[k].element}${tk}) cls=${cls} ka=${p.k_a.toFixed(4)} th0=${theta0.toFixed(2)} deg=${(theta * 180 / Math.PI).toFixed(2)} E=${e.toFixed(4)}`,
        );
      }
    }
  }
}
console.log(`our bend total: ${comp.angle_bend.toFixed(4)}`);
