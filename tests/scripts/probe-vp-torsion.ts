/**
 * Probe: our torsion resolution for the vinylphosphine fixture.
 * Run: npx tsx tests/scripts/probe-vp-torsion.ts
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_sdf } from '../../src/sdf.js';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types.js';
import { make_class_context, torsion_class, lookup_torsion, get_bond_order, is_aromatic_bond } from '../../src/mmff94/parameters/parameter-classes.js';
import { ATOM_TYPE_PROPERTIES } from '../../src/mmff94/parameters/index.js';
import { empirical_torsion } from '../../src/mmff94/parameters/empirical.js';
import { dihedral_angle } from '../../src/utils/vector.js';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const sdf = fs.readFileSync(path.join(HERE, '..', 'fixtures', 'sdf', 'vinylphosphine.sdf'), 'utf8');
const m = parse_sdf(sdf);
const typed = assign_atom_types(m);
const adj: number[][] = Array.from({ length: m.atoms.length }, () => []);
for (const b of m.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
const ctx = make_class_context(typed, adj);
const el = (i: number) => `${m.atoms[i].element}${i}`;
let total = 0;
for (let j = 0; j < m.atoms.length; j++) {
  for (const i of adj[j]) {
    for (const k of adj[j]) {
      if (k === i) continue;
      for (const l of adj[k]) {
        if (l === j || l === i) continue;
        const cls = torsion_class(ctx, i, j, k, l);
        const t = [typed.atom_types[i], typed.atom_types[j], typed.atom_types[k], typed.atom_types[l]];
        let v1 = 0, v2 = 0, v3 = 0, src = '?';
        const params = lookup_torsion(cls, t[0], t[1], t[2], t[3]);
        if (params) { v1 = params.v1; v2 = params.v2; v3 = params.v3; src = 'STORED'; }
        else {
          const emp = empirical_torsion(ATOM_TYPE_PROPERTIES[t[1]], ATOM_TYPE_PROPERTIES[t[2]],
            m.atoms[j].element, m.atoms[k].element, get_bond_order(ctx, j, k), is_aromatic_bond(ctx, j, k));
          if (!emp.skip) { v1 = emp.v1; v2 = emp.v2; v3 = emp.v3; src = 'EMP'; }
        }
        const tau = dihedral_angle(
          [m.atoms[i].x, m.atoms[i].y, m.atoms[i].z], [m.atoms[j].x, m.atoms[j].y, m.atoms[j].z],
          [m.atoms[k].x, m.atoms[k].y, m.atoms[k].z], [m.atoms[l].x, m.atoms[l].y, m.atoms[l].z]);
        const e = (v1 / 2) * (1 + Math.cos(1 * tau * Math.PI / 180))
          + (v2 / 2) * (1 + Math.cos(2 * tau * Math.PI / 180 - Math.PI))
          + (v3 / 2) * (1 + Math.cos(3 * tau * Math.PI / 180));
        total += e;
        console.log(`${el(i)}-${el(j)}-${el(k)}-${el(l)}  ${t.join('-')} cls ${cls} ${src} V1 ${v1} V2 ${v2} V3 ${v3}  τ ${tau.toFixed(1)} E ${e.toFixed(5)}`);
      }
    }
  }
}
console.log('torsion total:', total.toFixed(5));
