/**
 * Probe: our resolved torsion rows for a molecule vs the opti log's.
 *
 * Run: npx tsx tests/scripts/probe-torsions.ts <CODE>
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { make_class_context, torsion_class, lookup_torsion, get_bond_order, is_aromatic_bond } from '../../src/mmff94/parameters/parameter-classes.js';
import { ATOM_TYPE_PROPERTIES } from '../../src/mmff94/parameters/index.js';
import { empirical_torsion } from '../../src/mmff94/parameters/empirical.js';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const code = process.argv[2] ?? 'ERULE_03';
const SUITE = path.join(HERE, '..', 'fixtures', 'validation-suite', 'MMFF94.mmd');
const m = new Map(parse_mmd(fs.readFileSync(SUITE, 'utf8')).map(x => [x.name, x])).get(code)!;
const typed = assign_atom_types(m);
const adj: number[][] = Array.from({ length: m.atoms.length }, () => []);
for (const b of m.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
const ctx = make_class_context(typed, adj);
const el = (i: number) => `${m.atoms[i].element}${i}`;

// enumerate all 4-paths
for (let j = 0; j < m.atoms.length; j++) {
  for (const i of adj[j]) {
    for (const k of adj[j]) {
      if (k === i) continue;
      for (const l of adj[k]) {
        if (l === j) continue;
        const cls = torsion_class(ctx, i, j, k, l);
        const t = [typed.atom_types[i], typed.atom_types[j], typed.atom_types[k], typed.atom_types[l]];
        const params = lookup_torsion(cls, t[0], t[1], t[2], t[3]);
        if (params) {
          console.log(`${el(i)}-${el(j)}-${el(k)}-${el(l)}  ${t.join('-')}  cls ${cls}  STORED V1 ${params.v1} V2 ${params.v2} V3 ${params.v3}`);
        } else {
          const emp = empirical_torsion(
            ATOM_TYPE_PROPERTIES[t[1]], ATOM_TYPE_PROPERTIES[t[2]],
            m.atoms[j].element, m.atoms[k].element,
            get_bond_order(ctx, j, k), is_aromatic_bond(ctx, j, k));
          if (emp.skip) {
            console.log(`${el(i)}-${el(j)}-${el(k)}-${el(l)}  ${t.join('-')}  cls ${cls}  SKIP`);
          } else {
            console.log(`${el(i)}-${el(j)}-${el(k)}-${el(l)}  ${t.join('-')}  cls ${cls}  EMPIRICAL V1 ${emp.v1} V2 ${emp.v2} V3 ${emp.v3}`);
          }
        }
      }
    }
  }
}
