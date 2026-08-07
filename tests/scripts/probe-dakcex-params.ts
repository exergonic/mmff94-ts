/**
 * Probe: our resolved bond/angle rows for DAKCEX vs the new opti log's.
 *
 * Run: npx tsx tests/scripts/probe-dakcex-params.ts
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { make_class_context, bond_parameters, angle_parameters } from '../../src/mmff94/parameters/parameter-classes.js';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const NEW_MMD = path.join(HERE, '..', 'fixtures', 'validation-suite', 'MMFF94.mmd');
const m = new Map(parse_mmd(fs.readFileSync(NEW_MMD, 'utf8')).map(x => [x.name, x])).get('DAKCEX')!;
const typed = assign_atom_types(m);
const adj: number[][] = Array.from({ length: m.atoms.length }, () => []);
for (const b of m.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
const ctx = make_class_context(typed, adj);

console.log('bond    types cls  our r0    our kb');
for (let i = 0; i < m.bonds.length; i++) {
  const b = m.bonds[i];
  const t1 = typed.atom_types[b.atom1], t2 = typed.atom_types[b.atom2];
  const p = bond_parameters(ctx, b.atom1, b.atom2);
  console.log(`${i}: ${b.atom1}-${b.atom2} ${String(t1).padStart(3)}-${String(t2).padStart(3)}  ${p ? `${p.r0.toFixed(3)} ${p.k_b.toFixed(3)}` : 'NO ROW'}`);
}
console.log('\nangle       types      our theta0 our ka');
for (let j = 0; j < m.atoms.length; j++) {
  if (adj[j].length < 2) continue;
  for (let a = 0; a < adj[j].length; a++) for (let b = a + 1; b < adj[j].length; b++) {
    const i = adj[j][a], k = adj[j][b];
    const p = angle_parameters(ctx, i, j, k);
    console.log(`${i}-${j}-${k} ${String(typed.atom_types[i]).padStart(3)}-${String(typed.atom_types[j]).padStart(3)}-${String(typed.atom_types[k]).padStart(3)}  ${p ? `${p.theta0.toFixed(3)} ${p.k_a.toFixed(4)}` : 'NO ROW'}`);
  }
}
