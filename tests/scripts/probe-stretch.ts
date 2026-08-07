/**
 * Probe: our per-bond stretch vs the opti log's for a molecule.
 * Run: npx tsx tests/scripts/probe-stretch.ts <CODE>
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { make_class_context, bond_parameters, get_bond_order } from '../../src/mmff94/parameters/parameter-classes.js';
import { empirical_bond_parameters } from '../../src/mmff94/parameters/empirical.js';
import { distance } from '../../src/utils/vector.js';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const code = process.argv[2] ?? 'ERULE_06';
const SUITE = path.join(HERE, '..', 'fixtures', 'validation-suite', 'MMFF94.mmd');
const OPTI = 'C:/Users/mccan/AppData/Local/Temp/mmff94-761/MMFF94_opti.log';
const m = new Map(parse_mmd(fs.readFileSync(SUITE, 'utf8')).map(x => [x.name, x])).get(code)!;
const typed = assign_atom_types(m);
const adj: number[][] = Array.from({ length: m.atoms.length }, () => []);
for (const b of m.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
const ctx = make_class_context(typed, adj);

// reference rows: name -> (ideal, k)
const opti = fs.readFileSync(OPTI, 'utf8');
const i = opti.indexOf('New Structure Name/Conformational Index: ' + code);
const sec = opti.slice(i, i + 40000);
const j = sec.indexOf('OPTIMOL-ANALYZE>  # bonds');
const k = sec.indexOf('OPTIMOL-ANALYZE>  # ', j + 10);
const ref = new Map<string, [number, number]>();
for (const line of sec.slice(j, k).split('\n')) {
  const mm = line.match(/^ (\S+) #\d+\s+(\S+) #\d+\s+\d+\s+\d+\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/);
  if (mm) ref.set(`${mm[1]}-${mm[2]}`, [Number(mm[4]), Number(mm[6])]);
}

for (const b of m.bonds) {
  const a1 = m.atoms[b.atom1], a2 = m.atoms[b.atom2];
  const r = distance([a1.x, a1.y, a1.z], [a2.x, a2.y, a2.z]);
  let p = bond_parameters(ctx, b.atom1, b.atom2);
  if (!p) p = empirical_bond_parameters(a1, a2);
  const label = `${a1.element}${b.atom1}-${a2.element}${b.atom2}`;
  const key = `${a1.element}${b.atom1 + 1}-${a2.element}${b.atom2 + 1}`;
  const rf = ref.get(key);
  const ourE = p ? 0.5 * 143.9325 * p.k_b * (r - p.r0) ** 2 * (1 - 2 * (r - p.r0) + (7 / 12) * 4 * (r - p.r0) ** 2) : 0;
  const refE = rf ? 0.5 * 143.9325 * rf[1] * (r - rf[0]) ** 2 * (1 - 2 * (r - rf[0]) + (7 / 12) * 4 * (r - rf[0]) ** 2) : NaN;
  console.log(`${label.padEnd(10)} r ${r.toFixed(4)}  our r0 ${p?.r0.toFixed(4)} k ${p?.k_b.toFixed(4)} E ${ourE.toFixed(6)}  | ref r0 ${rf ? rf[0] : NaN} k ${rf ? rf[1] : NaN} E ${refE.toFixed(6)}  ΔE ${(ourE - refE).toExponential(2)}`);
}
