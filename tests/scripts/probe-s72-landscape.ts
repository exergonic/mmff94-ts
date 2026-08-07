/**
 * Survey: every S(72) atom in the suite — environment vs ours-vs-ref
 * charge, to derive the q⁰(72) environment rule.
 * Run: npx tsx tests/scripts/probe-s72-landscape.ts
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { assign_bci_charges } from '../../src/mmff94/charges';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SUITE = path.join(HERE, '..', 'fixtures', 'validation-suite', 'MMFF94.mmd');
const counts: Record<string, number> = {};
const mols = parse_mmd(fs.readFileSync(SUITE, 'utf8'));
for (const m of mols) {
  const typed = assign_atom_types(m);
  const charged = assign_bci_charges(typed);
  const adj: number[][] = Array.from({ length: m.atoms.length }, () => []);
  for (const b of m.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
  for (let i = 0; i < m.atoms.length; i++) {
    if (typed.atom_types[i] !== 72) continue;
    const ref = m.atoms[i].partial_charge ?? NaN;
    const ours = charged.partial_charges[i];
    const delta = Math.abs(ours - ref);
    // environment: the S's neighbors (with the S–P bond order)
    const env = adj[i].map(k => {
      const nb = m.atoms[k];
      const nbT = typed.atom_types[k];
      const order = m.bonds.find(b => (b.atom1 === i && b.atom2 === k) || (b.atom1 === k && b.atom2 === i))?.bond_order ?? 0;
      // does the neighbor carry an O2CM (type 32) or another terminal S?
      const o2cm = adj[k].some(b => typed.atom_types[b] === 32);
      const s2cm = adj[k].some(b => b !== i && typed.atom_types[b] === 72);
      return `${nb.element}${nbT}${order > 1 ? '=' : '-'}${o2cm ? '+O2CM' : ''}${s2cm ? '+S72' : ''}`;
    }).join(',');
    if (delta < 1e-3) counts[env] = (counts[env] ?? 0) + 1;
    const flag = delta >= 1e-3 ? '  <-- MISMATCH' : '';
    console.log(`${m.name.padEnd(8)} atom ${String(m.atoms[i].label ?? i).padEnd(4)}  ours ${ours.toFixed(4)}  ref ${ref.toFixed(4)}  Δ ${delta.toFixed(4)}  env[${env}]${flag}`);
  }
}
console.log('--- matching environments (count) ---');
for (const [env, n] of Object.entries(counts).sort((a, b) => b[1] - a[1])) console.log(`  ${env}: ${n}`);
console.log('done');
