/**
 * Probe: our BCI charges vs the reference .mmd pchg for one molecule.
 * The shallow-copy trap: Atom.partial_charge carries the PARSED pchg
 * through assign_bci_charges — always read the returned partial_charges
 * array, never Atom.partial_charge.
 * Run: npx tsx tests/scripts/probe-charges-one.ts <CODE>
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { assign_bci_charges } from '../../src/mmff94/charges';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const code = process.argv[2] ?? 'FAPLUD';
const SUITE = path.join(HERE, '..', 'fixtures', 'validation-suite', 'MMFF94.mmd');
const m = new Map(parse_mmd(fs.readFileSync(SUITE, 'utf8')).map(x => [x.name, x])).get(code)!;
const typed = assign_atom_types(m);
const charged = assign_bci_charges(typed);
console.log(`${code}  (ref pchg from the mmd column; ours from the partial_charges array)`);
let worst = 0;
for (let i = 0; i < m.atoms.length; i++) {
  const ref = m.atoms[i].partial_charge ?? NaN;
  const ours = charged.partial_charges[i];
  const d = Math.abs(ours - ref);
  worst = Math.max(worst, d);
  const flag = d > 1e-3 ? '  <--' : '';
  console.log(`  ${m.atoms[i].label ?? `${m.atoms[i].element}${i}`.padEnd(4)}  type ${String(typed.atom_types[i]).padStart(3)}  ours ${ours.toFixed(4)}  ref ${ref.toFixed(4)}  Δ ${d.toExponential(1)}${flag}`);
}
console.log('worst Δ:', worst.toExponential(3));
