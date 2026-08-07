/**
 * Probe: .mmd reference pchg (old vs new suite) and ours vs the new pchg
 * for the revision-affected molecules. The .mmd pchg column is what
 * BatchMin actually used for the electrostatic term.
 *
 * Run: npx tsx tests/scripts/probe-pchg.ts
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { assign_bci_charges } from '../../src/mmff94/charges';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const NEW_MMD = path.join(HERE, '..', 'fixtures', 'validation-suite', 'MMFF94.mmd');
const OLD_MMD = 'C:/Users/mccan/AppData/Local/Temp/mmff94-753-backup/MMFF94.mmd';
const CODES = ['CEWYIM30','DAKCEX','FAPLUD','GIGCEE','KEPKIZ','SAKGUG','TAPJUP','VEWZOM'];

const newMols = new Map(parse_mmd(fs.readFileSync(NEW_MMD, 'utf8')).map(m => [m.name, m]));
const oldMols = new Map(parse_mmd(fs.readFileSync(OLD_MMD, 'utf8')).map(m => [m.name, m]));

for (const code of CODES) {
  const oldP = oldMols.get(code)!.atoms.map(a => a.partial_charge ?? NaN);
  const newP = newMols.get(code)!.atoms.map(a => a.partial_charge ?? NaN);
  const moved = Math.max(...newP.map((p, i) => Math.abs(p - oldP[i])));
  const typed = assign_atom_types(newMols.get(code)!);
  const ch = assign_bci_charges(typed);
  const ours = ch.atoms.map(a => a.partial_charge ?? NaN);
  const d = Math.max(...ours.map((o, i) => Math.abs(o - newP[i])));
  console.log(`${code}: ref pchg moved max ${moved.toFixed(4)}, ours vs new pchg max ${d.toFixed(4)}`);
}
