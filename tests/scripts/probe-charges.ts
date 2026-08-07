/**
 * Probe: our BCI partial charges vs the reference's (opti log # char
 * section) for the 16 revision-affected molecules.
 *
 * Run: npx tsx tests/scripts/probe-charges.ts
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { assign_bci_charges } from '../../src/mmff94/charges';

const SUITE = path.join(path.dirname(fileURLToPath(import.meta.url)), '..', 'fixtures', 'validation-suite');
const OPTI = process.env.OPTI_LOG || 'C:/Users/mccan/AppData/Local/Temp/mmff94-761/MMFF94_opti.log';
const CODES = ['CEWYIM30','DAKCEX','FAPLUD','GIGCEE','KEPKIZ','SAKGUG','TAPJUP','VEWZOM',
  'ERULE_01','ERULE_02','ERULE_03','ERULE_04','ERULE_05','ERULE_06','ERULE_07','ERULE_08'];

const txt = fs.readFileSync(path.join(SUITE, 'MMFF94.mmd'), 'utf8');
const mols = new Map(parse_mmd(txt).map(m => [m.name, m]));
const opti = fs.readFileSync(OPTI, 'utf8');

for (const code of CODES) {
  const m = mols.get(code);
  if (!m) { console.log(code, 'not in mmd'); continue; }
  const typed = assign_atom_types(m);
  const ch = assign_bci_charges(typed);
  const ours = ch.atoms.map(a => a.partial_charge ?? NaN);
  // reference charges from the opti log # char section
  let i = opti.indexOf('New Structure Name/Conformational Index: ' + code);
  if (i < 0) i = opti.indexOf('Structure Name: ' + code);
  const sec = opti.slice(i, i + 20000);
  const j = sec.indexOf('# char');
  const k = sec.indexOf('# fchar');
  const block = sec.slice(j, k);
  const refs = [...block.matchAll(/#(\d+)\s+([-0-9.]+)/g)]
    .sort((a, b) => Number(a[1]) - Number(b[1]))
    .map(mm => Number(mm[2]));
  if (refs.length !== ours.length) { console.log(code, 'length mismatch', refs.length, ours.length); continue; }
  const maxd = Math.max(...ours.map((o, idx) => Math.abs(o - refs[idx])));
  console.log(code, 'max |d(q)| =', maxd.toFixed(4),
    '  ours:', ours.map(o => o.toFixed(3)).join(' '));
  if (maxd > 1e-3) {
    console.log('        ref :', refs.map(r => r.toFixed(3)).join(' '));
  }
}
