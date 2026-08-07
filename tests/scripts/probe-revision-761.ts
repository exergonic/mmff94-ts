/**
 * Focused probe: the 16 molecules affected by the suite revision
 * (8 corrected members + 8 ERULE additions). Prints, per term:
 *   old reference (May 1998 log) | new reference (Nov 1998 log) | ours on the new geometry
 * plus our delta vs the new reference.
 *
 * Run: npx tsx tests/scripts/probe-revision-761.ts
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_bmin_log } from './bmin-log';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { calc_energy } from '../../src/mmff94/energy/total';

const SUITE = path.join(path.dirname(fileURLToPath(import.meta.url)), '..', 'fixtures', 'validation-suite');

const CODES = [
  // corrected members (old reference values superseded)
  'CEWYIM30', 'DAKCEX', 'FAPLUD', 'GIGCEE', 'KEPKIZ', 'SAKGUG', 'TAPJUP', 'VEWZOM',
  // new molecules
  'ERULE_01', 'ERULE_02', 'ERULE_03', 'ERULE_04', 'ERULE_05', 'ERULE_06', 'ERULE_07', 'ERULE_08',
];

const newRef = parse_bmin_log(fs.readFileSync(path.join(SUITE, 'MMFF94_bmin.log'), 'utf8'));
const oldRef = parse_bmin_log(
  fs.readFileSync('C:/Users/mccan/AppData/Local/Temp/mmff94-753-backup/MMFF94_bmin.log', 'utf8'));

const mmdText = fs.readFileSync(path.join(SUITE, 'MMFF94.mmd'), 'utf8');
const molecules = new Map(parse_mmd(mmdText).map((m) => [m.name ?? '', m]));

const TERMS = ['stretch', 'bend', 'strbnd', 'torsion', 'oop', 'elec', 'vdw'] as const;
const OURS_KEY: Record<string, keyof typeof e> = {
  stretch: 'bond_stretch', bend: 'angle_bend', strbnd: 'stretch_bend',
  torsion: 'torsion', oop: 'out_of_plane', elec: 'electrostatic', vdw: 'van_der_waals',
};

for (const code of CODES) {
  const mol = molecules.get(code);
  if (!mol) { console.log(code, 'NOT IN MMD'); continue; }
  const typed = assign_atom_types(mol);
  const e = calc_energy(typed);
  const rNew = newRef.get(code);
  const rOld = oldRef.get(code);
  console.log(`\n=== ${code} ===`);
  for (const t of TERMS) {
    const ours = e[OURS_KEY[t]];
    const n = rNew ? rNew[t] : NaN;
    const o = rOld ? rOld[t] : NaN;
    const dn = (isFinite(n) && isFinite(ours)) ? ours - n : NaN;
    const do_ = (isFinite(o) && isFinite(ours)) ? ours - o : NaN;
    const moved = Math.abs(n - o) > 1e-4 ? `  ref moved ${(n - o).toFixed(4)}` : '';
    console.log(`  ${t.padEnd(15)} old ${String(o.toFixed(5)).padStart(11)}  new ${String(n.toFixed(5)).padStart(11)}  ours ${String(ours.toFixed(5)).padStart(11)}  d(new) ${String(dn.toFixed(5)).padStart(11)}${moved}`);
  }
  if (rNew && rOld) {
    const moved = Math.abs(rNew.total - rOld.total) > 1e-4;
    console.log(`  TOTAL            old ${rOld.total.toFixed(5).padStart(11)}  new ${rNew.total.toFixed(5).padStart(11)}  ours ${e.total.toFixed(5).padStart(11)}  d(new) ${(e.total - rNew.total).toFixed(5).padStart(11)}${moved ? '  <<< REF MOVED' : ''}`);
  }
}
