/**
 * Probe: our types vs the opti log's own # ty assignments for the
 * 16 revision-affected molecules (8 ERULE + 8 corrected).
 * Run: npx tsx tests/scripts/probe-opti-typing-761.ts
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SUITE = path.join(HERE, '..', 'fixtures', 'validation-suite', 'MMFF94.mmd');
const OPTI = 'C:/Users/mccan/AppData/Local/Temp/mmff94-761/MMFF94_opti.log';
const codes = ['ERULE_01', 'ERULE_02', 'ERULE_03', 'ERULE_04', 'ERULE_05', 'ERULE_06', 'ERULE_07', 'ERULE_08',
  'CEWYIM30', 'DAKCEX', 'FAPLUD', 'GIGCEE', 'KEPKIZ', 'SAKGUG', 'TAPJUP', 'VEWZOM'];
const opti = fs.readFileSync(OPTI, 'utf8');
const mols = new Map(parse_mmd(fs.readFileSync(SUITE, 'utf8')).map(x => [x.name, x]));
for (const code of codes) {
  const m = mols.get(code)!;
  const typed = assign_atom_types(m);
  // the opti log's # ty block: lines " NAME #n  TYPE" after the
  // "OPTIMOL-LIST>  # ty" header within the molecule's section.
  const i = opti.indexOf(`New Structure Name/Conformational Index: ${code}`);
  const sec = opti.slice(i, i + 8000);
  const j = sec.indexOf('OPTIMOL-LIST>  # ty');
  const k = sec.indexOf('OPTIMOL-LIST>  # ', j + 20);
  const tyLines = sec.slice(j, k < 0 ? j + 3000 : k).split('\n');
  const optiTypes = new Map<number, number>();
  for (const line of tyLines) {
    const mm = line.match(/^\s*(\S+)\s+#(\d+)\s+(-?\d+)/);
    if (mm) optiTypes.set(Number(mm[2]), Number(mm[3]));
  }
  // match by atom index (the mmd order == the opti order for these)
  let bad = 0;
  for (let a = 0; a < m.atoms.length; a++) {
    const ot = optiTypes.get(a + 1);
    if (ot !== undefined && ot !== typed.atom_types[a]) {
      bad++;
      console.log(`  ${code} atom ${a} (${m.atoms[a].element}): opti ${ot} vs ours ${typed.atom_types[a]}`);
    }
  }
  console.log(`${code}: ${m.atoms.length} atoms, ${bad} divergences`);
}
