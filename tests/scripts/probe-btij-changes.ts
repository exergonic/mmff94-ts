/**
 * Probe: old-vs-new opti bond sections for the 7 corrected molecules.
 * Lists every bond whose FF class / ideal / force constant changed
 * between revisions, with types.
 *
 * Run: npx tsx tests/scripts/probe-btij-changes.ts
 */
import * as fs from 'fs';

const OLD = 'C:/Users/mccan/AppData/Local/Temp/mmff94-753-backup/MMFF94_opti.log';
const NEW = 'C:/Users/mccan/AppData/Local/Temp/mmff94-761/MMFF94_opti.log';
const CODES = ['CEWYIM30', 'DAKCEX', 'GIGCEE', 'KEPKIZ', 'SAKGUG', 'TAPJUP', 'VEWZOM'];

function bondSection(txt: string, code: string): string {
  let i = txt.indexOf('New Structure Name/Conformational Index: ' + code);
  if (i < 0) i = txt.indexOf('Structure Name: ' + code);
  const sec = txt.slice(i, i + 40000);
  const j = sec.indexOf('OPTIMOL-ANALYZE>  # bonds');
  const k = sec.indexOf('OPTIMOL-ANALYZE>  # ', j + 10);
  return sec.slice(j, k);
}

interface BondRow { names: string; types: string; cls: string; ideal: string; k: string; }
function parse(sec: string): BondRow[] {
  const rows: BondRow[] = [];
  for (const line of sec.split('\n')) {
    const m = line.match(/^ ([A-Za-z0-9]+) #\d+\s+([A-Za-z0-9]+) #\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+([-0-9.]+)\s+([0-9.]+)\s+([0-9.]+)/);
    if (m) rows.push({ names: `${m[1]}-${m[2]}`, types: `${m[3]}-${m[4]}`, cls: m[5], ideal: m[7], k: m[9] });
  }
  return rows;
}

const oldT = fs.readFileSync(OLD, 'utf8');
const newT = fs.readFileSync(NEW, 'utf8');
for (const code of CODES) {
  const oldRows = parse(bondSection(oldT, code));
  const newRows = parse(bondSection(newT, code));
  console.log(`\n=== ${code} (${oldRows.length} old / ${newRows.length} new bonds) ===`);
  if (oldRows.length !== newRows.length) { console.log('  LENGTH MISMATCH'); continue; }
  let changed = 0;
  for (let i = 0; i < newRows.length; i++) {
    const o = oldRows[i], n = newRows[i];
    if (o.cls !== n.cls || o.ideal !== n.ideal || o.k !== n.k) {
      changed++;
      console.log(`  ${n.names.padEnd(12)} (${n.types.padEnd(6)})  old cls ${o.cls} r0 ${o.ideal} k ${o.k}  |  new cls ${n.cls} r0 ${n.ideal} k ${n.k}`);
    }
  }
  if (!changed) console.log('  no bond parameter changes');
}
