/**
 * Probe: the reference's bond-class table, single bonds only (bond order
 * from the mmd), with per-occurrence context for mixed-class pairs.
 *
 * Run: npx tsx tests/scripts/probe-btij-table.ts
 */
import * as fs from 'fs';

const OPTI = 'C:/Users/mccan/AppData/Local/Temp/mmff94-761/MMFF94_opti.log';
const MMD = 'C:/Users/mccan/Code/mmff94-ts/tests/fixtures/validation-suite/MMFF94.mmd';
const opti = fs.readFileSync(OPTI, 'utf8');
const mmd = fs.readFileSync(MMD, 'utf8');

// parse mmd: per molecule, name -> atom index, bonds (i,j,order)
interface MolInfo { names: Map<string, number>; bonds: Map<string, number>; }
function parseMmd(): Map<string, MolInfo> {
  const result = new Map<string, MolInfo>();
  const blocks = mmd.split(/(?=\[\w+,)/);
  for (const b of blocks) {
    const hdr = b.match(/^\[(\w+),/);
    if (!hdr) continue;
    const info: MolInfo = { names: new Map(), bonds: new Map() };
    const lines = b.split('\n');
    for (const line of lines) {
      const p = line.trim().split(/\s+/);
      if (p.length < 4 || !/^\d+$/.test(p[0])) continue;
      // format: <el-code> <nb1> <ord1> <nb2> <ord2> ... <x> <y> <z> <label> <fchg> <pchg> <pchg2> <name> <label> <num>
      // find the element label: token at -3 is the atom label (e.g. "S1"), -2 the molecule name
      const label = p[p.length - 3];
      const idx = info.names.size;
      info.names.set(label, idx);
      const bonds = p.slice(1);
      for (let bi = 0; bi + 1 < bonds.length; bi += 2) {
        const nb = Number(bonds[bi]);
        if (!nb) continue;
        const order = Number(bonds[bi + 1]);
        const other = info.names.get('') ?? undefined; // not available yet — resolve after
        // store pending: (neighborNum, order) with 1-based neighbor numbers
        // resolve after all names known — so store by number
        info.bonds.set(`${Math.min(0, nb)}`, 0); // placeholder
      }
    }
    result.set(hdr[1], info);
  }
  return result;
}

// simpler: parse atom lines -> (num, label); bonds as (num, nb, order); resolve later
function parseMmd2(): Map<string, { labels: string[]; bonds: [number, number, number][] }> {
  const result = new Map();
  const blocks = mmd.split(/(?=\[\w+,)/);
  for (const b of blocks) {
    const hdr = b.match(/^\[(\w+),/);
    if (!hdr) continue;
    const labels: string[] = [];
    const bonds: [number, number, number][] = [];
    for (const line of b.split('\n')) {
      const p = line.trim().split(/\s+/);
      if (p.length < 4 || !/^\d+$/.test(p[0])) continue;
      const num = labels.length + 1;
      const label = p[p.length - 2];
      labels.push(label);
      const rest = p.slice(1);
      for (let bi = 0; bi + 1 < rest.length; bi += 2) {
        const nb = Number(rest[bi]);
        if (!nb) continue;
        bonds.push([num, nb, Number(rest[bi + 1])]);
      }
    }
    result.set(hdr[1], { labels, bonds });
  }
  return result;
}

const mols = parseMmd2();
const molStarts = [...opti.matchAll(/New Structure Name\/Conformational Index: (\w+)|^ Structure Name: (\w+)/gm)];

interface Row { mol: string; a: string; b: string; ti: number; tj: number; cls: number; order: number; }
const rows: Row[] = [];
for (let m = 0; m < molStarts.length; m++) {
  const code = molStarts[m][1] ?? molStarts[m][2];
  const secStart = molStarts[m].index;
  const secEnd = m + 1 < molStarts.length ? molStarts[m + 1].index : opti.length;
  const sec = opti.slice(secStart, secEnd);
  const j = sec.indexOf('OPTIMOL-ANALYZE>  # bonds');
  if (j < 0) continue;
  const k = sec.indexOf('OPTIMOL-ANALYZE>  # ', j + 10);
  const block = sec.slice(j, k < 0 ? sec.length : k);
  const info = mols.get(code);
  const orderOf = (na: string, nb: string): number => {
    if (!info) return 1;
    const ia = info.labels.indexOf(na) + 1;
    const ib = info.labels.indexOf(nb) + 1;
    if (ia <= 0 || ib <= 0) return 1;
    const b = info.bonds.find(([x, y]) => (x === ia && y === ib) || (x === ib && y === ia));
    return b ? b[2] : 1;
  };
  for (const line of block.split('\n')) {
    const mm = line.match(/^ (\S+) #\d+\s+(\S+) #\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/);
    if (!mm) continue;
    rows.push({
      mol: code, a: mm[1], b: mm[2],
      ti: Number(mm[3]), tj: Number(mm[4]), cls: Number(mm[5]),
      order: orderOf(mm[1], mm[2]),
    });
  }
}
console.log('rows:', rows.length);
const single = rows.filter(r => r.order === 1);
console.log('single bonds:', single.length);
const pairClasses = new Map<string, Set<number>>();
for (const r of single) {
  const key = r.ti <= r.tj ? `${r.ti}-${r.tj}` : `${r.tj}-${r.ti}`;
  if (!pairClasses.has(key)) pairClasses.set(key, new Set());
  pairClasses.get(key)!.add(r.cls);
}
const mixed = [...pairClasses.entries()].filter(([, s]) => s.has(0) && s.has(1)).map(([k]) => k).sort();
console.log('mixed pairs:', mixed.length, mixed.join(' '));
// per-occurrence context for the mixed pairs
for (const key of mixed) {
  const [a, b] = key.split('-').map(Number);
  const occ = single.filter(r => (r.ti === a && r.tj === b) || (r.ti === b && r.tj === a));
  const cls1 = occ.filter(r => r.cls === 1);
  console.log(`\n${key}: ${occ.length} occurrences, class 1 in ${cls1.length}:`);
  for (const r of cls1) console.log(`   ${r.mol} ${r.a}-${r.b} cls 1`);
}
