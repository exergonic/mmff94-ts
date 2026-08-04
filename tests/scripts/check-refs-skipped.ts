// Tight reference check: does the bmin log have a Total Energy for the skipped?
// Also: what are the '?' elements (dump atom lines).
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';

const ref = JSON.parse(readFileSync('tests/fixtures/validation-suite/mmff94-atom-types.json', 'utf-8'));
const skipped = Object.keys(ref.skipped);
const log = readFileSync('tests/fixtures/validation-suite/MMFF94_bmin.log', 'utf-8');
const text = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');
const mols = parse_mmd(text);
const by_code = new Map(mols.map(m => [m.name, m]));

// Split the log into per-molecule chunks: header line through the next header
const lines = log.split('\n');
const chunks = new Map<string, string>();
let curCode: string | null = null;
const buf: string[] = [];
for (const line of lines) {
  const m = line.match(/\[\s*(\w+),/);
  if (m) {
    if (curCode) chunks.set(curCode, buf.join('\n'));
    curCode = m[1];
    buf.length = 0;
  }
  buf.push(line);
}
if (curCode) chunks.set(curCode, buf.join('\n'));

const has_total = (code: string) => /Total Energy \(kcal\/mol\)=/.test(chunks.get(code) ?? '');
const has_error = (code: string) => {
  const c = chunks.get(code) ?? '';
  return /error|failed|not (found|parametri)|unable/i.test(c) && !has_total(code);
};

let withRef = 0;
let failed = 0;
const noTotal: string[] = [];
for (const code of skipped) {
  if (has_total(code)) withRef++;
  else {
    failed++;
    noTotal.push(code);
  }
}
console.log(`skipped with BatchMin Total Energy: ${withRef} / ${skipped.length}`);
console.log(`skipped WITHOUT reference: ${failed}`);
console.log('no-reference sample:', noTotal.slice(0, 12).join(' '));

// What does the bmin log say for an organic skipped molecule?
for (const code of ['CORWUB10', 'ARGIND11', 'BEWCUB']) {
  const c = chunks.get(code) ?? '';
  const msg = c.split('\n').slice(0, 12).join(' | ');
  console.log(`\n=== ${code} log head: ${msg.slice(0, 400)}`);
}

// The '?' elements: dump the raw atom lines
const qMols = mols.filter(m => skipped.includes(m.name) && m.atoms.some(a => a.element === '?'));
console.log(`\n=== molecules with '?' elements: ${qMols.length}`);
for (const m of qMols.slice(0, 3)) {
  const rawLines = text.split('\n').filter(l => l.includes(`[${m.name},`));
  console.log(m.name, '— atoms:', m.atoms.map(a => a.element).join(','));
}
