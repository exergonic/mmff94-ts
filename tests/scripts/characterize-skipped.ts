// Characterize the skipped suite molecules with the real parser.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';

const ref = JSON.parse(readFileSync('tests/fixtures/validation-suite/mmff94-atom-types.json', 'utf-8'));
const skipped = Object.keys(ref.skipped);
const text = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');
const mols = parse_mmd(text);
const by_code = new Map(mols.map(m => [m.name, m]));

const counts = new Map<string, number>();
for (const code of skipped) {
  const m = by_code.get(code);
  if (!m) {
    console.log(`${code}: NO BLOCK`);
    continue;
  }
  const els = [...new Set(m.atoms.map(a => a.element))].sort();
  const key = els.join(',');
  counts.set(key, (counts.get(key) ?? 0) + 1);
}
console.log('=== element-set distribution of the 203 skipped ===');
for (const [els, n] of [...counts.entries()].sort((a, b) => b[1] - a[1])) {
  console.log(`  ${n.toString().padStart(3)}  ${els}`);
}
