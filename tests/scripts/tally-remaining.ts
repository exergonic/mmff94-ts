// Tally the ours-vs-ref type mismatches across the whole suite.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';

const ref = JSON.parse(readFileSync('tests/fixtures/validation-suite/mmff94-atom-types.json', 'utf-8'));
const text = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');
const mols = parse_mmd(text);
const by_code = new Map(mols.map(m => [m.name, m]));

const pairs = new Map<string, number>();
const molCount = new Map<string, number>();
for (const [code, types] of Object.entries(ref.molecules)) {
  const m = by_code.get(code);
  if (!m) continue;
  const ours = assign_atom_types(m);
  let worst = 0;
  const fam = new Set<string>();
  for (let i = 0; i < m.atoms.length; i++) {
    if (ours.atom_types[i] !== (types as number[])[i]) {
      const key = `${(types as number[])[i]}<-${ours.atom_types[i]}`;
      pairs.set(key, (pairs.get(key) ?? 0) + 1);
      fam.add(key);
      worst++;
    }
  }
  if (worst > 0) {
    molCount.set([...fam].join(' '), (molCount.get([...fam].join(' ')) ?? 0) + 1);
  }
}
console.log('=== type-pair mismatches (ref<-ours): count of atoms ===');
for (const [k, n] of [...pairs.entries()].sort((a, b) => b[1] - a[1]).slice(0, 25)) {
  console.log(`  ${k.padEnd(12)} ${n}`);
}
console.log('\n=== molecule family signatures (count of molecules) ===');
for (const [k, n] of [...molCount.entries()].sort((a, b) => b[1] - a[1]).slice(0, 20)) {
  console.log(`  ${n.toString().padStart(3)}  ${k}`);
}
