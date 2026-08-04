// List the non-exact molecules with their mismatch pairs, grouped by family.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';

const ref = JSON.parse(readFileSync('tests/fixtures/validation-suite/mmff94-atom-types.json', 'utf-8'));
const text = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');
const mols = parse_mmd(text);
const by_code = new Map(mols.map(m => [m.name, m]));

const want = process.argv[2] ?? ''; // e.g. '62<-8' or '53<-9'
for (const [code, types] of Object.entries(ref.molecules)) {
  const m = by_code.get(code);
  if (!m) continue;
  const ours = assign_atom_types(m);
  const pairs: string[] = [];
  for (let i = 0; i < m.atoms.length; i++) {
    if (ours.atom_types[i] !== (types as number[])[i]) {
      pairs.push(`${m.atoms[i].element}${i}:${(types as number[])[i]}<-${ours.atom_types[i]}`);
    }
  }
  if (pairs.length === 0) continue;
  if (want && !pairs.some(p => p.includes(want))) continue;
  console.log(`${code.padEnd(10)} ${pairs.join('; ')}`);
}
