// Dump any molecule: element, connectivity, our type vs reference type.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';

const name = process.argv[2] ?? 'COXZEU';
const suiteDir = 'tests/fixtures/validation-suite';
const reference = JSON.parse(readFileSync(`${suiteDir}/mmff94-atom-types.json`, 'utf-8')) as {
  molecules: Record<string, number[]>;
};
const text = readFileSync(`${suiteDir}/MMFF94.mmd`, 'utf-8');
const mols = parse_mmd(text);
const mol = mols.find(m => m.name === name);
if (!mol) { console.log(name + ' not found'); process.exit(1); }

const typed = assign_atom_types(mol);
const refTypes = reference.molecules[name] ?? [];

// Reference pchg from the mmd atom lines
const lines = text.split('\n');
let inM = false;
const refPchg: number[] = [];
for (const line of lines) {
  if (new RegExp(`\\[${name},`).test(line)) { inM = true; continue; }
  if (inM) {
    if (/^\s*\d+\s+\[/.test(line)) break;
    const parts = line.trim().split(/\s+/);
    if (parts.length < 20) continue;
    refPchg.push(parseFloat(parts[19]));
  }
}

console.log(`${name}: ${mol.atoms.length} atoms`);
for (let i = 0; i < mol.atoms.length; i++) {
  const a = mol.atoms[i];
  const nbrs = mol.bonds
    .filter(b => b.atom1 === i || b.atom2 === i)
    .map(b => (b.atom1 === i ? `${b.atom2}:${b.bond_order}` : `${b.atom1}:${b.bond_order}`))
    .join(',');
  const flag = typed.atom_types[i] === refTypes[i] ? '' : '  <<<';
  console.log(
    `${String(i).padStart(3)} ${a.element.padEnd(2)} ${nbrs.padEnd(18)} ` +
    `ours=${String(typed.atom_types[i]).padStart(3)} ref=${String(refTypes[i] ?? '?').padStart(3)}` +
    ` pchg=${(refPchg[i] ?? 0).toFixed(4)}${flag}`,
  );
}
