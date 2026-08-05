// Charge sweep over the full 753: worst per-atom deviation per molecule.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { assign_bci_charges } from '../../src/mmff94/charges';

const suiteDir = 'tests/fixtures/validation-suite';
const refs = JSON.parse(readFileSync(`${suiteDir}/mmff94-atom-types.json`, 'utf-8')) as {
  molecules: Record<string, number[]>;
};
const mmdText = readFileSync(`${suiteDir}/MMFF94.mmd`, 'utf-8');

function reference_charges(name: string, nAtoms: number): number[] {
  const pchg = new Array(nAtoms).fill(0);
  let inMol = false;
  for (const line of mmdText.split('\n')) {
    const head = line.match(/^\s*\d+\s+\[(\w+),/);
    if (head) { inMol = head[1] === name; continue; }
    if (!inMol) continue;
    const p = line.trim().split(/\s+/);
    if (p.length >= 20) {
      const serial = parseInt(p[p.length - 1], 10);
      if (!isNaN(serial)) pchg[serial - 1] = parseFloat(p[p.length - 4]);
    }
  }
  return pchg;
}

const bad: { code: string; worst: number; atoms: string }[] = [];
let checked = 0;
let worstAll = 0;
for (const mol of parse_mmd(mmdText)) {
  const refTypes = refs.molecules[mol.name];
  if (!refTypes || refTypes.length !== mol.atoms.length) continue;
  const typed = assign_atom_types(mol);
  if (!typed.atom_types.every((t, i) => t === refTypes[i])) continue;
  if (mol.name === 'JALSOE' || mol.name === 'SO18A') continue;
  const charged = assign_bci_charges(typed);
  const ref = reference_charges(mol.name, mol.atoms.length);
  checked++;
  let worst = 0;
  const worstAtoms: string[] = [];
  for (let i = 0; i < typed.atoms.length; i++) {
    const dev = Math.abs(charged.partial_charges![i] - ref[i]);
    if (dev > worst) worst = dev;
    if (dev > 1e-3) worstAtoms.push(`${i}:${typed.atoms[i].element}[${typed.atom_types[i]}] got ${charged.partial_charges![i].toFixed(3)} ref ${ref[i].toFixed(3)}`);
    if (dev > worstAll) worstAll = dev;
  }
  if (worst > 1e-3) bad.push({ code: mol.name, worst, atoms: worstAtoms.join('; ') });
}
bad.sort((a, b) => b.worst - a.worst);
console.log(`checked ${checked}, exact (all atoms <1e-3): ${checked - bad.length}, worst overall: ${worstAll}`);
for (const b of bad.slice(0, 25)) {
  console.log(`  ${b.code}: worst ${b.worst.toFixed(4)}`);
  if (process.argv[2] === '--atoms') console.log(`      ${b.atoms}`);
}
console.log(`total mismatched molecules: ${bad.length}`);
