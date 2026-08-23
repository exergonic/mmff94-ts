// Generate a Tinker xyz for trpcage (single molecule) so the local
// Tinker build can arbitrate the zwitterion electrostatics question.
// Same conventions as gen-tinker-fixtures.ts: class -> first original
// type mapping from mmff94.prm, \n endings, matching <name>.prm, and
// the same tinker.key (MMFF94 + MMFF-PIBOND; PIBOND is harmless here —
// no conjugated bonds — but keeping one key file for all fixtures
// avoids a second configuration to maintain).
//
// The SDF's N-terminal N carries formal charge +1 (V2000 charge code
// 4), which becomes the xyz line's optional formal-charge integer —
// without it TINKER reads the zwitterion as neutral and its charges
// (and electrostatics) diverge from ours.
//
// Run:  npx tsx tests/scripts/gen-trpcage-xyz.ts
import { readFileSync, writeFileSync, mkdirSync, copyFileSync } from 'fs';
import { join } from 'path';
import { parse_sdf } from '../../src/sdf';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';

const SDF = 'tests/fixtures/sdf/trpcage.sdf';
const TINKER_ROOT = 'C:/Users/mccan/Code/tinker';
const OUT = join(TINKER_ROOT, 'fixtures16');

function class_to_original(prmText: string): Map<number, number> {
  const map = new Map<number, number>();
  for (const line of prmText.split('\n')) {
    const m = line.match(/^atom\s+(\d+)\s+(\d+)\s/);
    if (!m) continue;
    const orig = parseInt(m[1], 10);
    const cls = parseInt(m[2], 10);
    if (!map.has(cls)) map.set(cls, orig);
  }
  return map;
}

const mol = parse_sdf(readFileSync(SDF, 'utf-8'));
mol.name = 'trpcage';
const typed = assign_atom_types(mol);

// Formal charge rides on the xyz line (the integer after the type).
const neighbors: number[][] = mol.atoms.map(() => []);
for (const b of mol.bonds) {
  neighbors[b.atom1].push(b.atom2 + 1);
  neighbors[b.atom2].push(b.atom1 + 1);
}

const prmText = readFileSync(join(TINKER_ROOT, 'params', 'mmff94.prm'), 'utf-8');
const map = class_to_original(prmText);

mkdirSync(OUT, { recursive: true });
const lines: string[] = [`${mol.atoms.length} trpcage`];
for (let i = 0; i < mol.atoms.length; i++) {
  const a = mol.atoms[i];
  const cls = typed.atom_types[i];
  const orig = map.get(cls);
  if (orig === undefined) throw new Error(`no original type for class ${cls} (atom ${i})`);
  const fc = a.formal_charge === 0 ? '' : `  ${a.formal_charge > 0 ? 1 : -1}`;
  const nb = neighbors[i].sort((x, y) => x - y).join(' ');
  lines.push(
    `${String(i + 1).padStart(5)}  ${a.element.padEnd(2)} ${a.x.toFixed(6).padStart(11)} ${a.y.toFixed(6).padStart(11)} ${a.z.toFixed(6).padStart(11)} ${String(orig).padStart(5)}${fc}  ${nb}`,
  );
}
writeFileSync(join(OUT, 'trpcage.xyz'), lines.join('\n') + '\n');
copyFileSync(join(TINKER_ROOT, 'params', 'mmff94.prm'), join(OUT, 'trpcage.prm'));
console.log(`wrote ${OUT}/trpcage.xyz (${mol.atoms.length} atoms, N-term fc=${mol.atoms[0].formal_charge})`);
