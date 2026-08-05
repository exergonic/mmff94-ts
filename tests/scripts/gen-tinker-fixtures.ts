// Generate Tinker input files (xyz + prm + key) for the 16 fixture
// molecules so the local Tinker build can cross-check them.
//
// Tinker's MMFF94 uses the ORIGINAL 1-178 type numbering in its xyz
// files, while this library assigns the class/OpenBabel numbering.
// The mapping comes from Tinker's own mmff94.prm atom section
// (`atom <original> <class> <name>`): the first original per class is
// canonical, and energy is class-keyed anyway (analyze maps type →
// class at read time), so any original with the right class is
// equivalent.
//
// Run:  npx tsx tests/scripts/gen-tinker-fixtures.ts
// Output: <TINKER_ROOT>/fixtures16/NAME.xyz + NAME.prm + tinker.key
import { readFileSync, writeFileSync, mkdirSync, copyFileSync, existsSync, readdirSync } from 'fs';
import { join } from 'path';
import { parse_sdf } from '../../src/sdf';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';

const SDF_DIR = 'tests/fixtures/sdf';
const REF_DIR = 'tests/references';
const TINKER_ROOT = 'C:/Users/mccan/Code/tinker';
const OUT = join(TINKER_ROOT, 'fixtures16');

// class (second column) → first original type (first column)
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

// The 16 molecules with obenergy reference logs (the *_non-optimized
// fixtures are for the geometry optimizers and have no reference).
function fixture_names(): string[] {
  return readdirSync(SDF_DIR)
    .filter((f) => f.endsWith('.sdf') && !f.includes('_non-optimized'))
    .map((f) => f.replace('.sdf', ''))
    .filter((name) => existsSync(join(REF_DIR, `${name}.mmff94.log`)))
    .sort();
}

const prmText = readFileSync(join(TINKER_ROOT, 'params', 'mmff94.prm'), 'utf-8');
const map = class_to_original(prmText);
mkdirSync(OUT, { recursive: true });

const names = fixture_names();
console.log(`fixtures: ${names.length} — ${names.join(', ')}`);

for (const name of names) {
  const mol = parse_sdf(readFileSync(join(SDF_DIR, `${name}.sdf`), 'utf-8'));
  const typed = assign_atom_types(mol);

  // Neighbor list per atom (1-based, from the SDF bonds).
  const neighbors: number[][] = mol.atoms.map(() => []);
  for (const b of mol.bonds) {
    neighbors[b.atom1].push(b.atom2 + 1);
    neighbors[b.atom2].push(b.atom1 + 1);
  }

  const lines: string[] = [`${mol.atoms.length} ${name}`];
  for (let i = 0; i < mol.atoms.length; i++) {
    const a = mol.atoms[i];
    const cls = typed.atom_types[i];
    const orig = map.get(cls);
    if (orig === undefined) {
      throw new Error(`${name}: no original type for class ${cls}`);
    }
    const nb = neighbors[i].sort((x, y) => x - y).join(' ');
    lines.push(
      `${String(i + 1).padStart(5)}  ${a.element.padEnd(2)} ${a.x.toFixed(6).padStart(11)} ${a.y.toFixed(6).padStart(11)} ${a.z.toFixed(6).padStart(11)} ${String(orig).padStart(5)}  ${nb}`,
    );
  }
  writeFileSync(join(OUT, `${name}.xyz`), lines.join('\n') + '\n');
  copyFileSync(join(TINKER_ROOT, 'params', 'mmff94.prm'), join(OUT, `${name}.prm`));
}

// MMFF-PIBOND: needed for the conjugated fixtures (benzene, pyridine,
// pyrrole, ethene, formamide, nicotine); harmless elsewhere.
writeFileSync(join(OUT, 'tinker.key'), 'MMFF94\nMMFF-PIBOND\n');
console.log(`wrote ${names.length} xyz + prm files to ${OUT}`);
