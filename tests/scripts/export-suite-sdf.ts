/**
 * Export the validation suite as individual SDF files — one per
 * molecule, `<CODE>.sdf` — for viewing in Avogadro or any SDF reader.
 *
 * Source: the suite's own MMFF94.mmd (the BatchMin hypervalent
 * representation — the same input the library's pipeline and the
 * BatchMin references use), MMFF94.titles (names), and
 * MMFF94.fc_hypervalent (formal charges, keyed by atom label).
 *
 * Each SDF carries:
 *   - the MOL V2000 block (atoms, bonds with orders, formal charges
 *     as M CHG lines);
 *   - data fields with the per-atom MMFF94 atom types and the
 *     reference partial charges (the .mmd pchg column), in atom
 *     order.
 *
 * Run: npx tsx tests/scripts/export-suite-sdf.ts
 * Output: tests/fixtures/validation-suite/sdf/<CODE>.sdf
 */
import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { parse_sdf } from '../../src/sdf';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const SUITE = path.join(HERE, '..', 'fixtures', 'validation-suite');
const OUT = path.join(SUITE, 'sdf');

/** V2000 charge codes: 3=+1, 2=+2, 1=+3, 5=−1, 6=−2, 7=−3. */
function charge_code(fc: number): number {
  return fc === 1 ? 3 : fc === 2 ? 2 : fc === 3 ? 1 : fc === -1 ? 5 : fc === -2 ? 6 : fc === -3 ? 7 : 0;
}

function atom_line(x: number, y: number, z: number, element: string, fc: number): string {
  const coords =
    x.toFixed(4).padStart(10) + y.toFixed(4).padStart(10) + z.toFixed(4).padStart(10);
  // V2000: element at 30-32, mass 33-35, charge 36 (1 char), then
  // the parity/H/valence/H0/mapping fields (37-68).
  const tail = '  0'.repeat(12);
  const code = String(charge_code(fc) || ' ');
  return coords + ' ' + element.padEnd(2) + tail.slice(0, 3) + code + tail.slice(4);
}

function bond_line(a1: number, a2: number, order: number): string {
  return String(a1 + 1).padStart(3) + String(a2 + 1).padStart(3) + String(order).padStart(3) + '  0  0  0  0  0  0';
}

function mol_to_sdf(code: string, name: string, m: ReturnType<typeof parse_mmd>[number], typed: ReturnType<typeof assign_atom_types>, fc: Map<number, number>): string {
  const lines: string[] = [];
  lines.push(name ? `${code}: ${name}` : code);
  lines.push('mmff94-ts export-suite-sdf — the MMFF94 validation suite (hypervalent representation)');
  lines.push('');
  const nA = m.atoms.length;
  const nB = m.bonds.length;
  lines.push(`${String(nA).padStart(3)}${String(nB).padStart(3)}  0  0  0  0  0  0  0  0999 V2000`);
  for (let i = 0; i < nA; i++) {
    const a = m.atoms[i];
    lines.push(atom_line(a.x, a.y, a.z, a.element, fc.get(i) ?? 0));
  }
  for (const b of m.bonds) {
    lines.push(bond_line(b.atom1, b.atom2, b.bond_order));
  }
  const charged = [...fc.entries()].filter(([, c]) => c !== 0);
  if (charged.length > 0) {
    // M CHG: 3-wide count, then 4-wide (atomid, charge) pairs.
    let chg = `M  CHG${String(charged.length).padStart(3)}`;
    for (const [i, c] of charged) chg += `${String(i + 1).padStart(4)}${String(c).padStart(4)}`;
    lines.push(chg);
  }
  lines.push('M  END');
  lines.push('> <MMFF94_atom_types>');
  lines.push(typed.atom_types.join(','));
  lines.push('');
  lines.push('> <MMFF94_reference_charges>');
  lines.push(m.atoms.map(a => (a.partial_charge ?? 0).toFixed(5)).join(','));
  lines.push('');
  lines.push('$$$$');
  return lines.join('\n') + '\n';
}

// --- inputs ---------------------------------------------------------
const mols = parse_mmd(fs.readFileSync(path.join(SUITE, 'MMFF94.mmd'), 'utf8'));

const titles = new Map<string, string>();
for (const line of fs.readFileSync(path.join(SUITE, 'MMFF94.titles'), 'utf8').split('\n')) {
  const i = line.indexOf(':');
  if (i > 0) titles.set(line.substring(0, i).trim(), line.substring(i + 1).trim());
}

// Formal charges by molecule + atom label (MMFF94.fc_hypervalent).
const fcByLabel = new Map<string, Map<string, number>>();
let current: string | null = null;
for (const line of fs.readFileSync(path.join(SUITE, 'MMFF94.fc_hypervalent'), 'utf8').split('\n')) {
  const mm = line.match(/^Molecule (\S+)/);
  if (mm) { current = mm[1]; fcByLabel.set(current, new Map()); continue; }
  const am = line.match(/^Atom (\S+), formal charge = (-?\d+)/);
  if (am && current) fcByLabel.get(current)!.set(am[1], parseInt(am[2], 10));
}

// --- export ---------------------------------------------------------
fs.mkdirSync(OUT, { recursive: true });
let missingTitles = 0;
let fcUnmatched = 0;
let fcAtoms = 0;
for (const m of mols) {
  const code = m.name;
  const name = titles.get(code);
  if (name === undefined) missingTitles++;
  const typed = assign_atom_types(m);
  const fc = new Map<number, number>();
  const byLabel = fcByLabel.get(code);
  if (byLabel) {
    const labels = new Set(m.atoms.map(a => a.label).filter(Boolean));
    for (const [lbl, c] of byLabel) {
      // find the atom with this label
      const idx = m.atoms.findIndex(a => a.label === lbl);
      if (idx >= 0) { fc.set(idx, c); fcAtoms++; }
      else fcUnmatched++;
    }
  }
  fs.writeFileSync(path.join(OUT, `${code}.sdf`), mol_to_sdf(code, name ?? '', m, typed, fc));
}

// --- verify: round-trip through our own SDF parser ------------------
let bad = 0;
for (const m of mols) {
  const sdf = fs.readFileSync(path.join(OUT, `${m.name}.sdf`), 'utf8');
  const back = parse_sdf(sdf);
  if (back.atoms.length !== m.atoms.length || back.bonds.length !== m.bonds.length) { bad++; console.log(`  ${m.name}: atom/bond count mismatch`); continue; }
  for (let i = 0; i < m.atoms.length; i++) {
    const a = m.atoms[i], b = back.atoms[i];
    if (a.element !== b.element || Math.abs(a.x - b.x) > 1e-3 || Math.abs(a.y - b.y) > 1e-3 || Math.abs(a.z - b.z) > 1e-3) {
      bad++; console.log(`  ${m.name}: atom ${i} mismatch`); break;
    }
  }
  for (let i = 0; i < m.bonds.length; i++) {
    const a = m.bonds[i], b = back.bonds[i];
    if (a.atom1 !== b.atom1 || a.atom2 !== b.atom2 || a.bond_order !== b.bond_order) {
      bad++; console.log(`  ${m.name}: bond ${i} mismatch`); break;
    }
  }
}

console.log(`exported ${mols.length} molecules to ${OUT}`);
console.log(`titles missing: ${missingTitles}; fc atoms mapped: ${fcAtoms}; fc lookups unmatched: ${fcUnmatched}`);
console.log(`round-trip mismatches: ${bad}`);
