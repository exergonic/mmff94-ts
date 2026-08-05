// Is the bend band (168 molecules at 1e-4..1e-3) a conversion constant?
// The paper's 0.043844 vs the exact 2xMM2-derived 0.04382836, and the
// cubic constant variants.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { calc_energy } from '../../src/mmff94/energy/total';
import { join } from 'path';

const suiteDir = 'tests/fixtures/validation-suite';
function parse_fortran(s: string): number {
  s = s.trim();
  if (s.includes('D')) s = s.replace('D', 'e');
  return parseFloat(s);
}
function parse_bmin_log(text: string): Map<string, Record<string, number>> {
  const result = new Map<string, Record<string, number>>();
  let currentCode = '';
  const lines = text.split('\n');
  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(/\[\s*(\w+),/);
    if (m) { currentCode = m[1]; continue; }
    if (!currentCode) continue;
    const s = lines[i].match(/^\s+Stretch\s*=\s*(\S+)/);
    if (s) {
      const e: Record<string, number> = {};
      e.stretch = parse_fortran(s[1]);
      e.bend = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
      e.torsion = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
      e.oop = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
      e.strbnd = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
      e.elec = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
      e.vdw = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
      result.set(currentCode, e);
      currentCode = '';
    }
  }
  return result;
}

const bmin = parse_bmin_log(readFileSync(join(suiteDir, 'MMFF94_bmin.log'), 'utf-8'));
const mols = parse_mmd(readFileSync(join(suiteDir, 'MMFF94.mmd'), 'utf-8'));

// Which molecules are in the bend band?
const band: { mol: string; d: number; bend: number; rel: number }[] = [];
for (const mol of mols) {
  const ref = bmin.get(mol.name);
  if (!ref) continue;
  const typed = assign_atom_types(mol);
  const got = calc_energy(typed);
  const d = Math.abs(got.angle_bend - ref.bend);
  if (d > 1e-4 && d <= 1e-2) {
    band.push({ mol: mol.name, d, bend: ref.bend, rel: d / Math.max(0.1, Math.abs(ref.bend)) });
  }
}
console.log(`bend band molecules: ${band.length}`);
// distribution of the relative deltas
const rels = band.map(b => b.rel).sort((a, b) => a - b);
const med = rels[Math.floor(rels.length / 2)];
console.log(`relative |d|/|E| median ${med.toFixed(5)}, p90 ${rels[Math.floor(rels.length*0.9)].toFixed(5)}, max ${rels[rels.length-1].toFixed(5)}`);
// GAKTAN details
const g = band.find(b => b.mol === 'GAKTAN')!;
console.log('GAKTAN: d', g.d.toFixed(5), 'bend', g.bend.toFixed(4), 'rel', g.rel.toFixed(5));
// sign pattern: do we over or under?
let over = 0, under = 0;
for (const mol of mols) {
  const ref = bmin.get(mol.name);
  if (!ref) continue;
  const got = calc_energy(assign_atom_types(mol));
  const d = got.angle_bend - ref.bend;
  if (Math.abs(d) > 1e-4) d > 0 ? over++ : under++;
}
console.log(`sign: over ${over}, under ${under}`);
