// Residual census: how far are we from exact per-term reproduction?
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';
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
const ANOMALY = new Set(['AN11A', 'DOZNIP', 'FE2PW3', 'CU1PW1', 'JALSOE', 'SO18A']);

const terms = ['stretch', 'bend', 'strbnd', 'torsion', 'oop', 'vdw', 'elec'] as const;
const gk = { stretch: 'bond_stretch', bend: 'angle_bend', strbnd: 'stretch_bend', torsion: 'torsion', oop: 'out_of_plane', vdw: 'van_der_waals', elec: 'electrostatic' } as const;

const rows: { mol: string; term: string; d: number }[] = [];
for (const mol of mols) {
  if (ANOMALY.has(mol.name)) continue;
  const ref = bmin.get(mol.name);
  if (!ref) continue;
  const typed = assign_atom_types(mol);
  const got = calc_energy(typed);
  for (const t of terms) {
    rows.push({ mol: mol.name, term: t, d: Math.abs(got[gk[t]] - ref[t]) });
  }
}

console.log('per-term bins (count of molecules with |d| in each range):');
for (const t of terms) {
  const rs = rows.filter(r => r.term === t);
  const b1e4 = rs.filter(r => r.d <= 1e-4).length;
  const b1e3 = rs.filter(r => r.d > 1e-4 && r.d <= 1e-3).length;
  const b1e2 = rs.filter(r => r.d > 1e-3 && r.d <= 1e-2).length;
  const b5e2 = rs.filter(r => r.d > 1e-2).length;
  const worst = rs.sort((a, b) => b.d - a.d).slice(0, 3);
  console.log(`  ${t.padEnd(8)} <=1e-4: ${b1e4}  <=1e-3: ${b1e3}  <=1e-2: ${b1e2}  >1e-2: ${b5e2}   worst: ${worst.map(w => `${w.mol} ${w.d.toExponential(1)}`).join(', ')}`);
}
const worst10 = rows.sort((a, b) => b.d - a.d).slice(0, 10);
console.log('\ntop-10 residuals overall:');
for (const w of worst10) console.log(`  ${w.mol.padEnd(10)} ${w.term.padEnd(8)} ${w.d.toExponential(2)}`);
