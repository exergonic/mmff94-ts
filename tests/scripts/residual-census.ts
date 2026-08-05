// Residual census: how far are we from exact per-term reproduction?
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { calc_energy } from '../../src/mmff94/energy/total';
import { join } from 'path';
import { load_bmin_log } from './bmin-log';

const suiteDir = 'tests/fixtures/validation-suite';
const bmin = load_bmin_log(suiteDir);
const mols = parse_mmd(readFileSync(join(suiteDir, 'MMFF94.mmd'), 'utf-8'));
// Reference anomalies — per-term, as documented in VALIDATION.md.
// JALSOE/SO18A left this set on 2026-08-05 (three-way verified).
const PER_TERM_EXCLUDED: Record<string, string[]> = {
  AN11A: ['elec'],
  DOZNIP: ['elec'],
  FE2PW3: ['vdw'],
  CU1PW1: ['vdw'],
};

const terms = ['stretch', 'bend', 'strbnd', 'torsion', 'oop', 'vdw', 'elec'] as const;
const gk = { stretch: 'bond_stretch', bend: 'angle_bend', strbnd: 'stretch_bend', torsion: 'torsion', oop: 'out_of_plane', vdw: 'van_der_waals', elec: 'electrostatic' } as const;

const rows: { mol: string; term: string; d: number }[] = [];
for (const mol of mols) {
  const excluded = PER_TERM_EXCLUDED[mol.name] ?? [];
  if (excluded.length === terms.length) continue;
  const ref = bmin.get(mol.name);
  if (!ref) continue;
  const typed = assign_atom_types(mol);
  const got = calc_energy(typed);
  for (const t of terms) {
    if (excluded.includes(t)) continue;
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
