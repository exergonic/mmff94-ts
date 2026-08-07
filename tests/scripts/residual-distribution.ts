/**
 * The residual DISTRIBUTION vs BatchMin — the honest accuracy claim
 * for the compliance statement. Extends the scoreboard's worst-only
 * tracking with per-term and total max/mean/RMS and counts at the
 * 1e-5 / 5e-5 (Wavefun's claim) / 1e-4 (our census) gates.
 *
 * Run: npx tsx tests/scripts/residual-distribution.ts
 */
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { calc_energy } from '../../src/mmff94/energy/total';
import { parse_bmin_log } from './bmin-log';

const suiteDir = 'tests/fixtures/validation-suite';
const reference = JSON.parse(readFileSync(`${suiteDir}/mmff94-atom-types.json`, 'utf-8')) as {
  molecules: Record<string, number[]>;
};
const bmin = parse_bmin_log(readFileSync(`${suiteDir}/MMFF94_bmin.log`, 'utf-8'));
const mols = parse_mmd(readFileSync(`${suiteDir}/MMFF94.mmd`, 'utf-8'));

const PER_TERM_EXCLUDED: Record<string, string[]> = {
  AN11A: ['elec'], DOZNIP: ['elec'], FE2PW3: ['vdw'], CU1PW1: ['vdw'],
};

const terms = ['stretch', 'bend', 'strbnd', 'torsion', 'oop', 'vdw', 'elec'] as const;

type Stat = { n: number; max: number; maxMol: string; sum: number; sumSq: number; le5: number; le5e5: number; le4: number };
const newStat = (): Stat => ({ n: 0, max: 0, maxMol: '', sum: 0, sumSq: 0, le5: 0, le5e5: 0, le4: 0 });
const perTerm = Object.fromEntries(terms.map(t => [t, newStat()])) as Record<string, Stat>;
const totals = newStat();
// how many molecules are fully within each gate on ALL comparable terms
const full = { le5: 0, le5e5: 0, le4: 0, n: 0 };

for (const mol of mols) {
  const excluded = PER_TERM_EXCLUDED[mol.name] ?? [];
  const refTypes = reference.molecules[mol.name];
  if (!refTypes || refTypes.length !== mol.atoms.length) continue;
  const typed = assign_atom_types(mol);
  const exact = typed.atom_types.every((t, i) => t === refTypes[i]);
  if (!exact) continue;
  const ref = bmin.get(mol.name);
  if (!ref) continue;
  const e = calc_energy(typed);
  const components = { stretch: e.bond_stretch, bend: e.angle_bend, strbnd: e.stretch_bend, torsion: e.torsion, oop: e.out_of_plane, vdw: e.van_der_waals, elec: e.electrostatic };
  const comparable = terms.filter(t => !excluded.includes(t));
  full.n++;
  let molLe5 = true, molLe5e5 = true, molLe4 = true;
  for (const t of comparable) {
    const s = perTerm[t];
    s.n++;
    const d = Math.abs(components[t] - ref[t]);
    s.sum += d; s.sumSq += d * d;
    if (d > s.max) { s.max = d; s.maxMol = mol.name; }
    if (d < 1e-5) s.le5++; else molLe5 = false;
    if (d < 5e-5) s.le5e5++; else molLe5e5 = false;
    if (d < 1e-4) s.le4++; else molLe4 = false;
  }
  const dTot = Math.abs(e.total - ref.total);
  // The totals census convention: a molecule with ANY excluded term
  // (the documented anomalies) is dropped from the total stats too.
  if (excluded.length > 0) { full.n++; continue; }
  totals.n++;
  totals.sum += dTot; totals.sumSq += dTot * dTot;
  if (dTot > totals.max) { totals.max = dTot; totals.maxMol = mol.name; }
  if (dTot < 1e-5) totals.le5++;
  if (dTot < 5e-5) totals.le5e5++;
  if (dTot < 1e-4) totals.le4++;
  if (molLe5) full.le5++;
  if (molLe5e5) full.le5e5++;
  if (molLe4) full.le4++;
}

const row = (label: string, s: Stat) => {
  const mean = s.sum / s.n;
  const rms = Math.sqrt(s.sumSq / s.n);
  console.log(
    `${label.padEnd(9)} n=${String(s.n).padStart(4)}  max=${s.max.toExponential(2).padStart(10)} (${s.maxMol})  mean=${mean.toExponential(2).padStart(9)}  rms=${rms.toExponential(2).padStart(9)}  ≤1e-5: ${String(s.le5).padStart(4)}  ≤5e-5: ${String(s.le5e5).padStart(4)}  ≤1e-4: ${String(s.le4).padStart(4)}`,
  );
};
console.log('Per-term residuals vs BatchMin (typing-exact molecules, documented exclusions applied):');
for (const t of terms) row(t, perTerm[t]);
console.log('');
row('TOTAL', totals);
console.log('');
console.log(`Molecules fully within the gate on ALL comparable terms: ≤1e-5: ${full.le5}/${full.n}  ≤5e-5: ${full.le5e5}/${full.n}  ≤1e-4: ${full.le4}/${full.n}`);
