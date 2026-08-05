// Per-term energy coverage scoreboard: for every typing-exact molecule,
// compare all seven terms against the BatchMin log components.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';
import { calc_energy } from '../../src/mmff94/energy/total';

const suiteDir = 'tests/fixtures/validation-suite';
const reference = JSON.parse(readFileSync(`${suiteDir}/mmff94-atom-types.json`, 'utf-8')) as {
  molecules: Record<string, number[]>;
};

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

const bmin = parse_bmin_log(readFileSync(`${suiteDir}/MMFF94_bmin.log`, 'utf-8'));
const mols = parse_mmd(readFileSync(`${suiteDir}/MMFF94.mmd`, 'utf-8'));

// Reference anomalies, excluded from the coverage counts (each is
// documented in tests/charges-suite.test.ts and VALIDATION.md):
// - AN11A and DOZNIP: the anionic 5-ring N⁻ (type 76) — their
//   reference charges are not reproducible from eq. (15) with a
//   uniform q⁰(76) (Halgren's caveat on "unsymmetrical but strongly
//   delocalized anions");
// - FE2PW3 and CU1PW1: BatchMin's van der Waals for the hydrated
//   metal cations predates the X94 metal parameters (the OB — which
//   matches this transcription exactly, 55.84481 vs 55.8448 — and
//   Tinker both use the X94 rows).
const ANOMALY_EXCLUDED = new Set(['AN11A', 'DOZNIP', 'FE2PW3', 'CU1PW1']);

const terms = ['stretch', 'bend', 'strbnd', 'torsion', 'oop', 'vdw', 'elec'] as const;
const totals = { exact: 0, n: 0, worst: 0, worstMol: '' };
const perTerm = Object.fromEntries(terms.map(t => [t, { exact: 0, n: 0, worst: 0, worstMol: '' }])) as Record<
  string, { exact: number; n: number; worst: number; worstMol: string }
>;

for (const mol of mols) {
  if (ANOMALY_EXCLUDED.has(mol.name)) continue;
  const refTypes = reference.molecules[mol.name];
  if (!refTypes || refTypes.length !== mol.atoms.length) continue;
  const typed = assign_atom_types(mol);
  const exact = typed.atom_types.every((t, i) => t === refTypes[i]);
  if (!exact) continue;
  const ref = bmin.get(mol.name);
  if (!ref) continue;
  const e = calc_energy(typed);
  const components = { stretch: e.bond_stretch, bend: e.angle_bend, strbnd: e.stretch_bend, torsion: e.torsion, oop: e.out_of_plane, vdw: e.van_der_waals, elec: e.electrostatic };
  totals.n++;
  let molExact = true;
  for (const t of terms) {
    perTerm[t].n++;
    const d = Math.abs(components[t] - ref[t]);
    if (d > perTerm[t].worst) { perTerm[t].worst = d; perTerm[t].worstMol = mol.name; }
    if (d < 0.05) perTerm[t].exact++;
    else molExact = false;
  }
  const molWorst = Math.max(...terms.map(t => Math.abs(components[t] - ref[t])));
  if (molWorst > totals.worst) { totals.worst = molWorst; totals.worstMol = mol.name; }
  if (molExact) totals.exact++;
  else {
    const bad = terms.filter(t => Math.abs(components[t] - ref[t]) >= 0.05)
      .map(t => `${t}=${Math.abs(components[t] - ref[t]).toFixed(3)}`)
      .join(' ');
    console.log(`  MISMATCH ${mol.name}: ${bad}`);
  }
}

console.log(`typing-exact molecules with bmin refs: ${totals.n}`);
console.log(`all-7-terms within 0.05: ${totals.exact} (${(100 * totals.exact / totals.n).toFixed(1)}%)  worst |d|: ${totals.worst.toFixed(4)} (${totals.worstMol})`);
for (const t of terms) {
  console.log(`  ${t.padEnd(8)} ${String(perTerm[t].exact).padStart(4)}/${perTerm[t].n} within 0.05  worst: ${perTerm[t].worst.toFixed(4)} (${perTerm[t].worstMol})`);
}
