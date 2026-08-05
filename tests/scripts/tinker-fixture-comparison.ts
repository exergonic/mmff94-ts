// Three-way comparison on the 16 fixture molecules: mmff94-ts vs
// OpenBabel (the committed obenergy logs) vs Tinker (the local build's
// analyze output, committed under tests/references/tinker/).
//
// Run:  npx tsx tests/scripts/tinker-fixture-comparison.ts
import { readFileSync, readdirSync, existsSync } from 'fs';
import { join } from 'path';
import { parse_sdf } from '../../src/sdf';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { calc_energy } from '../../src/mmff94/energy/total';

const SDF_DIR = 'tests/fixtures/sdf';
const REF_DIR = 'tests/references';
const TINKER_LOG_DIR = 'tests/references/tinker';

const TERMS = [
  ['bond_stretch', 'TOTAL BOND STRETCHING ENERGY'],
  ['angle_bend', 'TOTAL ANGLE BENDING ENERGY'],
  ['stretch_bend', 'TOTAL STRETCH BENDING ENERGY'],
  ['torsion', 'TOTAL TORSIONAL ENERGY'],
  ['out_of_plane', 'TOTAL OUT-OF-PLANE BENDING ENERGY'],
  ['van_der_waals', 'TOTAL VAN DER WAALS ENERGY'],
  ['electrostatic', 'TOTAL ELECTROSTATIC ENERGY'],
] as const;

function parse_ob_log(text: string): Record<string, number> {
  const out: Record<string, number> = {};
  for (const line of text.split('\n')) {
    const m = line.match(/TOTAL (\w[\w ()-]*?) ENERGY\s*=\s*([-\d.]+)/);
    if (m) out[m[1].trim()] = parseFloat(m[2]);
  }
  const t = text.match(/TOTAL ENERGY\s*=\s*([-\d.]+)/);
  if (t) out.total = parseFloat(t[1]);
  return out;
}

function parse_tinker_log(text: string): Record<string, number> {
  const out: Record<string, number> = {};
  const map: Record<string, string> = {
    'Bond Stretching': 'bond_stretch',
    'Angle Bending': 'angle_bend',
    'Stretch-Bend': 'stretch_bend',
    'Torsional Angle': 'torsion',
    'Out-of-Plane Bend': 'out_of_plane',
    'Van der Waals': 'van_der_waals',
    'Charge-Charge': 'electrostatic',
  };
  for (const line of text.split('\n')) {
    for (const [label, key] of Object.entries(map)) {
      const m = line.match(new RegExp(`^\\s*${label}\\s+([-\\d.]+)`));
      if (m) out[key] = parseFloat(m[1]);
    }
    const t = line.match(/^\s*Total Potential Energy\s*:\s*([-\d.]+)/);
    if (t) out.total = parseFloat(t[1]);
  }
  // Tinker omits zero terms from the listing; a missing term is 0.
  for (const key of Object.values(map)) if (!(key in out)) out[key] = 0;
  return out;
}

const names = readdirSync(SDF_DIR)
  .filter((f) => f.endsWith('.sdf') && !f.includes('_non-optimized'))
  .map((f) => f.replace('.sdf', ''))
  .filter((name) => existsSync(join(REF_DIR, `${name}.mmff94.log`)))
  .sort();

let totalMax = 0;
let totalWorst = '';
const termMax: Record<string, { d: number; mol: string }> = {};
for (const t of TERMS) termMax[t[0]] = { d: 0, mol: '' };
let termsWithin5e5 = 0;
let termsTotal = 0;

for (const name of names) {
  const mol = parse_sdf(readFileSync(join(SDF_DIR, `${name}.sdf`), 'utf-8'));
  const typed = assign_atom_types(mol);
  const got = calc_energy(typed);
  const ob = parse_ob_log(readFileSync(join(REF_DIR, `${name}.mmff94.log`), 'utf-8'));
  const tk = parse_tinker_log(readFileSync(join(TINKER_LOG_DIR, `${name}.log`), 'utf-8'));

  console.log(`\n${name}:`);
  console.log('  term       ours       OpenBabel   Δob       Tinker     Δtk');
  for (const [key, obLabel] of TERMS) {
    const ours = got[key] as number;
    const obKey = obLabel.replace('TOTAL ', '').replace(' ENERGY', '');
    const obv = ob[obKey];
    const tkv = tk[key];
    const dob = obv === undefined ? NaN : ours - obv;
    const dtk = tkv === undefined ? NaN : ours - tkv;
    if (!isNaN(dtk)) {
      termsTotal++;
      if (Math.abs(dtk) <= 5e-5) termsWithin5e5++;
      if (Math.abs(dtk) > termMax[key].d) termMax[key] = { d: Math.abs(dtk), mol: name };
    }
    console.log(
      `  ${key.padEnd(9)} ${ours.toFixed(5).padStart(10)} ${(obv ?? NaN).toFixed(5).padStart(10)} ${dob.toFixed(5).padStart(10)} ${(tkv ?? NaN).toFixed(5).padStart(10)} ${dtk.toFixed(5).padStart(10)}`,
    );
  }
  const dtot = Math.abs(got.total - tk.total);
  if (dtot > totalMax) { totalMax = dtot; totalWorst = name; }
  console.log(
    `  total      ${got.total.toFixed(5).padStart(10)} ${ob.total.toFixed(5).padStart(10)} ${(got.total - ob.total).toFixed(5).padStart(10)} ${tk.total.toFixed(5).padStart(10)} ${(got.total - tk.total).toFixed(5).padStart(10)}`,
  );
}

console.log('\n=== summary ===');
console.log(`totals vs Tinker: worst |Δ| = ${totalMax.toFixed(5)} kcal/mol (${totalWorst})`);
console.log(`per-term vs Tinker within 5e-5: ${termsWithin5e5}/${termsTotal}`);
for (const t of TERMS) {
  console.log(`  ${t[0].padEnd(8)} worst |Δ| vs Tinker: ${termMax[t[0]].d.toExponential(1)} (${termMax[t[0]].mol})`);
}
