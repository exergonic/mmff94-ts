/**
 * Hard per-term regression gate.
 *
 * Every typing-exact suite molecule's per-term residual must be ≤ 1e-4
 * kcal/mol. The exceptions are the documented coarse-precision
 * generated-bond rows of the ERULE fragments (ERULE_03/06 stretch,
 * ERULE_03 strbnd), where the reference prints the generated parameter
 * to 3 decimals — those pin at their measured residual with the reason
 * stated. The 863a70c torsion regression (0.29–0.41 kcal/mol on four
 * ERULE molecules) passed `npm test` for two days; this gate catches the
 * next one.
 *
 * The gate runs as part of `npm run test`. Detailed per-molecule
 * diagnostics live in tests/scripts/residual-distribution.ts and
 * energy-scoreboard.ts.
 */

import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { parse_mmd } from '../src/utils/mmd-parser.js';
import { assign_atom_types } from '../src/mmff94/assign-atom-types.js';
import { calc_energy } from '../src/mmff94/energy/total.js';
import { parse_bmin_log } from './scripts/bmin-log.js';

const suiteDir = join(__dirname, 'fixtures', 'validation-suite');
const molecules = parse_mmd(readFileSync(join(suiteDir, 'MMFF94.mmd'), 'utf-8'));
const refEnergies = parse_bmin_log(readFileSync(join(suiteDir, 'MMFF94_bmin.log'), 'utf-8'));
const reference = JSON.parse(
  readFileSync(join(suiteDir, 'mmff94-atom-types.json'), 'utf-8'),
) as { molecules: Record<string, number[]> };

const DEFAULT_TOL = 1e-4;

// Excluded molecules/terms: the reference itself is inconsistent for
// these (the type-76 anionic nitrogen, AN11A/DOZNIP — Halgren's own
// caveat: strongly delocalized anions have no uniform charge
// assignment; BatchMin, Tinker, OpenBabel, and us all disagree). Their
// other six terms are verified by the gate.
const EXCLUDED: Record<string, TermName[]> = {
  AN11A: ['electrostatic'],
  DOZNIP: ['electrostatic'],
};

// Coarse-precision generated-bond rows: the reference prints the generated
// parameter to 3 decimals, so the residual is bounded by the reference's
// own print precision, not the 1e-4 gate. The tolerance here pins the
// measured residual; a future drift fails.
const EXCEPTIONS: Record<string, Record<string, { tol: number; reason: string }>> = {
  ERULE_03: {
    bond_stretch: { tol: 2.0e-3, reason: 'generated P–Si bond at reference 3-dp print precision' },
    stretch_bend: { tol: 4.0e-4, reason: 'inherited from the P–Si generated bond' },
  },
  ERULE_06: {
    bond_stretch: { tol: 2.0e-3, reason: 'generated F–N bond at reference 3-dp print precision' },
  },
};

const termKey = {
  bond_stretch: 'stretch',
  angle_bend: 'bend',
  stretch_bend: 'strbnd',
  torsion: 'torsion',
  out_of_plane: 'oop',
  van_der_waals: 'vdw',
  electrostatic: 'elec',
} as const;

type TermName = keyof typeof termKey;

const TERM_ORDER: TermName[] = [
  'bond_stretch', 'angle_bend', 'stretch_bend',
  'torsion', 'out_of_plane', 'van_der_waals', 'electrostatic',
];

describe('per-term regression gate (suite, typing-exact molecules)', () => {
  const typed = molecules.map(m => ({
    mol: m,
    typed: assign_atom_types(m),
    refTypes: reference.molecules[m.name!],
  }));

  it('every typing-exact molecule is within the term tolerance', () => {
    const failures: string[] = [];
    let nChecked = 0;

    for (const { mol, typed: t, refTypes } of typed) {
      if (!refTypes || refTypes.length !== mol.atoms.length) continue;
      if (t.atom_types.some((ty, i) => ty !== refTypes[i])) continue;
      const ref = refEnergies.get(mol.name!);
      if (!ref) continue;
      nChecked++;

      const got = calc_energy(t);
      const excluded = EXCLUDED[mol.name!] ?? [];
      const exceptions = EXCEPTIONS[mol.name!] ?? {};

      for (const term of TERM_ORDER) {
        if (excluded.includes(term)) continue;
        const our = got[term];
        const rk = termKey[term] as keyof typeof ref;
        const refVal = ref[rk];
        const d = Math.abs(our - refVal);
        const ex = exceptions[term];
        const tol = ex?.tol ?? DEFAULT_TOL;
        if (d > tol) {
          failures.push(
            `  ${mol.name} ${term}: |${our.toFixed(5)} - ${refVal.toFixed(5)}| = ${d.toExponential(2)} > ${tol} (${ex?.reason ?? 'gate'})`,
          );
        }
      }
    }

    // The whole suite is typing-exact; the gate must cover all of it.
    expect(nChecked).toBe(molecules.length);
    if (failures.length > 0) {
      console.log(`\nPer-term regression gate FAILED (${failures.length} row(s)):\n` + failures.join('\n'));
    }
    expect(failures).toEqual([]);
  });
});

// Pinned-molecule regression rows: the ERULE fragments are the molecules
// that arbitrated the both-pilp torsion fix (2026-08-13). Their torsion
// totals sit at a stationary point where the suite comparison is
// structurally blind — every reading is green there — so a suite-wide
// 1e-4 gate cannot catch a rule change that affects them. These pins do:
// each asserts the torsion total against BatchMin's printed value.
// Tolerances reflect the reference's 5-decimal print precision.
describe('ERULE torsion regression pins (both-pilp arbitration)', () => {
  const cases: [string, number, number][] = [
    // [code, reference torsion, tolerance]
    ['ERULE_01', -2.41459, 5e-5],
    ['ERULE_02', -0.86876, 5e-5],
    ['ERULE_04',  0.43392, 5e-5],
    ['ERULE_08',  7.52598, 5e-5],
  ];

  for (const [code, refTorsion, tol] of cases) {
    it(`${code} torsion matches the reference (${refTorsion})`, () => {
      const mol = molecules.find(m => m.name === code)!;
      const got = calc_energy(assign_atom_types(mol));
      expect(Math.abs(got.torsion - refTorsion)).toBeLessThan(tol);
    });
  }
});
