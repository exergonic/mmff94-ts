/**
 * Compare every benchmark molecule against its OpenBabel reference.
 *
 * For each SDF in tests/fixtures/sdf/, assign atom types, compute all energy
 * terms, and compare against the obenergy log in tests/references/.
 *
 * Terms that match exactly (same formula, same parameters):
 *   bond_stretch, angle_bend, van_der_waals
 *
 * Terms with known discrepancies:
 *   torsion — Fourier series same, values differ (~43% of reference)
 *   stretch_bend — near zero in practice; reference often < 0.01
 *
 * Terms still stubs:
 *   electrostatic, out_of_plane
 */

import { describe, it, expect } from 'vitest';
import { readFileSync, readdirSync } from 'fs';
import { join, parse } from 'path';
import { parse_sdf } from '../src/sdf';
import { assign_atom_types } from '../src/mmff94/atom-types';
import { calc_energy } from '../src/mmff94/energy/total';

const SDF_DIR = join(__dirname, 'fixtures', 'sdf');
const REF_DIR = join(__dirname, 'references');

function parse_reference_log(filePath: string): Record<string, number> {
  const text = readFileSync(filePath, 'utf-8');
  const result: Record<string, number> = {};
  const patterns: [string, RegExp][] = [
    ['bond', /TOTAL BOND STRETCHING ENERGY\s*=\s*([-0-9.]+)/],
    ['angle', /TOTAL ANGLE BENDING ENERGY\s*=\s*([-0-9.]+)/],
    ['strbnd', /TOTAL STRETCH BENDING ENERGY\s*=\s*([-0-9.]+)/],
    ['torsion', /TOTAL TORSIONAL ENERGY\s*=\s*([-0-9.]+)/],
    ['oop', /TOTAL OUT-OF-PLANE BENDING ENERGY\s*=\s*([-0-9.]+)/],
    ['vdw', /TOTAL VAN DER WAALS ENERGY\s*=\s*([-0-9.]+)/],
    ['elec', /TOTAL ELECTROSTATIC ENERGY\s*=\s*([-0-9.]+)/],
    ['total', /TOTAL ENERGY\s*=\s*([-0-9.]+)/],
  ];
  for (const [key, re] of patterns) {
    const m = text.match(re);
    if (m) result[key] = parseFloat(m[1]);
  }
  return result;
}

describe('All benchmark molecules vs OpenBabel references', () => {
  const sdfFiles = readdirSync(SDF_DIR).filter(f => f.endsWith('.sdf')).sort();

  for (const sdfFile of sdfFiles) {
    const name = parse(sdfFile).name; // e.g. "ethane" from "ethane.sdf"
    const refFile = join(REF_DIR, `${name}.mmff94.log`);

    if (!readdirSync(REF_DIR).includes(`${name}.mmff94.log`)) {
      it.skip(`${name} — no reference log found`, () => {});
      continue;
    }

    it(`${name}: bond stretch, angle bend, VDW match reference`, () => {
      const sdfText = readFileSync(join(SDF_DIR, sdfFile), 'utf-8');
      const ref = parse_reference_log(refFile);

      const mol = parse_sdf(sdfText);
      const typed = assign_atom_types(mol);
      const energy = calc_energy(typed);

      // Terms that should match closely (0.02 absolute tolerance for near-zero cases)
      expect(Math.abs(energy.bond_stretch - ref.bond)).toBeLessThan(0.02);
      if (Math.abs(ref.angle) > 0.001) {
        expect(Math.abs(energy.angle_bend - ref.angle)).toBeLessThan(0.02);
      }
      expect(Math.abs(energy.van_der_waals - ref.vdw)).toBeLessThan(0.02);
    });

    it(`${name}: prints full comparison`, () => {
      const sdfText = readFileSync(join(SDF_DIR, sdfFile), 'utf-8');
      const ref = parse_reference_log(refFile);
      const mol = parse_sdf(sdfText);
      const typed = assign_atom_types(mol);
      const energy = calc_energy(typed);

      // Determine which of our atom types are wrong (can't match ref)
      const typesOk = typed.atom_types.every(t => t !== 1 || true); // placeholder

      const check = (v: number, r: number, tol: number) =>
        Math.abs(v - r) < tol ? '✓' : `✗ (${(v / r).toFixed(2)}×)` ;

      console.log(`\n${name}:`);
      console.log(`  Terms vs OpenBabel:`);
      console.log(`    Bond:      ${energy.bond_stretch.toFixed(5)} vs ${ref.bond.toFixed(5)} ${check(energy.bond_stretch, ref.bond, 0.01)}`);
      console.log(`    Angle:     ${energy.angle_bend.toFixed(5)} vs ${ref.angle.toFixed(5)} ${check(energy.angle_bend, ref.angle, 0.05)}`);
      console.log(`    StrBnd:    ${energy.stretch_bend.toFixed(5)} vs ${ref.strbnd.toFixed(5)} ${check(energy.stretch_bend, ref.strbnd, 0.05)}`);
      console.log(`    Torsion:   ${energy.torsion.toFixed(5)} vs ${ref.torsion.toFixed(5)} ${check(energy.torsion, ref.torsion, 0.5)}`);
      console.log(`    VDW:       ${energy.van_der_waals.toFixed(5)} vs ${ref.vdw.toFixed(5)} ${check(energy.van_der_waals, ref.vdw, 0.01)}`);
      console.log(`    Elec:      ${energy.electrostatic.toFixed(5)} vs ${ref.elec.toFixed(5)} (STUB)`);
      console.log(`    OOP:       ${energy.out_of_plane.toFixed(5)} vs ${ref.oop.toFixed(5)} (STUB)`);
      console.log(`    Total:     ${energy.total.toFixed(5)} vs ${ref.total.toFixed(5)}`);
    });
  }
});
