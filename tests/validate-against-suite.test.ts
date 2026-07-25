/**
 * Validation against Halgren's MMFF94 suite.
 *
 * Reads molecules from MMFF94.mmd (pre-assigned OPTIMOL types), computes our
 * energies, and compares per-component against BatchMin 5.5 references from
 * MMFF94_bmin.log.
 *
 * Currently informational: many terms are stubs (electrostatic, oop) and some
 * atom types fall outside our current parameter coverage. As terms are
 * implemented, tighten these tolerances toward the 0.01 kcal/mol target.
 */

import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { parse_mmd } from '../src/utils/mmd-parser';
import { calc_energy } from '../src/mmff94/energy/total';

const suiteDir = join(__dirname, 'fixtures', 'validation-suite');

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
      const t = lines[++i].match(/Total Energy\s*=\s*([-0-9.]+)/);
      if (t) e.total = parseFloat(t[1]);
      result.set(currentCode, e);
      currentCode = '';
    }
  }
  return result;
}

describe('Halgren MMFF94 suite validation', () => {
  const molecules = parse_mmd(readFileSync(join(suiteDir, 'MMFF94.mmd'), 'utf-8'));
  const refEnergies = parse_bmin_log(readFileSync(join(suiteDir, 'MMFF94_bmin.log'), 'utf-8'));

  it('parses the full suite', () => {
    expect(molecules.length).toBeGreaterThan(700);
    expect(refEnergies.size).toBeGreaterThan(700);
    expect(molecules[0].name).toBe('AGLYSL01');
  });

  it('reports per-component comparison for AGLYSL01', () => {
    const mol = molecules.find(m => m.name === 'AGLYSL01')!;
    const ref = refEnergies.get('AGLYSL01')!;
    const got = calc_energy(mol);

    const lines = [
      `\nAGLYSL01 — ammonium glycinium sulfate`,
      `  Term            Reference     Computed      Status`,
      `  ───────────────────────────────────────────────────`,
      `  Bond stretch    ${ref.stretch.toFixed(5).padStart(10)}  ${got.bond_stretch.toFixed(5).padStart(10)}  ${Math.abs(got.bond_stretch - ref.stretch) < 0.1 ? '✓' : '✗'}`,
      `  Angle bend      ${ref.bend.toFixed(5).padStart(10)}  ${got.angle_bend.toFixed(5).padStart(10)}  ${Math.abs(got.angle_bend - ref.bend) < 0.1 ? '✓' : '✗'}`,
      `  Stretch-bend    ${ref.strbnd.toFixed(5).padStart(10)}  ${got.stretch_bend.toFixed(5).padStart(10)}  ${Math.abs(got.stretch_bend - ref.strbnd) < 0.1 ? '✓' : '✗'}`,
      `  Torsion         ${ref.torsion.toFixed(5).padStart(10)}  ${got.torsion.toFixed(5).padStart(10)}  ${Math.abs(got.torsion - ref.torsion) < 0.1 ? '✓' : '✗'}`,
      `  VDW             ${ref.vdw.toFixed(5).padStart(10)}  ${got.van_der_waals.toFixed(5).padStart(10)}  ${Math.abs(got.van_der_waals - ref.vdw) < 0.1 ? '✓' : '✗'}`,
      `  Electrostatic   ${ref.elec.toFixed(5).padStart(10)}  ${got.electrostatic.toFixed(5).padStart(10)}  STUB`,
      `  OOP             ${ref.oop.toFixed(5).padStart(10)}  ${got.out_of_plane.toFixed(5).padStart(10)}  STUB`,
      `  ───────────────────────────────────────────────────`,
      `  TOTAL           ${ref.total.toFixed(5).padStart(10)}  ${got.total.toFixed(5).padStart(10)}`,
    ];
    console.log(lines.join('\n'));
  });
});
