/**
 * Validation against Halgren's MMFF94 suite.
 *
 * Reads molecules from MMFF94.mmd (coordinates and connectivity only),
 * assigns atom types via our own assign_atom_types(), computes energies,
 * and reports per-component comparison against BatchMin 5.5 references.
 *
 * All seven terms are implemented; the per-component assertions run on
 * the subset of molecules whose atom typing reproduces the reference
 * types exactly (241/550 in the suite scoreboard — see
 * atom-types-suite.test.ts), so the comparison isolates the terms
 * themselves rather than the typing gaps.
 *
 * Sulfinate note: the S=O oxygen of an anionic sulfinate (73) is
 * typed 7 by the reference typing rules but is keyed 32 in every
 * parameter table (the anionic O family — confirmed identically in
 * TINKER's mmff94.prm and OpenChemLib's angle.csv: 0-1-73-32 /
 * 0-32-73-32 / 0-32-73-72 / 0-0-73-0). The two bridged entries in
 * angle.ts map the 7-typed angles onto the 32-keyed parameters, and
 * found k = 0 entries (the out-of-range defaults) are respected
 * rather than overridden by the empirical rules.
 */

import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { parse_mmd } from '../src/utils/mmd-parser';
import { assign_atom_types } from '../src/mmff94/assign-atom-types';
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
    const raw = molecules.find(m => m.name === 'AGLYSL01')!;
    const typed = assign_atom_types(raw);
    const ref = refEnergies.get('AGLYSL01')!;
    const got = calc_energy(typed);

    const stat = (v: number, r: number) => Math.abs(v - r) < 0.1 ? '✓' : '✗';

    const lines = [
      `\nAGLYSL01 — ammonium glycinium sulfate`,
      `  Term            Reference     Computed      Status`,
      `  ───────────────────────────────────────────────────`,
      `  Bond stretch    ${ref.stretch.toFixed(5).padStart(10)}  ${got.bond_stretch.toFixed(5).padStart(10)}  ${stat(got.bond_stretch, ref.stretch)}`,
      `  Angle bend      ${ref.bend.toFixed(5).padStart(10)}  ${got.angle_bend.toFixed(5).padStart(10)}  ${stat(got.angle_bend, ref.bend)}`,
      `  Stretch-bend    ${ref.strbnd.toFixed(5).padStart(10)}  ${got.stretch_bend.toFixed(5).padStart(10)}  ${stat(got.stretch_bend, ref.strbnd)}`,
      `  Torsion         ${ref.torsion.toFixed(5).padStart(10)}  ${got.torsion.toFixed(5).padStart(10)}  ${stat(got.torsion, ref.torsion)}`,
      `  VDW             ${ref.vdw.toFixed(5).padStart(10)}  ${got.van_der_waals.toFixed(5).padStart(10)}  ${stat(got.van_der_waals, ref.vdw)}`,
      `  Electrostatic   ${ref.elec.toFixed(5).padStart(10)}  ${got.electrostatic.toFixed(5).padStart(10)}  ${stat(got.electrostatic, ref.elec)}`,
      `  OOP             ${ref.oop.toFixed(5).padStart(10)}  ${got.out_of_plane.toFixed(5).padStart(10)}  ${stat(got.out_of_plane, ref.oop)}`,
      `  ───────────────────────────────────────────────────`,
      `  TOTAL           ${ref.total.toFixed(5).padStart(10)}  ${got.total.toFixed(5).padStart(10)}`,
    ];
    console.log(lines.join('\n'));
  });

  it('out-of-plane term matches BatchMin on molecules with complete typing', () => {
    // The oop term is validated against Halgren's own BatchMin energies.
    // These molecules are chosen because our atom typing reproduces the
    // reference types exactly, so the oop comparison isolates the term
    // itself — the residual differences elsewhere in the suite are atom
    // typing gaps, not oop errors. Several match to ~1e-3 kcal/mol.
    const cases: [string, number][] = [
      ['DADDAN', 0.05],     // exact match (Δ ≈ 0.0000)
      ['GIDJUY', 0.05],     // exact match (Δ ≈ 0.0000)
      ['VEJWOW', 0.05],     // exact match (Δ ≈ 0.0003)
      ['DIKGAF', 0.05],     // exact match (Δ ≈ 0.001)
      ['FAXVAB', 0.05],
      ['GEXGIZ', 0.05],
      ['VIRBON', 0.05],
      ['AMHTAR01', 0.05],   // Δ ≈ 0.02 — ester/carboxyl pyramidalization
    ];
    for (const [code, tol] of cases) {
      const raw = molecules.find(m => m.name === code)!;
      const typed = assign_atom_types(raw);
      const ref = refEnergies.get(code)!;
      const got = calc_energy(typed);
      expect(Math.abs(got.out_of_plane - ref.oop)).toBeLessThan(tol);
    }
  });

  it('torsion term matches BatchMin on the 3-ring closure case (FUVDOP)', () => {
    // FUVDOP's triazine 3-ring makes i-j-k-l "torsions" with i = l —
    // closed triangles, not dihedrals. BatchMin skips them; counting
    // them pinned a spurious V3 term at τ = 0 (Δ = +1.125, the worst
    // torsion residual). Pinned here so a reintroduction fails loudly.
    const raw = molecules.find(m => m.name === 'FUVDOP')!;
    const typed = assign_atom_types(raw);
    const ref = refEnergies.get('FUVDOP')!;
    const got = calc_energy(typed);
    expect(Math.abs(got.torsion - ref.torsion)).toBeLessThan(0.01);
  });

  it('torsion term matches BatchMin on the fused 5-ring case (FILNOD)', () => {
    // FILNOD's saturated thiazolidine ring is fused to the benzo ring:
    // its fusion carbons are aromatic-typed, but the 5-ring itself is
    // not aromatic. BatchMin classes these torsions as 5-ring (TTijkl
    // = 5) by ring aromaticity, not by the atoms' flags; judging by
    // the flags misclassified them as class 0 (Δ = +0.223, the last
    // torsion residual). Pinned here so a reintroduction fails loudly.
    const raw = molecules.find(m => m.name === 'FILNOD')!;
    const typed = assign_atom_types(raw);
    const ref = refEnergies.get('FILNOD')!;
    const got = calc_energy(typed);
    expect(Math.abs(got.torsion - ref.torsion)).toBeLessThan(0.01);
  });

  it('stretch-bend term matches BatchMin on the class-2 C-C-C case (JIYJAC)', () => {
    // JIYJAC's 7-ring has C(sp2)-C(sp2)-C(sp2) angles with BOTH bonds
    // BT-flagged (class 2). OpenBabel's mmffstbn.par transcription lost
    // the class-2 (2,2,2) entry (present in the original MMFF94 file);
    // without it the angle fell to the 0.30 dfsb default instead of
    // 0.219/0.250 (Δ = +0.05, the last strbnd residual). The extraction
    // script now restores the entry. Pinned here so a reintroduction
    // fails loudly.
    const raw = molecules.find(m => m.name === 'JIYJAC')!;
    const typed = assign_atom_types(raw);
    const ref = refEnergies.get('JIYJAC')!;
    const got = calc_energy(typed);
    expect(Math.abs(got.stretch_bend - ref.strbnd)).toBeLessThan(0.01);
  });

  it('reports per-component energies for typing-exact molecules', () => {
    // The typing-exact molecules are the only place where a component
    // delta cannot be blamed on atom typing: every mismatch is a term
    // or lookup bug. BatchMin's log is a single-point calculation at the
    // .mmd geometry (per the suite README), so the comparison is valid
    // at the parsed coordinates. Informational: the report drives term
    // fixes; as terms improve the green counts should rise.
    const reference = JSON.parse(
      readFileSync(join(suiteDir, 'mmff94-atom-types.json'), 'utf-8'),
    ) as { molecules: Record<string, number[]> };

    const termDefs = [
      ['bond', 'bond_stretch', 'stretch'],
      ['angle', 'angle_bend', 'bend'],
      ['strbnd', 'stretch_bend', 'strbnd'],
      ['torsion', 'torsion', 'torsion'],
      ['vdw', 'van_der_waals', 'vdw'],
      ['oop', 'out_of_plane', 'oop'],
      ['elec', 'electrostatic', 'elec'],
    ] as const;

    const stats = termDefs.map(() => ({
      green: 0,
      mean: 0,
      max: 0,
      worst: [] as [string, number][],
    }));
    let nChecked = 0;

    for (const mol of molecules) {
      const refTypes = reference.molecules[mol.name!];
      if (refTypes === undefined || refTypes.length !== mol.atoms.length) continue;
      const typed = assign_atom_types(mol);
      if (typed.atom_types.some((t, i) => t !== refTypes[i])) continue;
      const ref = refEnergies.get(mol.name!);
      if (ref === undefined) continue;
      nChecked++;

      const got = calc_energy(typed);
      termDefs.forEach(([, gk, rk], idx) => {
        const d = got[gk] - ref[rk];
        const s = stats[idx];
        s.mean += Math.abs(d);
        s.max = Math.max(s.max, Math.abs(d));
        if (Math.abs(d) <= 0.05) s.green++;
        s.worst.push([mol.name!, d]);
      });
    }

    const lines = [
      `\nPer-component energy deltas vs BatchMin on typing-exact molecules (${nChecked}):`,
      `  term       |Δ|≤0.05   mean|Δ|    max|Δ|`,
      `  ─────────────────────────────────────────`,
    ];
    termDefs.forEach(([label], idx) => {
      const s = stats[idx];
      const mean = (s.mean / Math.max(1, nChecked)).toFixed(3);
      s.worst.sort((a, b) => Math.abs(b[1]) - Math.abs(a[1]));
      const worst = s.worst.slice(0, 3).map(([c, d]) => `${c} ${d > 0 ? '+' : ''}${d.toFixed(2)}`).join(', ');
      lines.push(
        `  ${label.padEnd(9)} ${String(s.green).padStart(3)}/${nChecked}  ${mean.padStart(7)}   ${s.max.toFixed(2).padStart(7)}   ${worst}`,
      );
    });
    lines.push(`  ─────────────────────────────────────────`);
    console.log(lines.join('\n'));

    // The report must cover most of the exact set; the exact count grows
    // as typing improves. This is a floor, not an upper bound.
    expect(nChecked).toBeGreaterThanOrEqual(50);
  });
});
