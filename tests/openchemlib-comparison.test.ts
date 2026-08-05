/**
 * Cross-validation against openchemlib-js — the only other pure-JS MMFF94.
 *
 * openchemlib-js is a GWT transpile of Actelion's Java MMFF94, the force
 * field behind DataWarrior. It implements all seven terms with the correct
 * functional forms (buffered 14-7 vdW, BCI electrostatics, full Fourier
 * torsion), but publishes no numeric validation against the original
 * MMFF94 program. This test runs it on the same Halgren-suite molecules
 * the BatchMin comparison uses, so two independent implementations can be
 * checked against each other and against BatchMin's totals.
 *
 * Charges: the typing reference (OpenBabel's mmd reader) assigns no formal
 * charges from the .mmd, so the molfile handed to openchemlib carries none
 * either — a like-for-like comparison.
 */

import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { Molecule as OCLMolecule, ForceFieldMMFF94, Resources } from 'openchemlib';
import { parse_mmd } from '../src/utils/mmd-parser';
import { assign_atom_types } from '../src/mmff94/assign-atom-types';
import { assign_bci_charges } from '../src/mmff94/charges';
import { calc_energy } from '../src/mmff94/energy/total';
import type { Molecule } from '../src/types';

// openchemlib's force-field parameter tables ship as dist/resources.json
// (the Actelion Java CSVs exported through GWT). They must be registered
// before the first ForceFieldMMFF94 is constructed.
Resources.registerFromNodejs();

const suiteDir = join(__dirname, 'fixtures', 'validation-suite');

function parse_bmin_totals(text: string): Map<string, number> {
  const result = new Map<string, number>();
  let currentCode = '';
  for (const line of text.split('\n')) {
    const m = line.match(/\[\s*(\w+),/);
    if (m) { currentCode = m[1]; continue; }
    if (!currentCode) continue;
    const t = line.match(/Total Energy\s*=\s*([-0-9.]+)/);
    if (t) { result.set(currentCode, parseFloat(t[1])); currentCode = ''; }
  }
  return result;
}

/**
 * Build a V2000 molfile from a parsed suite molecule.
 *
 * Two openchemlib parser quirks matter here: the second header line must
 * be non-empty (a blank one makes fromMolfile return an empty molecule),
 * and the counts line must carry the 3D flag so coordinates are kept.
 * Formal charges are deliberately absent — see the file header.
 */
function to_molfile(mol: Molecule): string {
  const lines: string[] = [];
  lines.push(mol.name ?? 'molecule');
  lines.push('  mmff94-ts cross-check');
  lines.push('');
  const n = mol.atoms.length;
  lines.push(`${String(n).padStart(3)}${String(mol.bonds.length).padStart(3)}  0  0  1  0            999 V2000`);
  for (const a of mol.atoms) {
    lines.push(
      `${a.x.toFixed(4).padStart(10)}${a.y.toFixed(4).padStart(10)}${a.z.toFixed(4).padStart(10)} ${a.element.padStart(3)} 0  0  0  0  0  0  0  0  0  0  0  0`,
    );
  }
  for (const b of mol.bonds) {
    lines.push(`${String(b.atom1 + 1).padStart(3)}${String(b.atom2 + 1).padStart(3)}${String(b.bond_order).padStart(3)}  0  0  0  0`);
  }
  lines.push('M  END');
  return lines.join('\n');
}

function our_total(mol: Molecule): number {
  const typed = assign_atom_types(mol);
  return calc_energy(assign_bci_charges(typed)).total;
}

function ocl_total(mol: Molecule): number | string {
  try {
    const omol = OCLMolecule.fromMolfile(to_molfile(mol));
    const ff = new ForceFieldMMFF94(omol, ForceFieldMMFF94.MMFF94);
    return ff.getTotalEnergy();
  } catch (e) {
    // GWT-transpiled exceptions surface as "Class$S63: message" — the
    // class prefix is internal mangling, the message is the useful part.
    return (e as Error).message.replace(/^Class\$S\d+: /, '').slice(0, 80);
  }
}

interface Row {
  code: string;
  our: number;
  ocl: number | string;
  bmin: number;
}

function compute_rows(): Row[] {
  const molecules = parse_mmd(readFileSync(join(suiteDir, 'MMFF94.mmd'), 'utf-8'));
  const refEnergies = parse_bmin_totals(readFileSync(join(suiteDir, 'MMFF94_bmin.log'), 'utf-8'));
  const reference = JSON.parse(
    readFileSync(join(suiteDir, 'mmff94-atom-types.json'), 'utf-8'),
  ) as { molecules: Record<string, number[]> };

  const rows: Row[] = [];
  for (const mol of molecules) {
    const refTypes = reference.molecules[mol.name!];
    if (refTypes === undefined || refTypes.length !== mol.atoms.length) continue;
    const typed = assign_atom_types(mol);
    if (typed.atom_types.some((t, i) => t !== refTypes[i])) continue;
    const bmin = refEnergies.get(mol.name!);
    if (bmin === undefined) continue;
    rows.push({ code: mol.name!, our: our_total(mol), ocl: ocl_total(mol), bmin });
  }
  return rows;
}

function stat(values: number[]): { mean: number; max: number; le01: number; le05: number; le10: number } {
  let mean = 0;
  let max = 0;
  let le01 = 0;
  let le05 = 0;
  let le10 = 0;
  for (const v of values) {
    const d = Math.abs(v);
    mean += d;
    max = Math.max(max, d);
    if (d <= 0.1) le01++;
    if (d <= 0.5) le05++;
    if (d <= 1.0) le10++;
  }
  return { mean: mean / values.length, max, le01, le05, le10 };
}

describe('openchemlib-js cross-validation', () => {
  const rows = compute_rows();
  const oclOk = rows.filter(r => typeof r.ocl === 'number') as { code: string; our: number; ocl: number; bmin: number }[];
  const oclFailed = rows.filter(r => typeof r.ocl === 'string') as { code: string; ocl: string }[];

  it('parses the suite and produces BatchMin totals', () => {
    expect(rows.length).toBeGreaterThanOrEqual(50);
  });

  it('reports total-energy comparison on typing-exact molecules', () => {
    const lines = [
      `\nTotal energies (kcal/mol) vs BatchMin on typing-exact molecules (${rows.length}):`,
      `  code         ours      openchemlib    BatchMin   Δours   Δocl`,
      `  ─────────────────────────────────────────────────────────────`,
    ];
    for (const r of rows) {
      const dOurs = (r.our - r.bmin).toFixed(3);
      const dOcl = typeof r.ocl === 'number' ? (r.ocl - r.bmin).toFixed(3) : '  n/a';
      lines.push(
        `  ${r.code.padEnd(9)} ${r.our.toFixed(5).padStart(10)} ${(typeof r.ocl === 'number' ? r.ocl : NaN).toFixed(5).padStart(12)} ${r.bmin.toFixed(5).padStart(11)}  ${dOurs.padStart(6)} ${dOcl.padStart(7)}`,
      );
    }
    lines.push(`  ─────────────────────────────────────────────────────────────`);

    const ourStat = stat(oclOk.map(r => r.our - r.bmin));
    const oclStat = stat(oclOk.map(r => r.ocl - r.bmin));
    const both = oclOk.filter(r => Math.abs(r.our - r.bmin) <= 0.1 && Math.abs(r.ocl - r.bmin) <= 0.1).length;

    // Disagreement candidates: one implementation matches BatchMin while
    // the other is off by more than half a kcal/mol. Where openchemlib
    // agrees with BatchMin but we do not, that is our bug; the reverse
    // direction is an openchemlib quirk worth documenting.
    const ourBug = oclOk.filter(r => Math.abs(r.our - r.bmin) > 0.5 && Math.abs(r.ocl - r.bmin) <= 0.1);
    const oclQuirk = oclOk.filter(r => Math.abs(r.ocl - r.bmin) > 0.5 && Math.abs(r.our - r.bmin) <= 0.1);

    lines.push(
      `  ours:       mean|Δ| ${ourStat.mean.toFixed(3)}  max|Δ| ${ourStat.max.toFixed(3)}  |Δ|≤0.1: ${ourStat.le01}/${oclOk.length}`,
      `  openchemlib: mean|Δ| ${oclStat.mean.toFixed(3)}  max|Δ| ${oclStat.max.toFixed(3)}  |Δ|≤0.1: ${oclStat.le01}/${oclOk.length}`,
      `  both within 0.1 of BatchMin: ${both}/${oclOk.length}`,
      `  openchemlib failed on ${oclFailed.length}: ${oclFailed.slice(0, 5).map(r => `${r.code} (${r.ocl})`).join(', ')}`,
    );
    if (ourBug.length > 0) {
      lines.push(`  OUR BUG CANDIDATES (ocl≈BatchMin, we differ): ${ourBug.slice(0, 6).map(r => `${r.code} ${(r.our - r.bmin).toFixed(2)}`).join(', ')}`);
    }
    if (oclQuirk.length > 0) {
      lines.push(`  OCL QUIRK CANDIDATES (we≈BatchMin, ocl differs): ${oclQuirk.slice(0, 6).map(r => `${r.code} ${(r.ocl - r.bmin).toFixed(2)}`).join(', ')}`);
    }
    lines.push(
      `  ─ note: BAOXLM01/CAMALD03 are the documented formal-charge salts (elec 89/91` +
      ` in the suite test); the .mmd carries the reference partial charges and neither` +
      ` implementation feeds them. FUVDOP's and FILNOD's torsion residuals and` +
      ` JIYJAC's strbnd residual were real bugs — degenerate i = l "torsions" in` +
      ` a 3-ring, 5-ring torsions classed by atom flags instead of ring aromaticity,` +
      ` and an entry OpenBabel's strbnd transcription lost — all fixed (torsion` +
      ` 91/91, strbnd 91/91).`,
    );
    console.log(lines.join('\n'));
  });

  it('reports how much of the full suite openchemlib can type at all', () => {
    // The typing-exact set above is where our typing is proven right. The
    // full suite is the harder question: of 753 molecules, how many can the
    // only other JS MMFF94 even construct a force field for?
    const molecules = parse_mmd(readFileSync(join(suiteDir, 'MMFF94.mmd'), 'utf-8'));
    let accepted = 0;
    const rejected: [string, string][] = [];
    for (const mol of molecules) {
      const e = ocl_total(mol);
      if (typeof e === 'number') accepted++;
      else if (rejected.length < 8) rejected.push([mol.name!, e]);
    }
    console.log(
      `\nopenchemlib full-suite acceptance: ${accepted}/${molecules.length}` +
      `\n  rejected examples: ${rejected.map(([c, m]) => `${c} (${m})`).join('\n  ')}`,
    );
    expect(accepted).toBeGreaterThanOrEqual(400);
  }, 120000);
});
