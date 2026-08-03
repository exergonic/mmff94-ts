/**
 * The formal-charge model (part V eq. 15) against the suite's own
 * reference partial charges.
 *
 * Every .mmd file carries the OPTIMOL/BatchMin partial charge per atom
 * (the pchg column). For each typing-exact molecule we compute our
 * charges — BCI sum + the primary formal charges q⁰ with their sharing
 * (charges.ts) — and compare per atom.
 *
 * Two molecules are excluded, both thiosulfinate anions (JALSOE,
 * SO18A): BatchMin adjusts their S–S bond to the "MMFF dative
 * representation" (its log says so explicitly), and the reference
 * partial charges there are not reproducible from eq. (15) with the
 * transcribed parameters — the same family that carries the
 * documented angle-bend residual in validate-against-suite.test.ts.
 */

import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { parse_mmd } from '../src/utils/mmd-parser';
import { assign_atom_types } from '../src/mmff94/atom-types';
import { compute_bci_charges } from '../src/mmff94/charges';

const suiteDir = join(__dirname, 'fixtures', 'validation-suite');

/** The per-atom reference partial charges (the .mmd pchg column). */
function reference_charges(mmdText: string, name: string, nAtoms: number): number[] {
  const pchg: number[] = new Array(nAtoms).fill(0);
  let inMol = false;
  for (const line of mmdText.split('\n')) {
    const head = line.match(/^\s*\d+\s+\[(\w+),/);
    if (head) { inMol = head[1] === name; continue; }
    if (!inMol) continue;
    const p = line.trim().split(/\s+/);
    // atom line: bmin_type + 6 neighbor pairs + x y z label idx fchg pchg name subname serial
    if (p.length >= 20) {
      const serial = parseInt(p[p.length - 1], 10);
      if (!isNaN(serial)) pchg[serial - 1] = parseFloat(p[p.length - 4]);
    }
  }
  return pchg;
}

describe('partial charges vs the validation-suite reference (pchg)', () => {
  it('reproduces the reference per-atom charges for every typing-exact molecule', () => {
    const refs = JSON.parse(readFileSync(join(suiteDir, 'mmff94-atom-types.json'), 'utf-8')) as {
      molecules: Record<string, number[]>;
    };
    const mmdText = readFileSync(join(suiteDir, 'MMFF94.mmd'), 'utf-8');
    const molecules = parse_mmd(mmdText);

    // BatchMin's dative-adjusted thiosulfinate pair (see header).
    const DATIVE_EXCLUDED = new Set(['JALSOE', 'SO18A']);

    let checked = 0;
    let worst = 0;
    for (const mol of molecules) {
      const refTypes = refs.molecules[mol.name];
      if (!refTypes || refTypes.length !== mol.atoms.length) continue;
      const typed = assign_atom_types(mol);
      if (!typed.atom_types.every((t, i) => t === refTypes[i])) continue;
      if (DATIVE_EXCLUDED.has(mol.name)) continue;

      compute_bci_charges(typed);
      const ref = reference_charges(mmdText, mol.name, mol.atoms.length);
      checked++;
      for (let i = 0; i < typed.atoms.length; i++) {
        const dev = Math.abs(typed.partial_charges![i] - ref[i]);
        worst = Math.max(worst, dev);
        expect(dev, `${mol.name} atom ${i}: got ${typed.partial_charges![i]}, ref ${ref[i]}`)
          .toBeLessThan(1e-3);
      }
    }
    // The model is exact on the clean set: the worst deviation across
    // all atoms is far below the 1e-3 assertion (measured ~1e-4).
    expect(checked).toBeGreaterThan(100);
    expect(worst).toBeLessThan(1e-3);
  });
});
