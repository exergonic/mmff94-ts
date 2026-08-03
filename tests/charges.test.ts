/**
 * The MMFF94 bond-charge-increment (BCI) partial-charge model.
 *
 * Reference values come from the PARTIAL CHARGES sections of the
 * obenergy logs in tests/references/ — OpenBabel's own MMFF94 output.
 */

import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { parse_sdf } from '../src/sdf';
import { assign_atom_types } from '../src/mmff94/atom-types';
import { compute_bci_charges } from '../src/mmff94/charges';

const SDF_DIR = join(__dirname, 'fixtures', 'sdf');
const REF_DIR = join(__dirname, 'references');

/** The per-atom partial charges from an obenergy log. */
function parse_reference_charges(filePath: string): number[] {
  const text = readFileSync(filePath, 'utf-8');
  const lines = text.split('\n');
  const charges: number[] = [];
  let inPartial = false;
  for (const line of lines) {
    // The obenergy logs letter-space their section headers
    // ("P A R T I A L   C H A R G E S") — match the spaced form.
    if (line.includes('P A R T I A L')) {
      inPartial = true;
      continue;
    }
    if (!inPartial) continue;
    const parts = line.split('\t');
    if (parts.length >= 2 && /^\d+$/.test(parts[0].trim())) {
      charges.push(parseFloat(parts[1]));
    } else if (!/^IDX/.test(parts[0].trim()) && parts[0].trim() !== '') {
      break; // left the charge table
    }
  }
  return charges;
}

describe('BCI partial charges', () => {
  // Molecule → expected per-atom charges, taken from the reference
  // logs (the numbers ARE the OpenBabel BCI sum).
  const cases: Record<string, number[]> = {
    'ammonia.sdf': [-1.08, 0.36, 0.36, 0.36],           // N(8)–H(23): −0.36 × 3
    'water.sdf': [-0.86, 0.43, 0.43],                    // O(70)–H(31): '0-31-70' = −0.43
    'benzene.sdf': [-0.15, -0.15, -0.15, -0.15, -0.15, -0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15],
    'pyridine.sdf': [-0.62, 0.16, -0.15, -0.15, -0.15, 0.16, 0.15, 0.15, 0.15, 0.15, 0.15],
    'formaldehyde.sdf': [0.45, -0.57, 0.06, 0.06],
  };

  for (const [filename, expected] of Object.entries(cases)) {
    it(`matches the reference charges for ${filename}`, () => {
      const sdf = readFileSync(join(SDF_DIR, filename), 'utf-8');
      const typed = assign_atom_types(parse_sdf(sdf));
      // Pure step: the ORIGINAL molecule is untouched; the charges
      // arrive on the returned copy.
      const charged = compute_bci_charges(typed);

      expect(typed.partial_charges).toBeUndefined();
      expect(charged.partial_charges).toBeDefined();
      expect(charged.partial_charges!.length).toBe(expected.length);
      for (let i = 0; i < expected.length; i++) {
        expect(charged.partial_charges![i]).toBeCloseTo(expected[i], 5);
      }
    });
  }

  it('matches the reference log charges for every fixture', () => {
    for (const file of ['ammonia', 'water', 'benzene', 'pyridine', 'formaldehyde']) {
      const sdf = readFileSync(join(SDF_DIR, `${file}.sdf`), 'utf-8');
      const typed = assign_atom_types(parse_sdf(sdf));
      const charged = compute_bci_charges(typed);
      const ref = parse_reference_charges(join(REF_DIR, `${file}.mmff94.log`));
      expect(charged.partial_charges!.length).toBe(ref.length);
      charged.partial_charges!.forEach((q, i) => {
        expect(Math.abs(q - ref[i])).toBeLessThan(1e-4);
      });
    }
  });
});
