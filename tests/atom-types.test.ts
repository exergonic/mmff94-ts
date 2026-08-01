import { describe, it, expect } from 'vitest';
import { parse_sdf } from '../src/sdf';
import { assign_atom_types } from '../src/mmff94/atom-types';
import { readFileSync, readdirSync } from 'fs';
import { join } from 'path';

// Expected types from OpenBabel for our benchmark SDFs
const EXPECTED_TYPES: Record<string, number[]> = {
  'methane.sdf':    [1, 5, 5, 5, 5],
  'ethane.sdf':     [1, 1, 5, 5, 5, 5, 5, 5],
  'propane.sdf':    [1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5],
  'butane.sdf':     [1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
  'cyclohexane.sdf': [1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
  'ethene.sdf':     [2, 2, 5, 5, 5, 5],
  'benzene.sdf':    [37, 37, 37, 37, 37, 37, 5, 5, 5, 5, 5, 5],
  'formaldehyde.sdf': [3, 7, 5, 5],
  'water.sdf':      [70, 31, 31],   // MMFF94 dedicated water types: O = 70, H = 31
};

const SDF_DIR = join(__dirname, '..', 'tests', 'fixtures', 'sdf');

describe('Atom typing against OpenBabel reference', () => {
  for (const [filename, expected] of Object.entries(EXPECTED_TYPES)) {
    it(`assigns correct MMFF94 types for ${filename}`, () => {
      const sdf = readFileSync(join(SDF_DIR, filename), 'utf-8');
      const mol = parse_sdf(sdf);
      const typed = assign_atom_types(mol);
      expect(typed.atom_types).toEqual(expected);
    });
  }
});
