/**
 * Unit tests for the SDF parser.
 *
 * Tests parse_sdf() with real SDF blocks to ensure atoms, bonds,
 * and coordinates are read correctly.
 */

import { describe, it, expect } from 'vitest';
import { parse_sdf } from '../src/sdf';

const ETHANE_SDF = `
  -ISIS-  0123456789

  8  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0900    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8900   -0.6300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3400    0.6300    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4300    0.6300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8800   -0.8900    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8800    0.0000   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
M  END
`;

describe('parse_sdf', () => {
  it('parses ethane correctly — 8 atoms, 7 bonds', () => {
    const mol = parse_sdf(ETHANE_SDF);
    expect(mol.atoms).toHaveLength(8);
    expect(mol.bonds).toHaveLength(7);
  });

  it('assigns correct element symbols', () => {
    const mol = parse_sdf(ETHANE_SDF);
    expect(mol.atoms[0].element).toBe('C');
    expect(mol.atoms[1].element).toBe('C');
    expect(mol.atoms[2].element).toBe('H');
  });

  it('reads bond orders correctly', () => {
    const mol = parse_sdf(ETHANE_SDF);
    // C–C bond is index 0, should be order 1 (single)
    expect(mol.bonds[0].bond_order).toBe(1);
  });

  it('returns empty molecule for malformed input', () => {
    const mol = parse_sdf('not valid');
    expect(mol.atoms).toHaveLength(0);
    expect(mol.bonds).toHaveLength(0);
  });

  it('returns empty molecule for empty input', () => {
    const mol = parse_sdf('');
    expect(mol.atoms).toHaveLength(0);
    expect(mol.bonds).toHaveLength(0);
  });

  it('remaps bond endpoints across a skipped atom line, drops bonds to the skipped atom', () => {
    // Ethane with atom line 3 corrupted: the atom is skipped and every
    // later atom compacts up one slot, while bond lines still address
    // the DECLARED positions. The parser must translate through the
    // declared→compacted map: "1 4" keeps pointing at the same
    // physical H (now array index 2, not 3), and a bond naming the
    // skipped atom ("3 5") is dropped rather than re-pointed.
    const BROKEN = `broken
  test

  8  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  GARBAGE-LINE-SHORT
   -0.8900   -0.6300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3400    0.6300    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4300    0.6300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8800   -0.8900    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8800    0.0000   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  4  1  0
  3  5  1  0
  2  6  1  0
M  END
`;
    const mol = parse_sdf(BROKEN);
    expect(mol.atoms).toHaveLength(7); // atom line 3 skipped, indices compacted
    // "1 2" → (0,1); "1 4" → the H at (−0.89, −0.63, 0) is now index 2;
    // "2 6" → that H is now index 4; "3 5" names the skipped atom → gone.
    // (The pre-fix parser emitted raw declared indices — C1 bonded to
    // whichever H happened to land in the shifted slot.)
    expect(mol.bonds.map(b => [b.atom1, b.atom2]).sort((a, b) => a[0] - b[0] || a[1] - b[1]))
      .toEqual([[0, 1], [0, 2], [1, 4]]);
  });
});
