import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { parse_mmd } from '../src/utils/mmd-parser';

describe('MMD Parser', () => {
  it('parses first molecule (AGLYSL01)', () => {
    const mmd = readFileSync(
      join(__dirname, 'fixtures', 'validation-suite', 'MMFF94.mmd'), 'utf-8'
    );
    const molecules = parse_mmd(mmd);
    expect(molecules.length).toBeGreaterThan(0);

    const mol = molecules[0];
    expect(mol.name).toBe('AGLYSL01');
    expect(mol.atoms.length).toBe(10);
    expect(mol.bonds.length).toBeGreaterThan(0);
    expect(mol.atom_types.length).toBe(10);
    expect(mol.partial_charges?.length).toBe(10);

    // First atom: should be a carbon (type 3, carbonyl C)
    expect(mol.atom_types[0]).toBe(3);
    expect(mol.atoms[0].element).toBe('C');
  });

  it('parses second molecule (AMHTAR01)', () => {
    const mmd = readFileSync(
      join(__dirname, 'fixtures', 'validation-suite', 'MMFF94.mmd'), 'utf-8'
    );
    const molecules = parse_mmd(mmd);
    expect(molecules.length).toBeGreaterThan(1);

    const mol = molecules[1];
    expect(mol.name).toBe('AMHTAR01');
    expect(mol.atoms.length).toBe(15);
  });

  it('handles empty text', () => {
    expect(parse_mmd('')).toEqual([]);
  });

  it('handles text with no headers', () => {
    expect(parse_mmd('some random text\nwith no structure')).toEqual([]);
  });
});
