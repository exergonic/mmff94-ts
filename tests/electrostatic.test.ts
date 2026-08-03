/**
 * Electrostatic (Coulombic) energy — eq. (6) of MMFF part III.
 *
 * The hand-computed case is a 4-atom chain: only the 1-4 pair (0, 3)
 * is evaluated (1-2 and 1-3 pairs are excluded from MMFF94
 * electrostatics), scaled by 0.75, with the buffering constant
 * S = 0.05 Å added to the distance:
 *
 *   E = 332.0716 · q₀ · q₃ · 0.75 / (r + 0.05)
 */

import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { parse_sdf } from '../src/sdf';
import { assign_atom_types } from '../src/mmff94/atom-types';
import { compute_bci_charges } from '../src/mmff94/charges';
import { calc_electrostatic_energy } from '../src/mmff94/energy/electrostatic';
import { calc_energy } from '../src/mmff94/energy/total';

/** A linear chain of n atoms along x, 1 Å apart, all type 1 (CR). */
function chain(n: number): ReturnType<typeof parse_sdf> {
  return {
    name: 'chain',
    atoms: Array.from({ length: n }, (_, i) => ({
      index: i,
      element: 'C',
      x: i,
      y: 0,
      z: 0,
    })),
    bonds: Array.from({ length: n - 1 }, (_, i) => ({
      atom1: i,
      atom2: i + 1,
      bond_order: 1,
    })),
  };
}

describe('Electrostatic energy', () => {
  it('evaluates only the 1-4 pair of a 4-atom chain, scaled by 0.75', () => {
    const mol = chain(4);
    const typed = { ...mol, atom_types: [1, 1, 1, 1], partial_charges: [1, 0, 0, -1] };

    // Pairs: (0,1) 1-2, (1,2) 1-2, (2,3) 1-2, (0,2) 1-3, (1,3) 1-3 — all
    // excluded. Only (0,3): 1-4, r = 3 Å, ×0.75, buffered by +0.05.
    const expected = (332.0716 * 1 * -1 * 0.75) / (3 + 0.05);
    expect(calc_electrostatic_energy(typed)).toBeCloseTo(expected, 10);
  });

  it('includes pairs beyond 1-4 unscaled', () => {
    const mol = chain(5);
    const typed = { ...mol, atom_types: [1, 1, 1, 1, 1], partial_charges: [1, 0, 0, -1, -1] };

    // (0,3): 1-4 ×0.75 at r = 3; (0,4): 1-5 unscaled at r = 4.
    const expected =
      (332.0716 * -0.75) / (3 + 0.05) + (332.0716 * -1) / (4 + 0.05);
    expect(calc_electrostatic_energy(typed)).toBeCloseTo(expected, 10);
  });

  it('computes charges on the fly when none are stored', () => {
    // Ammonia without precomputed charges: the term computes them
    // internally (all pairs 1-2/1-3) and still gives zero — the
    // input molecule stays untouched (pure contract).
    const sdf = readFileSync(join(__dirname, 'fixtures', 'sdf', 'ammonia.sdf'), 'utf-8');
    const typed = assign_atom_types(parse_sdf(sdf));
    expect(typed.partial_charges).toBeUndefined();
    expect(calc_electrostatic_energy(typed)).toBe(0);
    expect(typed.partial_charges).toBeUndefined(); // input never mutated
  });

  it('matches the reference logs for every fixture (with charges)', () => {
    const refs: Record<string, number> = {
      benzene: 3.0781,
      ethene: 8.0530,
      pyridine: 2.0939,
      pyrrole: 3.0720,
      nicotine: -2.2135,
      ammonia: 0.0,
    };
    for (const [name, ref] of Object.entries(refs)) {
      const sdf = readFileSync(join(__dirname, 'fixtures', 'sdf', `${name}.sdf`), 'utf-8');
      const typed = assign_atom_types(parse_sdf(sdf));
      const charged = compute_bci_charges(typed);
      expect(Math.abs(calc_electrostatic_energy(charged) - ref)).toBeLessThan(0.001);
    }
  });

  it('completes the total energy — fixture totals now match exactly', () => {
    // Electrostatic was the last term to land; with it live, every
    // typed fixture's total matches its obenergy log.
    const refs: Record<string, number> = {
      benzene: 16.22697,
      pyridine: 15.5234,
      pyrrole: 3.2868,
      nicotine: 30.2543,
    };
    for (const [name, ref] of Object.entries(refs)) {
      const sdf = readFileSync(join(__dirname, 'fixtures', 'sdf', `${name}.sdf`), 'utf-8');
      const typed = assign_atom_types(parse_sdf(sdf));
      const charged = compute_bci_charges(typed);
      expect(Math.abs(calc_energy(charged).total - ref)).toBeLessThan(0.001);
    }
  });
});
