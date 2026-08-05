// The part V empirical bond rules (eqs. 18-19, empirical.ts), pinned
// against hand-computed values. The suite exercises the generation
// exactly once (OHMW1's hydroxide O–H, residual 1.4e-6); these tests
// pin the formulas directly, including the measured plain form (no
// δ = 0.008, no mltb/BOij corrections) and the Table V inverse-sixth
// force-constant rule.
import { describe, it, expect } from 'vitest';
import {
  COVALENT_RADII,
  ELECTRONEGATIVITY,
  empirical_bond_length,
  empirical_bond_parameters,
} from '../src/mmff94/parameters/empirical';
import type { Atom } from '../src/types';

const A = (element: string): Atom => ({ index: 0, element, x: 0, y: 0, z: 0 });

describe('eq. (18) — empirical reference bond length', () => {
  it('the OHMW1 O–H: r0 = 0.72 + 0.33 − 0.05·1.3^1.4 = 0.977805', () => {
    const r0 = empirical_bond_length(A('O'), A('H'))!;
    const expected = 0.72 + 0.33 - 0.05 * Math.pow(1.3, 1.4);
    expect(r0).toBeCloseTo(expected, 9);
    expect(r0).toBeCloseTo(0.977805, 5);
  });

  it('the C–C pair uses c = 0.085 (no H) and the same-χ radius sum: 1.54', () => {
    const r0 = empirical_bond_length(A('C'), A('C'))!;
    expect(r0).toBeCloseTo(0.77 + 0.77, 9);
  });

  it('a hetero pair with the 0.085 constant: N–O = 0.73 + 0.72 − 0.085·0.43^1.4', () => {
    const r0 = empirical_bond_length(A('N'), A('O'))!;
    const expected = 0.73 + 0.72 - 0.085 * Math.pow(Math.abs(3.07 - 3.5), 1.4);
    expect(r0).toBeCloseTo(expected, 9);
  });

  it('returns undefined when the element data is missing (Fe)', () => {
    expect(COVALENT_RADII.Fe).toBeUndefined();
    expect(empirical_bond_length(A('Fe'), A('O'))).toBeUndefined();
  });
});

describe('eqs. (18)-(19) — empirical bond parameters', () => {
  it('the OHMW1 O–H: k = 9.1·(0.947/0.977808)^6 = 7.5097', () => {
    const p = empirical_bond_parameters(A('O'), A('H'))!;
    const expectedR0 = 0.72 + 0.33 - 0.05 * Math.pow(1.3, 1.4);
    expect(p.r0).toBeCloseTo(expectedR0, 9);
    const expectedK = 9.1 * Math.pow(0.947 / p.r0, 6);
    expect(p.k_b).toBeCloseTo(expectedK, 9);
    expect(p.k_b).toBeCloseTo(7.5097, 3);
  });

  it('the C–C pair: k = 3.8·(1.512/1.54)^6 = 3.4036', () => {
    const p = empirical_bond_parameters(A('C'), A('C'))!;
    expect(p.r0).toBeCloseTo(1.54, 9);
    expect(p.k_b).toBeCloseTo(3.8 * Math.pow(1.512 / 1.54, 6), 9);
    expect(p.k_b).toBeCloseTo(3.4036, 3);
  });

  it('the key is the element pair by ascending atomic number', () => {
    // O–H and H–O give the same parameters.
    const pHO = empirical_bond_parameters(A('H'), A('O'))!;
    const pOH = empirical_bond_parameters(A('O'), A('H'))!;
    expect(pHO.k_b).toBe(pOH.k_b);
    expect(pHO.r0).toBe(pOH.r0);
  });

  it('returns undefined for pairs outside Table V (H–H has no row)', () => {
    expect(empirical_bond_parameters(A('H'), A('H'))).toBeUndefined();
  });

  it('the published tables are used as-is: r(O) = 0.72, χ(O) = 3.5', () => {
    expect(COVALENT_RADII.O).toBe(0.72);
    expect(ELECTRONEGATIVITY.O).toBe(3.5);
  });
});
