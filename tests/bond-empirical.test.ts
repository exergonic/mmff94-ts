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
  HL_BADGER_D,
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

  it('a hetero pair with the 0.085 constant: N–O = 0.73 + 0.72 − 0.085·|χ_N−χ_O|^1.4', () => {
    const r0 = empirical_bond_length(A('N'), A('O'))!;
    const expected = 0.73 + 0.72 - 0.085 * Math.pow(Math.abs(ELECTRONEGATIVITY.N - ELECTRONEGATIVITY.O), 1.4);
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

describe("Badger's-rule fallback (the paper's 'should a case arise')", () => {
  it('the d_ij table reproduces the E94 rows of Table V (k = 1.86/(r − d)³)', () => {
    // The E94 force constants were generated from the tabulated
    // lengths by exactly this relation; each row pins its own d.
    const cases: [number, number, number, number][] = [
      // [z1, z2, r0_ref, k_ref] from mmffbndk.par (E94 rows)
      [1, 9, 0.92, 10.6],    // H–F
      [1, 17, 1.28, 4.3],    // H–Cl
      [14, 14, 2.32, 1.3],   // Si–Si
      [16, 53, 2.4, 1.7],    // S–I
      [17, 17, 1.99, 3.5],   // Cl–Cl
      [35, 35, 2.28, 2.4],   // Br–Br
    ];
    for (const [z1, z2, r0, k] of cases) {
      const key = z1 < z2 ? `${z1}-${z2}` : `${z2}-${z1}`;
      const d = HL_BADGER_D[key];
      expect(d).toBeDefined();
      const kBadger = 1.86 / Math.pow(r0 - d!, 3);
      // k_ref is printed at 1 decimal in the par file.
      expect(Math.abs(kBadger - k)).toBeLessThan(0.06);
    }
  });

  it('B–O (outside Table V): r0 = 0.81 + 0.72 − 0.085·1.49^1.4, k = 1.86/(r0 − d(1,1))³', () => {
    const p = empirical_bond_parameters(A('B'), A('O'))!;
    const expectedR0 = 0.81 + 0.72 - 0.085 * Math.pow(Math.abs(2.01 - 3.5), 1.4);
    expect(p.r0).toBeCloseTo(expectedR0, 9);
    expect(p.r0).toBeCloseTo(1.38145, 4);
    // B and O are both row 1 (ELEMENT_ROW): d = 0.679
    expect(p.k_b).toBeCloseTo(1.86 / Math.pow(expectedR0 - 0.679, 3), 9);
    expect(p.k_b).toBeCloseTo(5.3663, 3);
  });

  it('B–N: the row-pair fallback with the same d(1,1), different length', () => {
    const p = empirical_bond_parameters(A('B'), A('N'))!;
    const expectedR0 = 0.81 + 0.73 - 0.085 * Math.pow(Math.abs(2.01 - ELECTRONEGATIVITY.N), 1.4);
    expect(p.r0).toBeCloseTo(expectedR0, 9);
    expect(p.k_b).toBeCloseTo(1.86 / Math.pow(expectedR0 - 0.679, 3), 9);
  });

  it('the fallback key is element-symmetric (B–O ≡ O–B)', () => {
    const pBO = empirical_bond_parameters(A('B'), A('O'))!;
    const pOB = empirical_bond_parameters(A('O'), A('B'))!;
    expect(pBO.k_b).toBe(pOB.k_b);
    expect(pBO.r0).toBe(pOB.r0);
  });

  it('pairs with no element data at all still fall through (Fe–O)', () => {
    expect(empirical_bond_parameters(A('Fe'), A('O'))).toBeUndefined();
  });
});
