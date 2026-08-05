// Part V empirical angle rules (eq. 20 + the θ₀ default protocol):
// the last-resort fallback when the angle lookup chain misses
// entirely. The suite never exercises it (the par's wildcards always
// resolve, validated 747/747), so these tests pin the formulas against
// the paper's spec — the θ₀ protocol (part V p. 626) and the eq. (20)
// force-constant form with D = ((r̃_ij − r̃_jk)/(r̃_ij + r̃_jk))² (both
// numerator and denominator squared — the paper's form, also Tinker's
// kbond.f; OpenBabel's GetAngleParam omits the numerator square).
import { describe, it, expect } from 'vitest';
import { ATOM_TYPE_PROPERTIES } from '../src/mmff94/parameters';
import { empirical_theta0, empirical_ka } from '../src/mmff94/parameters/empirical';
import type { Atom } from '../src/types';

const A = (element: string): Atom => ({ index: 0, element, x: 0, y: 0, z: 0 });

describe('part V θ₀ default protocol (empirical_theta0)', () => {
  it('starts at 120° and assigns 109.45° for a tetrahedral center', () => {
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[1], 'C', 0)).toBeCloseTo(109.45, 6);
  });

  it('assigns 105° for a divalent oxygen center', () => {
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[6], 'O', 0)).toBeCloseTo(105.0, 6);
  });

  it('assigns 95° for a divalent center with Z > 10 (sulfur)', () => {
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[15], 'S', 0)).toBeCloseTo(95.0, 6);
  });

  it('assigns 180° for a linear center (lin flag)', () => {
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[4], 'C', 0)).toBeCloseTo(180.0, 6);
  });

  it('assigns 107° for a tricoordinate amine nitrogen', () => {
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[8], 'N', 0)).toBeCloseTo(107.0, 6);
  });

  it('assigns 92° for a tricoordinate center that is not nitrogen', () => {
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[26], 'P', 0)).toBeCloseTo(92.0, 6);
  });

  it('does NOT apply the 95° rule outside the crd = 2 branch', () => {
    // A heavy element with crd 1 (halide type) keeps the 120° default —
    // the atomic-number branch lives inside crd = 2.
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[13], 'I', 0)).toBeCloseTo(120.0, 6);
  });

  it('forces 60°/90° for three- and four-membered ring angles', () => {
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[1], 'C', 3)).toBeCloseTo(60.0, 6);
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[1], 'C', 4)).toBeCloseTo(90.0, 6);
    expect(empirical_theta0(ATOM_TYPE_PROPERTIES[1], 'C', 7)).toBeCloseTo(90.0, 6);
  });
});

describe('part V eq. (20) angle-bending force constant (empirical_ka)', () => {
  it('computes the alkane C–C–C value by hand', () => {
    // k = β·Z_C·C_C·Z_C / (rr·θ₀²·rad²), β = 1.75, r0(C–C) = 1.508
    // (the par's 0-1-1 row), θ₀ = 109.47°:
    //   1.75·2.494·1.016·2.494 / (3.016·(109.47·π/180)²) = 1.0045
    const theta0 = 109.47;
    const k = empirical_ka(A('C'), A('C'), A('C'), 1.508, 1.508, theta0, 0);
    const expected =
      (1.75 * 2.494 * 1.016 * 2.494) /
      (3.016 * Math.pow(theta0 * Math.PI / 180, 2));
    expect(k).toBeCloseTo(expected, 9);
    expect(k).toBeCloseTo(1.0045, 3);
  });

  it('applies the ring-scale factors (×0.85 for 4-rings, ×0.05 for 3-rings)', () => {
    const theta0 = 109.47;
    const k0 = empirical_ka(A('C'), A('C'), A('C'), 1.508, 1.508, theta0, 0);
    expect(empirical_ka(A('C'), A('C'), A('C'), 1.508, 1.508, theta0, 4)).toBeCloseTo(0.85 * k0, 9);
    expect(empirical_ka(A('C'), A('C'), A('C'), 1.508, 1.508, theta0, 3)).toBeCloseTo(0.05 * k0, 9);
  });

  it('uses the SQUARED ratio D = ((r0ab−r0bc)/(r0ab+r0bc))² with exp(−2D)', () => {
    // Asymmetric center O–C–C with the par's r0(O–C) = 1.451 and
    // r0(C–C) = 1.508 — D = ((1.451−1.508)/rr)², k scaled by exp(−2D).
    const theta0 = 109.47;
    const k = empirical_ka(A('O'), A('C'), A('C'), 1.451, 1.508, theta0, 0);
    const rr = 1.451 + 1.508;
    const D = ((1.451 - 1.508) * (1.451 - 1.508)) / (rr * rr);
    const expected =
      (1.75 * 3.045 * 1.016 * 2.494 * Math.exp(-2 * D)) /
      (rr * Math.pow(theta0 * Math.PI / 180, 2));
    expect(k).toBeCloseTo(expected, 9);
    // And the exp(−2D) factor is active: recomputing the same product
    // with D = 0 would give a larger k.
    const D0 = 0;
    const k_noD =
      (1.75 * 3.045 * 1.016 * 2.494 * Math.exp(-2 * D0)) /
      (rr * Math.pow(theta0 * Math.PI / 180, 2));
    expect(k).toBeLessThan(k_noD);
  });
});
