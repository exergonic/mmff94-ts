// The part V empirical torsion rules (empirical.ts), pinned against the
// paper (part V pp. 631-632). The suite never exercises the rules (the
// par's wildcards always resolve, validated 747/747), so these tests
// pin the protocol's cases (a)-(h) against hand-computed values —
// including the three paper-arbitrated fixes: rule (c) gates on the
// j-k bond order (not the i-j bond, and not the mltb flags alone),
// rule (g) cases (2)/(3) read the PILP atom's own mltb, and rule (h)
// divides by N_bc = (crd(J)−1)·(crd(K)−1).
//
// Table X constants used below: U_C = U_N = U_O = 2.0 (so √(U·U) = 2),
// V_C = 2.12, V_N = 1.5.
import { describe, it, expect } from 'vitest';
import { ATOM_TYPES, ATOM_TYPE_PROPERTIES } from '../src/mmff94/parameters';
import { empirical_torsion } from '../src/mmff94/parameters/empirical';

const P = (t: number) => ATOM_TYPE_PROPERTIES[t];
const el = (t: number) => ATOM_TYPES[t].element;
const tor = (tj: number, tk: number, order_jk: number, aromatic_jk = false) =>
  empirical_torsion(P(tj), P(tk), el(tj), el(tk), order_jk, aromatic_jk);

describe('rule (a) — linear centers carry no torsion', () => {
  it('skips when j is linear (type 4, alkyne C)', () => {
    expect(tor(4, 1, 1).skip).toBe(true);
  });
});

describe('rule (b) — aromatic central bond → V2', () => {
  it('benzene C–C (37-37): π = 0.5, β = 6 → V2 = 6.0', () => {
    const r = tor(37, 37, 1, true);
    expect(r.skip).toBe(false);
    expect(r.v2).toBeCloseTo(6.0, 9);
    expect(r.v1).toBe(0);
    expect(r.v3).toBe(0);
  });

  it('pyridine C–N (37-38): val 4/3 → β = 3 → V2 = 3.0', () => {
    expect(tor(38, 37, 1, true).v2).toBeCloseTo(3.0, 9);
  });

  it('pyrrole N–C (39-37): one lone pair → π = 0.3, β = 3 → V2 = 1.8', () => {
    expect(tor(39, 37, 1, true).v2).toBeCloseTo(1.8, 9);
  });
});

describe('rule (c) — double-bonded central bond → V2', () => {
  it('the alkene C=C (2-2, order 2): full double bond → π = 1.0 → V2 = 12.0', () => {
    expect(tor(2, 2, 2).v2).toBeCloseTo(12.0, 9);
  });

  it('the guanidinium-type C=N (2-10, order 2): mltb 2/1 → π = 0.4 → V2 = 4.8', () => {
    // The paper's own 0.4 example: the formal double bond with
    // mltb(J) = 2 but mltb(K) = 1.
    expect(tor(2, 10, 2).v2).toBeCloseTo(4.8, 9);
  });
});

describe('the order-1 non-aromatic fall-through (the rule-(c) gate)', () => {
  it('the butadiene central single (2-2, order 1) is NOT V2 = 12 — rule (g) case (5) gives π = 0.15 → V2 = 1.8', () => {
    // The paper's cross-reference: the central C–C of butadiene gets
    // π = 0.15 from rule (g) — reachable only because rule (c) gates
    // on the j-k bond order.
    expect(tor(2, 2, 1).v2).toBeCloseTo(1.8, 9);
  });

  it('the alkane C–C (1-1, order 1) is NOT V2 = 4.8 — rule (d) gives V3 = 0.2356', () => {
    const r = tor(1, 1, 1);
    expect(r.v2).toBe(0);
    expect(r.v3).toBeCloseTo(2.12 / 9, 9);
  });
});

describe('rule (d) — both sp3 → V3 = √(V_i·V_j)/9', () => {
  it('the alkane C–C: V3 = 2.12/9 = 0.23556', () => {
    expect(tor(1, 1, 1).v3).toBeCloseTo(0.2355556, 6);
  });
});

describe('rules (e)/(f) — sp3-sp2 rotors are skipped', () => {
  it('alkyl-vinyl (1-2): the vinylic mltb ≠ 0 → zero torsion', () => {
    expect(tor(1, 2, 1).skip).toBe(true);
  });
});

describe('rule (g) — resonant single bonds → V2', () => {
  it('case (1): both atoms carry pi lone pairs (10-10) → zero torsion', () => {
    expect(tor(10, 10, 1).skip).toBe(true);
  });

  it('case (2): j = amide N (10, pilp 1 AND mltb 1), k = vinylic C (2) → π = 0.5 → V2 = 6.0', () => {
    // The paper: the 0.5 reads the PILP atom's own mltb (j here) —
    // the amide-N example; OpenBabel reads k's mltb instead.
    expect(tor(10, 2, 1).v2).toBeCloseTo(6.0, 9);
  });

  it('case (2) fallback: j = amine N (8, pilp 1, mltb 0) → π = 0.3 (carbon-row pair) → V2 = 3.6', () => {
    expect(tor(8, 2, 1).v2).toBeCloseTo(3.6, 9);
  });

  it('case (3) mirror: j = vinylic C (2), k = amide N (10) → π = 0.5 → V2 = 6.0', () => {
    expect(tor(2, 10, 1).v2).toBeCloseTo(6.0, 9);
  });

  it('case (4): mltb-1 with a non-carbon (81-2) → π = 0.4 → V2 = 4.8', () => {
    expect(tor(81, 2, 1).v2).toBeCloseTo(4.8, 9);
  });

  it('case (5): mltb 2/2 without lone pairs (2-2) → π = 0.15 → V2 = 1.8', () => {
    expect(tor(2, 2, 1).v2).toBeCloseTo(1.8, 9);
  });
});

describe('rule (h) — the saturated remainder', () => {
  it('the O–O pair: negative V2 = −√(2·2) = −2.0', () => {
    expect(tor(6, 6, 1).v2).toBeCloseTo(-2.0, 9);
  });

  it('the S–S pair: negative V2 = −√(8·8) = −8.0', () => {
    expect(tor(15, 15, 1).v2).toBeCloseTo(-8.0, 9);
  });

  it('two tricoordinate amines (8-8): V3 = √(1.5·1.5)/N_bc with N_bc = (3−1)(3−1) = 4 → 0.375', () => {
    // The eq. (22) N_bc definition — OpenBabel's crd·crd = 9 would
    // give 0.1667 here.
    expect(tor(8, 8, 1).v3).toBeCloseTo(0.375, 9);
  });

  it('mixed sp3/sp3-N (1-8): V3 = √(2.12·1.5)/((4−1)(3−1)) = 0.29721', () => {
    expect(tor(1, 8, 1).v3).toBeCloseTo(Math.sqrt(2.12 * 1.5) / 6, 9);
  });
});
