// The part V empirical torsion rules (empirical.ts), pinned against the
// paper (part V pp. 631-632) AND the corroborated implementations. The
// suite never exercises the rules (the par's wildcards always resolve,
// validated 753/753), so these tests pin the protocol against
// hand-computed values.
//
// The 2026-08-06 vinyl-phosphine arbitration: the paper's rule (c)
// text gates on the j-k formal bond order of 2, but BOTH reference
// implementations (OpenBabel, Tinker) treat it as the else of the
// aromatic rule — every non-aromatic central bond gets eq. (21) with
// π = 1.0 when both atoms carry mltb 2, else π = 0.4. Measured on
// vinyl phosphine's C–P dihedrals: V2 = 6·0.4·√(U_C·U_P) = 3.795 in
// both references (the paper's rule (g) case (5) would give 1.423).
// The corroboration wins; the paper's rules (d)-(h) are unreachable
// for non-aromatic bonds and no longer pinned.
//
// Table X constants used below: U_C = U_N = U_O = 2.0, U_S = 1.25 (so
// √(U·U) = 2 for the first-row pairs, 1.25 for S–S).
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

describe('rule (c) — the corroborated universal non-aromatic → V2', () => {
  it('the alkene C=C (2-2, order 2): full double bond → π = 1.0 → V2 = 12.0', () => {
    expect(tor(2, 2, 2).v2).toBeCloseTo(12.0, 9);
  });

  it('the guanidinium-type C=N (2-10, order 2): mltb 2/1 → π = 0.4 → V2 = 4.8', () => {
    // The paper's own 0.4 example: the formal double bond with
    // mltb(J) = 2 but mltb(K) = 1.
    expect(tor(2, 10, 2).v2).toBeCloseTo(4.8, 9);
  });

  it('the butadiene central single (2-2, order 1): mltb 2/2 → π = 1.0 → V2 = 12.0', () => {
    // The corroborated reading: both references give π = 1.0 for the
    // central single of butadiene (the paper's rule (g) case (5)
    // assigns 0.15, but the implementations agree and the suite never
    // exercises the rules).
    expect(tor(2, 2, 1).v2).toBeCloseTo(12.0, 9);
  });

  it('the alkane C–C (1-1, order 1): π = 0.4 → V2 = 4.8', () => {
    const r = tor(1, 1, 1);
    expect(r.skip).toBe(false);
    expect(r.v2).toBeCloseTo(4.8, 9);
    expect(r.v3).toBe(0);
  });

  it('the amide N–vinyl (10-2): mltb 1/2 → π = 0.4 → V2 = 4.8', () => {
    // The paper's rule (g) case (2) would give π = 0.5 (mltb(J) = 1);
    // the corroborated rule (c) reads only the mltb-2 pair.
    expect(tor(10, 2, 1).v2).toBeCloseTo(4.8, 9);
  });

  it('the vinyl C–P (2-26, order 1): second-row lone pair → π = 0.4 → V2 = 3.795', () => {
    // The measured arbitration case: both references resolve the
    // vinyl phosphine P-dihedrals to 6·0.4·√(2.0·1.25) = 3.795.
    expect(tor(2, 26, 1).v2).toBeCloseTo(3.7947, 4);
  });

  it('the O–O pair (6-6): π = 0.4 → V2 = 4.8', () => {
    expect(tor(6, 6, 1).v2).toBeCloseTo(4.8, 9);
  });

  it('the S–S pair (15-15): π = 0.4, U_S = 1.25 → V2 = 3.0', () => {
    expect(tor(15, 15, 1).v2).toBeCloseTo(3.0, 9);
  });

  it('the amine N–N (8-8): π = 0.4 → V2 = 4.8', () => {
    expect(tor(8, 8, 1).v2).toBeCloseTo(4.8, 9);
  });
});
