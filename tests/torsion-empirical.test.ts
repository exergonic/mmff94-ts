// The part V empirical torsion rules (empirical.ts), pinned against the
// paper (part V pp. 631-632) AND the 761-suite's generated rows. The
// ERULE fragments (added in the suite's Nov 1998 revision) exercise the
// rules the old suite never reached: ERULE_03's P–Si resolves through
// eq. (22) (V3 = 0.285 — NOT rule (c)'s V2 = 3.0), proving rule (c) is
// gated on the formal bond order of 2 as the paper's text says, the
// (8,1)/(15,1) rows pin rules (d)-(h), and the both-pilp dihedrals of
// ERULE_01/02/04/08 (central (15,8)/(8,15)/(8,8), τ = 7.8-22.7° in the
// reference geometries) pin rule (g) case (1): NO V2 from rule (g) —
// but rule (h) still assigns the V3 = √(V_b·V_c)/N_bc (measured
// 2026-08-13: a both-pilp skip left the reference's torsion 0.29-0.41
// kcal/mol higher on those four molecules; the earlier "τ ≈ 60°, green
// both ways" claim was never checked). Tinker's ktors zeroes both-pilp
// torsions entirely — a Tinker deviation from the reference.
//
// The 2026-08-06 vinyl-phosphine arbitration is SUPERSEDED: the old
// "universal reading" (rule (c) as the else of the aromatic rule — π =
// 1.0/0.4 for every non-aromatic bond) was measured against OpenBabel
// and Tinker, which both apply eq. (21) to order-1 bonds; the suite's
// own generated rows show the reference follows the paper's order-2
// gate (the vinyl-phosphine C–P now resolves through rule (g) case (3)
// → π = 0.15 → V2 = 1.423, and OB's 3.795 is a documented OB deviation).
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

describe('rule (c) — formal bond order 2 → V2', () => {
  it('the alkene C=C (2-2, order 2): full double bond → π = 1.0 → V2 = 12.0', () => {
    expect(tor(2, 2, 2).v2).toBeCloseTo(12.0, 9);
  });

  it('the guanidinium-type C=N (2-10, order 2): mltb 2/1 → π = 0.4 → V2 = 4.8', () => {
    // The paper's own 0.4 example: the formal double bond with
    // mltb(J) = 2 but mltb(K) = 1.
    expect(tor(2, 10, 2).v2).toBeCloseTo(4.8, 9);
  });

  it('the order-1 cases fall through to rules (d)-(h), never to rule (c)', () => {
    // ERULE_03's P–Si (order 1) resolves V3 = √(2.4·1.22)/6 = 0.285 in
    // the reference — the suite's own generated row — NOT rule (c)'s
    // V2 = 6·0.4·√(U_P·U_Si) = 3.0. The rule is gated on the bond
    // order, exactly as the paper's text says.
    const r = tor(26, 19, 1);
    expect(r.v3).toBeCloseTo(0.2852, 4);
    expect(r.v2).toBe(0);
  });
});

describe('rules (d)-(h) — order-1 non-aromatic central bonds', () => {
  it('the butadiene central single (2-2, order 1): rule (g) case (5) → π = 0.15 → V2 = 1.8', () => {
    // The paper's own example for case (5): "the central C-C bond in
    // butadiene" — mltb 2/2, no lone pairs, no mltb-1 → π = 0.15.
    // (OpenBabel and Tinker apply rule (c) here and give 12.0 — a
    // documented deviation from the paper and the suite reference.)
    expect(tor(2, 2, 1).v2).toBeCloseTo(1.8, 9);
  });

  it('the alkane C–C (1-1, order 1): rule (d) both tetracoordinate → V3 = √(2.12·2.12)/9', () => {
    const r = tor(1, 1, 1);
    expect(r.skip).toBe(false);
    expect(r.v3).toBeCloseTo(2.12 / 9, 9);
    expect(r.v2).toBe(0);
  });

  it('the amide N–vinyl (10-2): rule (g) case (2) → π = 0.5 → V2 = 6.0', () => {
    // pilp(10) = 1, mltb(2) = 2, and the PILP atom's own mltb = 1
    // (amide N) → the paper's 0.5 (the old universal reading gave
    // 4.8 via case (4); the paper's case (2) fires now).
    expect(tor(10, 2, 1).v2).toBeCloseTo(6.0, 9);
  });

  it('the vinyl C–P (2-26, order 1): rule (g) case (3) → π = 0.15 → V2 = 1.423', () => {
    // pilp(26) = 1 (the tricoordinate P lone pair), mltb(2) = 2; the
    // PILP atom's mltb = 0, and P is not in carbon's row → π = 0.15.
    // (The old universal reading gave 3.795 — OpenBabel's current
    // value; the suite's order-2 gating makes the reference follow
    // the paper here.)
    expect(tor(2, 26, 1).v2).toBeCloseTo(6 * 0.15 * Math.sqrt(2.0 * 1.25), 9);
  });

  it('the O–O pair (6-6): rule (g) case (1) — both pilp, no mltb → no V2; rule (h) gives V2 = −√(2·2) = −2', () => {
    // Case (1) suppresses only rule (g)'s V2; rule (h)'s O/S special
    // still assigns the negative V2 (the paper's W = 2 for O). No
    // suite case exercises this pair — the ERULE both-pilp pairs are
    // (15,8)/(8,15)/(8,8) — the value follows the paper. (The
    // 2026-08-10 Tinker-arbitrated "no torsion" reading was
    // disproven by the suite; see the header.)
    const r = tor(6, 6, 1);
    expect(r.skip).toBe(false);
    expect(r.v2).toBeCloseTo(-2.0, 9);
    expect(r.v3).toBe(0);
  });

  it('the S–S pair (15-15): rule (h) gives V2 = −√(8·8) = −8', () => {
    // Same arbitration as (6,6): rule (h)'s O/S special applies.
    const r = tor(15, 15, 1);
    expect(r.skip).toBe(false);
    expect(r.v2).toBeCloseTo(-8.0, 9);
    expect(r.v3).toBe(0);
  });

  it('the amine N–N (8-8): rule (h) gives V3 = √(1.5·1.5)/4 = 0.375', () => {
    // Suite-pinned: ERULE_08's both-pilp (8,8) dihedral sits at
    // τ = 14.9° in the reference geometry, where the V3 term is near
    // its maximum — the BatchMin torsion total matches
    // 0.375·(1+cos3τ)/2 to 5 decimals.
    const r = tor(8, 8, 1);
    expect(r.skip).toBe(false);
    expect(r.v3).toBeCloseTo(1.5 / 4, 9);
    expect(r.v2).toBe(0);
  });

  it('the amine N–S (8-15): rule (h) gives V3 = √(1.5·0.48)/2 = 0.4243', () => {
    // Suite-pinned: ERULE_01/02/04's both-pilp (15,8)/(8,15)
    // dihedrals (τ = 7.8-22.7° in the reference geometries) match
    // √(1.5·0.48)/((3−1)·(2−1)) = 0.4243 to 5 decimals in the
    // BatchMin torsion totals.
    const r = tor(8, 15, 1);
    expect(r.skip).toBe(false);
    expect(r.v3).toBeCloseTo(Math.sqrt(1.5 * 0.48) / 2, 9);
    expect(r.v2).toBe(0);
  });
});
