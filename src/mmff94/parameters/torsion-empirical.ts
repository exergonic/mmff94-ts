/**
 * The part-IV empirical torsion rules (MMFF part IV, p. 631-632),
 * used when the TTijkl step-down chain misses entirely.
 *
 * Rules (b)-(h) always set something; a few cases SKIP the torsion
 * (linear centers, or specific crd/val/mltb combinations in rules
 * (e)/(f)/(g)). The constants U_i, V_i and W_i are element-based
 * (part V, table X and the W values given in part IV).
 */

import { ATOM_TYPE_PROPERTIES } from './index';
import {
  is_aromatic_bond,
  get_bond_order,
  ELEMENT_ROW,
  type ClassContext,
} from './parameter-classes';

// MMFF part V, table X: U_i.
const ELEMENT_U: Record<string, number> = {
  C: 2.0, N: 2.0, O: 2.0, Si: 1.25, P: 1.25, S: 1.25,
};
// MMFF part V, table X: V_i.
const ELEMENT_V: Record<string, number> = {
  C: 2.12, N: 1.5, O: 0.2, Si: 1.22, P: 2.4, S: 0.49,
};

/** The empirical torsion barrier heights for i-j-k-l, or skip = true
 *  when the rules say this torsion does not exist (linear centers and
 *  the rule (e)/(f)/(g) exclusions). */
export function empirical_torsion(
  ctx: ClassContext,
  i: number,
  j: number,
  k: number,
): { v1: number; v2: number; v3: number; skip: boolean } {
  const { mol } = ctx;
  const tj = mol.atom_types[j];
  const tk = mol.atom_types[k];
  const pj = ATOM_TYPE_PROPERTIES[tj];
  const pk = ATOM_TYPE_PROPERTIES[tk];
  const ej = mol.atoms[j].element;
  const ek = mol.atoms[k].element;
  const ub = ELEMENT_U[ej] ?? 0;
  const uc = ELEMENT_U[ek] ?? 0;
  const vb = ELEMENT_V[ej] ?? 0;
  const vc = ELEMENT_V[ek] ?? 0;
  const v = { v1: 0, v2: 0, v3: 0, skip: false };
  let found = false;

  // rule (a): linear centers carry no torsion
  if (pj?.lin || pk?.lin) { v.skip = true; return v; }

  const aromatic_central = is_aromatic_bond(ctx, j, k);

  // rules (b)/(c): V2 from the U parameters
  if (aromatic_central) {
    // rule (b): aromatic central bond — π = 0.5 without lone pairs
    // on either atom, 0.3 with; β = 3 for the val 3/4 combination,
    // 6 otherwise
    const pi_bc = (!pj?.pilp && !pk?.pilp) ? 0.5 : 0.3;
    const beta =
      (pj?.val === 3 && pk?.val === 4) || (pj?.val === 4 && pk?.val === 3) ? 3.0 : 6.0;
    v.v2 = beta * pi_bc * Math.sqrt(ub * uc);
    found = true;
  } else {
    // rule (c): π = 1.0 only for a double-bonded sp2 pair whose own
    // i-j bond is the double bond, else 0.4
    const pi_bc =
      pj?.mltb === 2 && pk?.mltb === 2 &&
      get_bond_order(ctx, i, j) === 2 && !is_aromatic_bond(ctx, i, j)
        ? 1.0
        : 0.4;
    v.v2 = 6.0 * pi_bc * Math.sqrt(ub * uc);
    found = true;
  }

  // rule (d): both sp3 → V3
  if (!found && pj?.crd === 4 && pk?.crd === 4) {
    v.v3 = Math.sqrt(vb * vc) / 9.0;
    found = true;
  }

  // rules (e)/(f): sp3/sp2 mixed — no torsion for the unsaturated
  // sp2 combinations listed (a saturated sp3-sp2 bond falls through
  // to the remaining rules)
  const unsaturated_sp2 = (t: number): boolean => {
    const p = ATOM_TYPE_PROPERTIES[t];
    if (!p) return false;
    if (p.crd === 3) return p.val === 4 || p.val === 34 || p.mltb !== 0;
    if (p.crd === 2) return p.val === 3 || p.mltb !== 0;
    return false;
  };
  if (!found && pj?.crd === 4 && pk?.crd !== 4) {
    if (unsaturated_sp2(tk)) { v.skip = true; return v; }
  } else if (!found && pk?.crd === 4 && pj?.crd !== 4) {
    if (unsaturated_sp2(tj)) { v.skip = true; return v; }
  }

  // rule (g): order-1 central bond between mltb/pilp-carrying types
  // → V2 (the π value depends on which side carries what)
  if (!found) {
    const central_single =
      get_bond_order(ctx, j, k) === 1 && !aromatic_central;
    if (
      central_single &&
      ((pj?.mltb && pk?.mltb) || (pj?.mltb && pk?.pilp) || (pk?.mltb && pj?.pilp))
    ) {
      if (pj?.pilp && pk?.pilp) { v.skip = true; return v; } // case (1)
      let pi_bc = 0.15;
      if (pj?.pilp && pk?.mltb) { // case (2)
        if (pk.mltb === 1) pi_bc = 0.5;
        else if (ELEMENT_ROW[ej] === 1 && ELEMENT_ROW[ek] === 1) pi_bc = 0.3;
        else pi_bc = 0.15;
        found = true;
      }
      if (pk?.pilp && pj?.mltb) { // case (3)
        if (pj.mltb === 1) pi_bc = 0.5;
        else if (ELEMENT_ROW[ej] === 1 && ELEMENT_ROW[ek] === 1) pi_bc = 0.3;
        else pi_bc = 0.15;
        found = true;
      }
      if (!found && (pj?.mltb === 1 || pk?.mltb === 1) && (ej !== 'C' || ek !== 'C')) {
        pi_bc = 0.4;
        found = true;
      }
      if (!found) pi_bc = 0.15;
      v.v2 = 6.0 * pi_bc * Math.sqrt(ub * uc);
      found = true;
    }
  }

  // rule (h): O/S central pair → negative V2 (W = 2 for O, 8 for S);
  // otherwise V3 from the V parameters
  if (!found) {
    const o_s = (e: string) => e === 'O' || e === 'S';
    if (o_s(ej) && o_s(ek)) {
      v.v2 = -Math.sqrt((ej === 'O' ? 2.0 : 8.0) * (ek === 'O' ? 2.0 : 8.0));
    } else {
      v.v3 = Math.sqrt(vb * vc) / ((pj?.crd ?? 1) * (pk?.crd ?? 1));
    }
  }
  return v;
}
