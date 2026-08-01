/**
 * The BTij bond-type flag and the parameter classes it selects
 * (MMFF part V, p. 620).
 *
 * BTij = 1 when:
 *   a) the bond is single and both atoms are non-aromatic types with
 *      the sbmb flag set in mmffprop.par (e.g. the central single
 *      bond of a conjugated diene), or
 *   b) both atoms are aromatic types but the bond itself is not
 *      aromatic (e.g. the connecting bond of biphenyl).
 * BTij = 0 otherwise — including double/triple bonds and bonds INSIDE
 * an aromatic ring (BatchMin flags ring bonds aromatic, so they read 0
 * even in Kekulé input files).
 *
 * BTij selects the bond-stretch class (class 1 = conjugated single
 * bonds) and feeds the angle class ATijk: classes 1/2 for one/two
 * BT-flagged bonds, 3/5/6 for 3-rings and 4/7/8 for 4-rings by the sum
 * of the two BT values, else 0. Each class is looked up separately —
 * the class-0 entry is never a stand-in for a ring or BT-flagged angle
 * (a miss falls to the empirical rules below, Halgren part II).
 */

import type { TypedMolecule } from '../../types';
import {
  ATOM_PROPERTIES,
  ANGLE_PARAMS,
  BOND_PARAMS,
  lookup_param,
  type AngleParams,
} from '../parameters';

/** Per-term call context: molecule, adjacency, and the memoized
 *  in-ring checks shared by every class query in one energy pass. */
export interface ClassContext {
  mol: TypedMolecule;
  adj: number[][];
  ringMemo: Map<string, boolean>;
}

export function make_class_context(mol: TypedMolecule, adj: number[][]): ClassContext {
  return { mol, adj, ringMemo: new Map() };
}

/** Is bond (i, j) part of any ring? BFS for an alternate path, capped
 *  at ring size 8. Only consulted for aromatic-type pairs, so the
 *  answer distinguishes an aromatic ring bond (BTij = 0) from an
 *  external aromatic-aromatic bond like biphenyl's (BTij = 1). */
function in_ring(ctx: ClassContext, i: number, j: number): boolean {
  const key = i < j ? `${i}-${j}` : `${j}-${i}`;
  const memoized = ctx.ringMemo.get(key);
  if (memoized !== undefined) return memoized;

  const { adj } = ctx;
  let found = false;
  const queue: [number, number][] = [];
  const seen = new Set<number>([i]);
  // Start from i's neighbors, skipping the direct edge to j.
  for (const nb of adj[i]) {
    if (nb !== j && !seen.has(nb)) {
      seen.add(nb);
      queue.push([nb, 1]);
    }
  }
  while (queue.length > 0 && !found) {
    const [node, depth] = queue.shift()!;
    if (depth >= 8) continue;
    for (const nb of adj[node]) {
      if (nb === j) {
        found = true;
        break;
      }
      if (!seen.has(nb)) {
        seen.add(nb);
        queue.push([nb, depth + 1]);
      }
    }
  }
  ctx.ringMemo.set(key, found);
  return found;
}

/** Is bond (i, j) part of an aromatic ring? Both endpoints are
 *  aromatic-typed and the bond lies in a ring (Kekulé input files
 *  carry no aromatic flags; BatchMin's perception marks ring bonds). */
export function is_aromatic_bond(ctx: ClassContext, i: number, j: number): boolean {
  const { mol } = ctx;
  const pi = ATOM_PROPERTIES[mol.atom_types[i]];
  const pj = ATOM_PROPERTIES[mol.atom_types[j]];
  if (!pi || !pj || !pi.arom || !pj.arom) return false;
  return in_ring(ctx, i, j);
}

/** Element → periodic-table row (empirical rules): 0 = Z ≤ 2,
 *  1 = 3–10, 2 = 11–18, 3 = 19–36, 4 = 37–54. */
export const ELEMENT_ROW: Record<string, number> = {
  H: 0, B: 1, C: 1, N: 1, O: 1, F: 1,
  Si: 2, P: 2, S: 2, Cl: 2, Br: 3, I: 4,
};

/** GetBondType — the BTij flag for bond (i, j). */
export function bond_type_flag(ctx: ClassContext, i: number, j: number): number {
  const { mol } = ctx;
  const bond = mol.bonds.find(
    b => (b.atom1 === i && b.atom2 === j) || (b.atom1 === j && b.atom2 === i),
  );
  if (!bond || bond.bond_order !== 1) return 0;

  const pi = ATOM_PROPERTIES[mol.atom_types[i]];
  const pj = ATOM_PROPERTIES[mol.atom_types[j]];
  if (!pi || !pj) return 0;

  if (pi.arom && pj.arom) {
    // Aromatic pair: 0 inside an aromatic ring (BatchMin's ring bonds
    // are aromatic), 1 for a non-aromatic bond like biphenyl's.
    return in_ring(ctx, i, j) ? 0 : 1;
  }
  // Case (a): both types have the sbmb flag — the code returns 1 for
  // ANY sbmb pair (the "not aromatic" clause in the part-V comment is
  // not enforced; aromaticity only enters through the ring-bond check
  // above). This is what makes C(=O)-C(ar) and C=C-C(ar) angles
  // class 1/2 and conjugated single bonds class 1.
  if (pi.sbmb && pj.sbmb) return 1;
  return 0;
}

/** GetAngleType — the ATijk angle class for i-j-k. */
export function angle_class(ctx: ClassContext, i: number, j: number, k: number): number {
  const sum = bond_type_flag(ctx, i, j) + bond_type_flag(ctx, j, k);

  // 3-ring: i and k are directly bonded (the triangle i-j-k).
  if (ctx.adj[i].includes(k)) {
    switch (sum) {
      case 0: return 3;
      case 1: return 5;
      default: return 6;
    }
  }
  // 4-ring: i and k share a common neighbor other than j (cycle i-j-k-x).
  const common = ctx.adj[i].find(x => x !== j && ctx.adj[k].includes(x));
  if (common !== undefined) {
    switch (sum) {
      case 0: return 4;
      case 1: return 7;
      default: return 8;
    }
  }
  return sum;
}

// ── Empirical fallbacks (Halgren part II) ─────────────────────────────

// MMFF part V, table VI: Z_i (used in the empirical angle constant).
const ELEMENT_Z: Record<string, number> = {
  H: 1.395, C: 2.494, N: 2.711, O: 3.045, F: 2.847,
  Si: 2.350, P: 2.350, S: 2.980, Cl: 2.909, Br: 3.017, I: 3.086,
};

// MMFF part V, table X: C_i (the central atom's constant; H has none,
// so an H-centered empirical angle gets k_a = 0).
const ELEMENT_C: Record<string, number> = {
  B: 0.704, C: 1.016, N: 1.113, O: 1.337,
  Si: 0.811, P: 1.068, S: 1.249, Cl: 1.078, As: 0.825,
};

const ELEMENT_ATOMIC_NUMBER: Record<string, number> = {
  H: 1, B: 5, C: 6, N: 7, O: 8, F: 9, Si: 14, P: 15, S: 16, Cl: 17, Br: 35, I: 53,
};

function empirical_theta0(ctx: ClassContext, j: number, cls: number): number {
  const { mol } = ctx;
  const prop = ATOM_PROPERTIES[mol.atom_types[j]];
  let theta0 = 120.0;
  if (prop) {
    if (prop.crd === 4) theta0 = 109.45;
    if (prop.crd === 2 && mol.atoms[j].element === 'O') theta0 = 105.0;
    if ((ELEMENT_ATOMIC_NUMBER[mol.atoms[j].element] ?? 0) > 10) theta0 = 95.0;
    if (prop.lin) theta0 = 180.0;
    if (prop.crd === 3 && prop.val === 3 && !prop.mltb) {
      theta0 = mol.atoms[j].element === 'N' ? 107.0 : 92.0;
    }
  }
  // Small-ring angles are forced to the ring geometry.
  if (cls === 3 || cls === 5 || cls === 6) theta0 = 60.0;
  if (cls === 4 || cls === 7 || cls === 8) theta0 = 90.0;
  return theta0;
}

function empirical_ka(ctx: ClassContext, i: number, j: number, k: number, theta0: number, cls: number): number {
  const { mol } = ctx;
  const Za = ELEMENT_Z[mol.atoms[i].element] ?? 0;
  const Cb = ELEMENT_C[mol.atoms[j].element] ?? 0; // 0 for H — spec
  const Zc = ELEMENT_Z[mol.atoms[k].element] ?? 0;

  let beta = 1.75;
  if (cls === 4 || cls === 7 || cls === 8) beta *= 0.85; // 4-ring
  if (cls === 3 || cls === 5 || cls === 6) beta *= 0.05; // 3-ring

  const r0ab = bond_parameters(ctx, i, j)?.r0 ?? 1.5;
  const r0bc = bond_parameters(ctx, j, k)?.r0 ?? 1.5;
  const rr = r0ab + r0bc;
  const D = (r0ab - r0bc) / (rr * rr);
  const rad2 = (Math.PI / 180) * (Math.PI / 180);
  return (beta * Za * Cb * Zc * Math.exp(-2 * D)) / (rr * theta0 * theta0 * rad2);
}

// ── Public lookups ────────────────────────────────────────────────────

/** Bond stretch parameters, class-aware: a BTij=1 bond (conjugated
 *  single bond) uses the class-1 entry when one exists — the class-0
 *  entry for the same pair is often the DOUBLE-bond parameter (e.g.
 *  '0-2-2' is C=C while '1-2-2' is the conjugated single bond). */
export function bond_parameters(
  ctx: ClassContext,
  i: number,
  j: number,
): { k_b: number; r0: number } | undefined {
  const { mol } = ctx;
  const t_min = Math.min(mol.atom_types[i], mol.atom_types[j]);
  const t_max = Math.max(mol.atom_types[i], mol.atom_types[j]);
  if (bond_type_flag(ctx, i, j) === 1) {
    const class1 = BOND_PARAMS[`1-${t_min}-${t_max}`];
    if (class1) return class1;
  }
  return lookup_param(BOND_PARAMS, [t_min, t_max]);
}

/** Step-down angle lookup (MMFF part I, p. 513): exact terminal types,
 *  then the EqLvl3/4/5 equivalence levels of the terminals (the central
 *  type is never reduced; level 5 is 0 = wildcard for most types). */
function lookup_angle_entry(
  table: Record<string, AngleParams>,
  cls: number,
  ti: number,
  tj: number,
  tk: number,
): AngleParams | undefined {
  const key = (a: number, b: number) => `${cls}-${Math.min(a, b)}-${tj}-${Math.max(a, b)}`;
  const exact = table[key(ti, tk)];
  if (exact) return exact;
  for (const lvl of ['lvl3', 'lvl4', 'lvl5'] as const) {
    const pi = ATOM_PROPERTIES[ti]?.[lvl] ?? ti;
    const pk = ATOM_PROPERTIES[tk]?.[lvl] ?? tk;
    const p = table[key(pi, pk)];
    if (p) return p;
  }
  return undefined;
}

/** GetStrBndType — the stretch-bend class, a REMAP of the angle class
 *  (part IV, p. 609): class 1/2 split by which side carries the BT flag
 *  (and the i-vs-k type order), 2→3, 3→5, 4→4, 5→6/7, 6→8, 7→9/10,
 *  8→11. The strbnd par only carries classes 0/1/2/4/5; the remapped
 *  classes miss and fall to the default-fsb table. */
export function strbnd_type(ctx: ClassContext, i: number, j: number, k: number): number {
  const { mol } = ctx;
  const atabc = angle_class(ctx, i, j, k);
  const btab = bond_type_flag(ctx, i, j);
  const btbc = bond_type_flag(ctx, j, k);
  const inverse = mol.atom_types[i] > mol.atom_types[k];
  switch (atabc) {
    case 0: return 0;
    case 1:
      if (btab) return inverse ? 2 : 1;
      if (btbc) return inverse ? 1 : 2;
      return 0; // unreachable: atabc 1 implies one BT flag
    case 2: return 3;
    case 3: return 5;
    case 4: return 4;
    case 5:
      if (btab) return inverse ? 7 : 6;
      return inverse ? 6 : 7;
    case 6: return 8;
    case 7:
      if (btab) return inverse ? 10 : 9;
      return inverse ? 9 : 10;
    case 8: return 11;
    default: return 0;
  }
}

/**
 * Full angle-parameter resolution: class-scoped step-down lookup (exact
 * → EqLvl3 → EqLvl4 → EqLvl5), falling to the empirical θ₀ rules and the
 * empirical force constant (Halgren part II) when the lookup misses or
 * returns k_a = 0.
 */
export function angle_parameters(
  ctx: ClassContext,
  i: number,
  j: number,
  k: number,
): { k_a: number; theta0: number; linear: boolean } {
  const { mol } = ctx;
  const ti = mol.atom_types[i];
  const tj = mol.atom_types[j];
  const tk = mol.atom_types[k];
  const cls = angle_class(ctx, i, j, k);

  // The class is part of the key; class 0 is never consulted for a
  // ring or BT-flagged angle (a miss means empirical rules).
  const params = lookup_angle_entry(ANGLE_PARAMS, cls, ti, tj, tk);

  let k_a: number;
  let theta0: number;
  if (params) {
    k_a = params.k_a;
    theta0 = params.theta0;
  } else {
    k_a = 0;
    theta0 = empirical_theta0(ctx, j, cls);
  }
  if (k_a === 0) k_a = empirical_ka(ctx, i, j, k, theta0, cls);

  const linear = ATOM_PROPERTIES[tj]?.lin === 1;
  return { k_a, theta0, linear };
}
