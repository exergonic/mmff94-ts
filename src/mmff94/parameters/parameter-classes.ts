/**
 * The MMFF94 parameter-class system: BTij, ATijk, TTijkl, STijk
 * (parts IV and V) and the class-scoped parameter resolution built
 * on them.
 *
 * Every MMFF94 parameter table carries a class column (first field):
 * class 0 is the general entry, and higher classes hold
 * context-specific alternatives — conjugated single bonds (bond
 * class 1), angles with BT-flagged bonds (1/2), 3-ring angles (3/5/6),
 * 4-ring angles (4/7/8), and the torsion classes below. Which class
 * applies is a CHEMICAL question, decided from the bond-type flag
 * BTij (part V, p. 620):
 *
 *   BTij = 1 when:
 *     a) the bond is single and both atom types carry the sbmb flag
 *        (mmffprop.par) — e.g. the central single bond of a
 *        conjugated diene, or a C(=O)-C(aromatic) bond. The
 *        "not aromatic" clause in the part-V comment is NOT enforced
 *        by the code; aromaticity enters only through the ring-bond
 *        check below.
 *     b) both atoms are aromatic types but the bond itself is not
 *        aromatic (e.g. the connecting bond of biphenyl).
 *   BTij = 0 otherwise — including double/triple bonds and bonds
 *   INSIDE an aromatic ring (BatchMin flags ring bonds aromatic, so
 *   they read 0 even in Kekulé input files; the in_ring check
 *   reproduces that).
 *
 * The class then selects the lookup: ATijk = BT(i,j) + BT(j,k) with
 * ring overrides (classes 3/5/6 and 4/7/8), TTijkl for torsions
 * (1 = central BT, 2 = terminal BT, 4 = 4-ring, 5 = non-aromatic
 * 5-ring with an sp3 C), and STijk is a REMAP of ATijk for
 * stretch-bend. Each class is looked up separately — the class-0
 * entry is never a stand-in for a ring or BT-flagged angle; a miss
 * falls to the empirical rules.
 *
 * The step-down chain (part I, p. 513): exact terminal types first,
 * then the EqLvl3/4/5 equivalence levels of the terminals (from
 * mmffdef.par; the central type is never reduced; level 5 is 0 =
 * wildcard for most types).
 */

import type { TypedMolecule } from '../../types';
import {
  ATOM_TYPE_PROPERTIES,
  ANGLE_PARAMS,
  BOND_PARAMS,
  TORSION_PARAMS,
  lookup_param,
  type AngleParams,
} from './index';

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

/** Bond order of (i, j) — 0 when the pair is not bonded. */
export function get_bond_order(ctx: ClassContext, i: number, j: number): number {
  const { mol } = ctx;
  const bond = mol.bonds.find(
    b => (b.atom1 === i && b.atom2 === j) || (b.atom1 === j && b.atom2 === i),
  );
  return bond ? bond.bond_order : 0;
}

/** Is bond (i, j) part of any ring? BFS for an alternate path, capped
 *  at ring size 8 (MMFF ring classes only reach 5-membered rings and
 *  the validation suite's largest rings are 8-membered; the cap just
 *  bounds the search). Only consulted for aromatic-type pairs, so the
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
  const pi = ATOM_TYPE_PROPERTIES[mol.atom_types[i]];
  const pj = ATOM_TYPE_PROPERTIES[mol.atom_types[j]];
  if (!pi || !pj || !pi.arom || !pj.arom) return false;
  return in_ring(ctx, i, j);
}

/** Element → periodic-table row (empirical rules): 0 = Z ≤ 2,
 *  1 = 3–10, 2 = 11–18, 3 = 19–36, 4 = 37–54. */
export const ELEMENT_ROW: Record<string, number> = {
  H: 0, B: 1, C: 1, N: 1, O: 1, F: 1,
  Si: 2, P: 2, S: 2, Cl: 2, Br: 3, I: 4,
};

/** GetBondType — the BTij flag for bond (i, j), see the header. */
export function bond_type_flag(ctx: ClassContext, i: number, j: number): number {
  const { mol } = ctx;
  if (get_bond_order(ctx, i, j) !== 1) return 0;

  const pi = ATOM_TYPE_PROPERTIES[mol.atom_types[i]];
  const pj = ATOM_TYPE_PROPERTIES[mol.atom_types[j]];
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

/** GetAngleType — the ATijk angle class for i-j-k: the sum of the two
 *  BT values, with ring overrides for 3- and 4-membered rings. */
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

/** GetTorsionType — the TTijkl torsion class for i-j-k-l:
 *   1 = central bond BT-flagged (e.g. the central single bond of a
 *       conjugated diene — V2 ≈ 1.8, not the alkene's 12)
 *   2 = central bond single but a terminal bond BT-flagged
 *   4 = all four atoms in the same 4-membered ring
 *   5 = all four in a non-aromatic 5-ring with at least one sp3 C
 *   0 = everything else (alkanes, alkenes, aromatic rings) */
export function torsion_class(ctx: ClassContext, i: number, j: number, k: number, l: number): number {
  const { mol } = ctx;
  const bt_ab = bond_type_flag(ctx, i, j);
  const bt_bc = bond_type_flag(ctx, j, k);
  const bt_cd = bond_type_flag(ctx, k, l);

  if (bt_bc === 1) return 1;

  // 4-ring: the closing bond i-l makes i-j-k-l a 4-cycle.
  if (ctx.adj[i].includes(l)) return 4;

  if (get_bond_order(ctx, j, k) === 1 && !is_aromatic_bond(ctx, j, k)) {
    if (bt_ab || bt_cd) return 2;
  }

  // 5-ring: a common neighbor x of i and l closes the cycle i-j-k-l-x.
  // Only non-aromatic rings qualify, and only with an sp3 carbon. The
  // aromaticity test is on the RING, not the torsion atoms: a fused
  // non-aromatic 5-ring can share aromatic-typed carbons with a benzo
  // ring (FILNOD's thiazolidine ring), and BatchMin's GetTorsionType
  // skips rings by ring->IsAromatic(), not by the atoms' own flags.
  const types = [i, j, k, l].map(a => mol.atom_types[a]);
  if (types.includes(1)) {
    for (const x of ctx.adj[i]) {
      if (x === j || x === k || !ctx.adj[l].includes(x)) continue;
      // The ring i-j-k-l-x is aromatic iff every one of its bonds is
      // aromatic; try each closing atom x — a fused system can offer
      // both an aromatic and a non-aromatic 5-ring through i and l.
      const ring_aromatic =
        is_aromatic_bond(ctx, i, j) &&
        is_aromatic_bond(ctx, j, k) &&
        is_aromatic_bond(ctx, k, l) &&
        is_aromatic_bond(ctx, l, x) &&
        is_aromatic_bond(ctx, x, i);
      if (!ring_aromatic) return 5;
    }
  }

  return 0;
}

/** GetStrBndType — the stretch-bend class, a REMAP of the angle class
 *  (part IV, p. 609): class 1/2 split by which side carries the BT flag
 *  (and the i-vs-k type order), 2→3, 3→5, 4→4, 5→6/7, 6→8, 7→9/10,
 *  8→11. The strbnd par only carries classes 0/1/2/4/5; the remapped
 *  classes miss and fall to the default stretch-bend table. */
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

// ── Part-II empirical angle fallbacks ─────────────────────────────────
// When the class-scoped lookup misses (or its entry has k_a = 0), the
// reference angle comes from the atom-type coordination numbers and the
// force constant from the part-II empirical formula (with the part-V
// element constants Z_i and C_i; C_i has no H entry, so an H-centered
// empirical angle gets k_a = 0 — the spec).

// MMFF part V, table VI: Z_i.
const ELEMENT_Z: Record<string, number> = {
  H: 1.395, C: 2.494, N: 2.711, O: 3.045, F: 2.847,
  Si: 2.350, P: 2.350, S: 2.980, Cl: 2.909, Br: 3.017, I: 3.086,
};

// MMFF part V, table X: C_i (the central atom's constant).
const ELEMENT_C: Record<string, number> = {
  B: 0.704, C: 1.016, N: 1.113, O: 1.337,
  Si: 0.811, P: 1.068, S: 1.249, Cl: 1.078, As: 0.825,
};

const ELEMENT_ATOMIC_NUMBER: Record<string, number> = {
  H: 1, B: 5, C: 6, N: 7, O: 8, F: 9, Si: 14, P: 15, S: 16, Cl: 17, Br: 35, I: 53,
};

function empirical_theta0(ctx: ClassContext, j: number, cls: number): number {
  const { mol } = ctx;
  const prop = ATOM_TYPE_PROPERTIES[mol.atom_types[j]];
  let theta0 = 120.0;
  if (prop) {
    if (prop.crd === 4) theta0 = 109.45;
    if (prop.crd === 2 && mol.atoms[j].element === 'O') theta0 = 105.0;
    // The rule is "atomic number > 10" (part II); the handful of
    // MMFF-typed elements is mapped explicitly.
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

  // Ring strain scales the constant down (part II): ×0.85 for
  // 4-rings, ×0.05 for 3-rings.
  let beta = 1.75;
  if (cls === 4 || cls === 7 || cls === 8) beta *= 0.85;
  if (cls === 3 || cls === 5 || cls === 6) beta *= 0.05;

  const r0ab = bond_parameters(ctx, i, j)?.r0 ?? 1.5;
  const r0bc = bond_parameters(ctx, j, k)?.r0 ?? 1.5;
  const rr = r0ab + r0bc;
  const D = (r0ab - r0bc) / (rr * rr);
  const rad2 = (Math.PI / 180) * (Math.PI / 180);
  return (beta * Za * Cb * Zc * Math.exp(-2 * D)) / (rr * theta0 * theta0 * rad2);
}

// ── Class-scoped lookups ──────────────────────────────────────────────

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

/** Step-down angle lookup (part I, p. 513): exact terminal types (the
 *  table stores terminals sorted, so the key is min/max), then the
 *  EqLvl3/4/5 equivalence levels of the terminals — the central type
 *  is never reduced; level 5 is 0 = wildcard for most types. */
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
    const pi = ATOM_TYPE_PROPERTIES[ti]?.[lvl] ?? ti;
    const pk = ATOM_TYPE_PROPERTIES[tk]?.[lvl] ?? tk;
    const p = table[key(pi, pk)];
    if (p) return p;
  }
  return undefined;
}

/** Full angle-parameter resolution: class-scoped step-down lookup,
 *  falling to the empirical θ₀ rules and the empirical force constant
 *  (part II) when the lookup misses or returns k_a = 0. */
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

  const linear = ATOM_TYPE_PROPERTIES[tj]?.lin === 1;
  return { k_a, theta0, linear };
}

/** Torsion step-down chain (part I, p. 513): exact, then the
 *  asymmetric (EqLvl3, EqLvl5) terminal reductions. The par file
 *  stores each entry in ONE canonical direction; the order index
 *  decides which direction that is, and only that direction's chain
 *  is consulted — the other direction would hit wildcard defaults
 *  (e.g. '0-0-1-1-0') before the exact reversed entry. */
export function lookup_torsion(
  cls: number,
  ti: number,
  tj: number,
  tk: number,
  tl: number,
): { v1: number; v2: number; v3: number } | undefined {
  const lvl3 = (t: number) => ATOM_TYPE_PROPERTIES[t]?.lvl3 ?? t;
  const lvl5 = (t: number) => ATOM_TYPE_PROPERTIES[t]?.lvl5 ?? t;
  const chain = (a: number, b: number, c: number, d: number) => {
    const keys = [
      `${cls}-${a}-${b}-${c}-${d}`,
      `${cls}-${lvl3(a)}-${b}-${c}-${lvl5(d)}`,
      `${cls}-${lvl5(a)}-${b}-${c}-${lvl3(d)}`,
      `${cls}-${lvl5(a)}-${b}-${c}-${lvl5(d)}`,
    ];
    for (const key of keys) {
      const p = TORSION_PARAMS[key];
      if (p) return p;
    }
    return undefined;
  };

  // The order index is the canonical-direction discriminator (M = 136,
  // the type-code base in the par's key encoding).
  const M = 136;
  const order =
    tk * M ** 3 + tj * M ** 2 + tl * M + ti - (tj * M ** 3 + tk * M ** 2 + ti * M + tl);
  const entry = order >= 0 ? chain(ti, tj, tk, tl) : chain(tl, tk, tj, ti);
  if (!entry) return undefined;
  let v1 = 0;
  let v2 = 0;
  let v3 = 0;
  for (const term of entry.terms) {
    if (term.periodicity === 1) v1 = term.V;
    else if (term.periodicity === 2) v2 = term.V;
    else v3 = term.V;
  }
  return { v1, v2, v3 };
}
