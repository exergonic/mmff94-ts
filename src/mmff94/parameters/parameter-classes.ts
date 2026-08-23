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

import type { TypedMolecule } from '../../types.js';
import {
  ATOM_TYPE_PROPERTIES,
  ANGLE_PARAMS,
  BOND_PARAMS,
  TORSION_PARAMS,
  lookup_param,
  type AngleParams,
} from './index.js';
import { empirical_bond_length, empirical_ka, empirical_theta0 } from './empirical.js';
import {
  find_aromatic_rings,
  find_ring_atoms,
  type AromaticRing,
} from '../assign-atom-types.js';

/** Per-term call context: molecule, adjacency, and the aromatic-ring
 *  perception shared by every class query in one energy pass. The
 *  BTij rule (part V p. 620) distinguishes "same aromatic ring"
 *  (BTij = 0) from "different aromatic rings" (BTij = 1), so the
 *  context carries the per-atom aromatic RING sets, not just
 *  membership. */
export interface ClassContext {
  mol: TypedMolecule;
  adj: number[][];
  /** The aromatic rings containing each atom (shared ring objects —
   *  identity comparison is ring equality). */
  aromatic_rings: Map<number, AromaticRing[]>;
}

export function make_class_context(mol: TypedMolecule, adj: number[][]): ClassContext {
  // The aromatic perception walks the neighbor list with bond orders;
  // rebuild it from the molecule's bonds (the passed adjacency is
  // order-free and shared with the pair-enumeration loops).
  const n = mol.atoms.length;
  const adj_ordered: { nbr: number; order: number }[][] = Array.from({ length: n }, () => []);
  for (const bond of mol.bonds) {
    adj_ordered[bond.atom1].push({ nbr: bond.atom2, order: bond.bond_order });
    adj_ordered[bond.atom2].push({ nbr: bond.atom1, order: bond.bond_order });
  }
  const is_ring = find_ring_atoms(adj_ordered, n);
  const aromatic = find_aromatic_rings(adj_ordered, mol, is_ring);
  return { mol, adj, aromatic_rings: aromatic.rings_of };
}

// Class contexts are pure functions of TOPOLOGY (bonds → ring/aromatic
// perception), never of coordinates. Every term rebuilds one per call —
// seven energy + seven gradient perceptions per optimization oracle
// step on data that cannot change while atoms move. Cache per molecule
// object; the mutation contract is the optimizer's own (coordinates
// move; chemistry does not — see nonbonded-context.ts). Keyed by the
// adjacency IDENTITY as well: a caller may pass a freshly built adj
// for an untyped variant, and the cache entry must not outlive the
// molecule it was built from.
const class_context_cache = new WeakMap<TypedMolecule, Map<number[][], ClassContext>>();

export function class_context_for(mol: TypedMolecule, adj: number[][]): ClassContext {
  let by_adj = class_context_cache.get(mol);
  if (!by_adj) {
    by_adj = new Map();
    class_context_cache.set(mol, by_adj);
  }
  let ctx = by_adj.get(adj);
  if (!ctx) {
    ctx = make_class_context(mol, adj);
    by_adj.set(adj, ctx);
  }
  return ctx;
}

/** Bond order of (i, j) — 0 when the pair is not bonded. */
export function get_bond_order(ctx: ClassContext, i: number, j: number): number {
  const { mol } = ctx;
  const bond = mol.bonds.find(
    b => (b.atom1 === i && b.atom2 === j) || (b.atom1 === j && b.atom2 === i),
  );
  return bond ? bond.bond_order : 0;
}

/** Is bond (i, j) part of an aromatic ring? Pure ring perception —
 *  the bond's ring is in the aromatic set (Kekulé input files carry
 *  no aromatic flags; BatchMin's perception marks ring bonds). This
 *  is the MMFF meaning of an "aromatic bond": the torsion class-2
 *  branch must not fire on one (TAJSUS's triazole C(80)–N(81)
 *  resolves class 2 / V2 = 4.8 where the reference skips to class 0 /
 *  V2 = 4.0), and the class-5 ring test checks ring aromaticity by
 *  every bond of the ring (FILNOD's thiazolidine). */
export function is_aromatic_bond(ctx: ClassContext, i: number, j: number): boolean {
  return in_aromatic_ring(ctx, i, j);
}

/** Is bond (i, j) an edge of an aromatic ring? */
function in_aromatic_ring(ctx: ClassContext, i: number, j: number): boolean {
  const rings = ctx.aromatic_rings.get(i);
  if (!rings) return false;
  return rings.some(r => is_edge_of_ring(r, i, j));
}

/** Are i and j consecutive atoms of the ring's cycle? */
function is_edge_of_ring(ring: AromaticRing, i: number, j: number): boolean {
  const p = ring.path;
  const idx = p.indexOf(i);
  if (idx < 0) return false;
  const prev = p[(idx - 1 + p.length) % p.length];
  const next = p[(idx + 1) % p.length];
  return j === prev || j === next;
}

/** Do i and j belong to DIFFERENT sets of aromatic rings? The part-V
 *  p. 620 clause (b) — "between pairs of atoms belonging to different
 *  aromatic rings" — is a SET comparison, not disjointness: an atom
 *  shared with the partner's ring is fine (naphthalene's fusion bond
 *  stays BTij = 0 because both atoms sit in the same two rings), and
 *  one side with NO aromatic ring counts as different (DAKCEX's
 *  37–63 bond: C8 in the benzene ring, C9 in none → conjugated).
 *  The ring objects are shared across atoms, so identity compares
 *  rings. */
function different_aromatic_rings(ctx: ClassContext, i: number, j: number): boolean {
  const ri = ctx.aromatic_rings.get(i) ?? [];
  const rj = ctx.aromatic_rings.get(j) ?? [];
  if (ri.length !== rj.length) return true;
  return ri.some(r => !rj.includes(r)) || rj.some(r => !ri.includes(r));
}

/** GetBondType — the BTij flag for bond (i, j), see the header.
 *
 *  Part V p. 620: BTij = 1 for a single bond (formal bond order 1)
 *  (a) between atoms of types that are not both aromatic and for
 *  which the sbmb flag is set in Table I, or (b) between pairs of
 *  atoms belonging to different aromatic rings (biphenyl's central
 *  bond). A bond inside an aromatic ring is an aromatic bond and
 *  reads 0 — including fused rings, where both atoms share the same
 *  aromatic-ring set (naphthalene). */
export function bond_type_flag(ctx: ClassContext, i: number, j: number): number {
  const { mol } = ctx;
  if (get_bond_order(ctx, i, j) !== 1) return 0;

  const pi = ATOM_TYPE_PROPERTIES[mol.atom_types[i]];
  const pj = ATOM_TYPE_PROPERTIES[mol.atom_types[j]];
  if (!pi || !pj) return 0;

  // CIM+ (80) is an aromatic ring type whose par entry lacks the arom
  // flag — treat it as aromatic so its ring bonds read class 0 (the
  // BCI table keys them 0-80-81).
  const a_i = pi.arom || mol.atom_types[i] === 80;
  const a_j = pj.arom || mol.atom_types[j] === 80;

  // An aromatic ring bond is never a conjugated single bond.
  if (in_aromatic_ring(ctx, i, j)) return 0;

  if (a_i && a_j) {
    // Both aromatic-flagged: conjugated only across DIFFERENT
    // aromatic rings (biphenyl's central bond, biphenylene's
    // 4-ring bonds). Same ring, no rings, or one side ringless:
    // 0 — e.g. two acyclic aromatic-typed N's (DAKCEX's N1–N2).
    return different_aromatic_rings(ctx, i, j) ? 1 : 0;
  }
  // Case (a): both types have the sbmb flag — what makes C(=O)-C(ar)
  // and C=C-C(ar) angles class 1/2 and conjugated single bonds
  // class 1.
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

  // 4-ring: the closing bond i-l makes i-j-k-l a 4-cycle. Checked
  // BEFORE the BT class: the reference's class for a 4-ring torsion
  // with a conjugated central bond is 4 (CEWYIM30's C1-C6-C12-C7 —
  // BT(C6-C12) = 1, class 4, V2 = 6.0), not 1.
  if (ctx.adj[i].includes(l)) return 4;

  if (bt_bc === 1) return 1;

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
    if (params.k_a === 0) {
      // A found k_a = 0 row is an out-of-range default (e.g. the
      // '0-0-26-0' phosphine wildcard): the references do NOT treat
      // it as "no contribution". OpenBabel's setup log ("USING
      // EMPIRICAL RULE FOR ANGLE BENDING FORCE CONSTANT 1-3-7") and
      // Tinker's parameter listing (KB = 0.661 @ 98.1° for the vinyl
      // phosphine C–P–H) both apply the part V empirical force
      // constant (eq. 20) with the row's θ₀. The old "found k = 0 is
      // final" policy (from the sulfinate S=O angles, BatchMin's 0)
      // zeroed the C–P–H and C–C–P terms entirely, so nothing held
      // the P pyramidal and it flattened to trigonal planar.
      theta0 = params.theta0;
      const r0ab =
        bond_parameters(ctx, i, j)?.r0 ??
        empirical_bond_length(mol.atoms[i], mol.atoms[j]) ??
        1.5;
      const r0bc =
        bond_parameters(ctx, j, k)?.r0 ??
        empirical_bond_length(mol.atoms[j], mol.atoms[k]) ??
        1.5;
      k_a = empirical_ka(mol.atoms[i], mol.atoms[j], mol.atoms[k], r0ab, r0bc, theta0, cls);
    } else {
      k_a = params.k_a;
      theta0 = params.theta0;
    }
  } else {
    // Total miss: the part V empirical θ₀ and force-constant rules
    // (empirical.ts). The eq. (20) reference bond lengths fall back
    // from the par to the eq. (18) empirical values (the Tinker
    // behavior) before the 1.5 Å stand-in.
    theta0 = empirical_theta0(ATOM_TYPE_PROPERTIES[tj], mol.atoms[j].element, cls);
    const r0ab =
      bond_parameters(ctx, i, j)?.r0 ??
      empirical_bond_length(mol.atoms[i], mol.atoms[j]) ??
      1.5;
    const r0bc =
      bond_parameters(ctx, j, k)?.r0 ??
      empirical_bond_length(mol.atoms[j], mol.atoms[k]) ??
      1.5;
    k_a = empirical_ka(mol.atoms[i], mol.atoms[j], mol.atoms[k], r0ab, r0bc, theta0, cls);
  }

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

  // The order index is the canonical-direction discriminator. M = 136
  // is one more than the largest MMFF94 atom type number (135), so
  // each type occupies a digit in a base-136 mixed-radix encoding of
  // the i-j-k-l direction; order >= 0 selects the lexicographically
  // canonical direction, matching the direction under which the par
  // stores its asymmetric torsion entries.
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
