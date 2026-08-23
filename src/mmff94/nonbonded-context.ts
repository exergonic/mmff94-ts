/**
 * Cached nonbonded pair context.
 *
 * The van der Waals and electrostatic terms iterate the same pair set —
 * every i < j pair except 1-2 (bonded) and 1-3 (angle) pairs — and every
 * quantity they need beyond the live coordinates is TOPOLOGY, not
 * geometry: the exclusion sets, the exactly-three-bonds flag (the
 * electrostatic 1-4 ×0.75), the combined vdW parameters (Waldman-Hagler
 * R*ij / Slater-Kirkwood εij depend only on the two atoms' parameter
 * rows), and the charge products. Rebuilding all of that inside every
 * energy/gradient call — including a per-pair BFS for the 1-4 test —
 * was the dominant cost of optimization, where the terms run thousands
 * of times on a molecule whose atoms/bonds/types/charges never change.
 *
 * This module builds that context ONCE per molecule object and caches
 * it (WeakMap, garbage-collected with the molecule). The cache contract
 * mirrors the optimizer's: COORDINATES may move between calls; the
 * atom/bond list, atom types, formal charges, and partial charges may
 * not. Mutating any of those in place invalidates nothing here — pass a
 * fresh molecule (or drop this module's entry by letting the molecule
 * go out of scope) after changing chemistry.
 *
 * The pair order produced by build() matches the historical inline
 * loops exactly (i ascending, j ascending within i, zero-charge pairs
 * omitted for electrostatics), so summation order — and therefore the
 * floating-point result bit-for-bit — is unchanged.
 */

import type { TypedMolecule } from '../types.js';
import { vdw_parameters_for, vdw_pair_parameters } from './energy/van-der-waals.js';

/** One nonbonded pair, flattened into parallel arrays for the hot loop. */
export interface NonbondedContext {
  /** Atom count (for gradient array sizing). */
  n_atoms: number;
  /** Number of listed pairs. */
  n_pairs: number;
  pair_i: Int32Array;
  pair_j: Int32Array;
  /** Exactly-three-bonds flag — the electrostatic ×0.75 set (vdW never scales). */
  pair_is_14: Uint8Array;
  /** Combined vdW parameters per pair (NaN rows where either atom lacks parameters). */
  pair_R_ij: Float64Array;
  pair_epsilon_ij: Float64Array;
  /** Pre-fused constants of the buffered 14-7 form: a = 1.07·R*, b = 0.07·R*, C = 1.12·R*⁷, D = 0.12·R*⁷. */
  pair_vdw_a: Float64Array;
  pair_vdw_b: Float64Array;
  pair_vdw_C: Float64Array;
  pair_vdw_D: Float64Array;
}

const cache = new WeakMap<TypedMolecule, NonbondedContext>();

/**
 * The context for `molecule`, built on first use. See the module header
 * for the mutation contract.
 */
export function nonbonded_context_for(molecule: TypedMolecule): NonbondedContext {
  let ctx = cache.get(molecule);
  if (!ctx) {
    ctx = build_nonbonded_context(molecule);
    cache.set(molecule, ctx);
  }
  return ctx;
}

function build_nonbonded_context(molecule: TypedMolecule): NonbondedContext {
  const n = molecule.atoms.length;

  // Adjacency (1-2)
  const adj: number[][] = Array.from({ length: n }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  // 1-3 pairs: atoms sharing a common neighbor
  const pairs_1_3: Set<number>[] = Array.from({ length: n }, () => new Set());
  for (let i = 0; i < n; i++) {
    for (const n1 of adj[i]) {
      for (const n2 of adj[n1]) {
        if (n2 !== i) pairs_1_3[i].add(n2);
      }
    }
  }

  // Exactly-depth-3 reachability, ONE bounded BFS per atom instead of
  // one per pair (this replaces the hot-path is_1_4_pair calls; same
  // shortest-path-wins semantics: first arrival at depth 3).
  const reach_1_4: Set<number>[] = Array.from({ length: n }, () => new Set());
  for (let i = 0; i < n; i++) {
    const seen = new Set<number>([i]);
    let frontier = adj[i].slice();
    for (let depth = 1; depth <= 3 && frontier.length > 0; depth++) {
      const next: number[] = [];
      for (const node of frontier) {
        if (seen.has(node)) continue;
        seen.add(node);
        if (depth === 3) reach_1_4[i].add(node);
        else for (const nb of adj[node]) next.push(nb);
      }
      frontier = next;
    }
  }

  // Per-atom vdW parameter rows (topology-fixed: the sulfinate and
  // metal-oxidation bridges read types/formal charges/adjacency).
  const vdw_rows: ReturnType<typeof vdw_parameters_for>[] = new Array(n);
  for (let i = 0; i < n; i++) {
    vdw_rows[i] = vdw_parameters_for(molecule, adj, i);
  }

  // Pair enumeration — the same order as the historical i<j double loops.
  const pairs: number[] = [];   // flattened i, j
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (adj[i].includes(j)) continue;
      if (pairs_1_3[i].has(j)) continue;
      pairs.push(i, j);
    }
  }
  const n_pairs = pairs.length / 2;

  const ctx: NonbondedContext = {
    n_atoms: n,
    n_pairs,
    pair_i: new Int32Array(n_pairs),
    pair_j: new Int32Array(n_pairs),
    pair_is_14: new Uint8Array(n_pairs),
    pair_R_ij: new Float64Array(n_pairs),
    pair_epsilon_ij: new Float64Array(n_pairs),
    pair_vdw_a: new Float64Array(n_pairs),
    pair_vdw_b: new Float64Array(n_pairs),
    pair_vdw_C: new Float64Array(n_pairs),
    pair_vdw_D: new Float64Array(n_pairs),
  };

  for (let p = 0; p < n_pairs; p++) {
    const i = pairs[2 * p];
    const j = pairs[2 * p + 1];
    ctx.pair_i[p] = i;
    ctx.pair_j[p] = j;
    ctx.pair_is_14[p] = reach_1_4[i].has(j) ? 1 : 0;

    // Combined vdW parameters — pure functions of the two rows.
    const row_i = vdw_rows[i];
    const row_j = vdw_rows[j];
    if (row_i && row_j) {
      const { R_ij, epsilon_ij } = vdw_pair_parameters(row_i, row_j);
      ctx.pair_R_ij[p] = R_ij;
      ctx.pair_epsilon_ij[p] = epsilon_ij;
      ctx.pair_vdw_a[p] = 1.07 * R_ij;
      ctx.pair_vdw_b[p] = 0.07 * R_ij;
      ctx.pair_vdw_C[p] = 1.12 * Math.pow(R_ij, 7);
      ctx.pair_vdw_D[p] = 0.12 * Math.pow(R_ij, 7);
    } else {
      ctx.pair_R_ij[p] = NaN;
      ctx.pair_epsilon_ij[p] = NaN;
      ctx.pair_vdw_a[p] = NaN;
      ctx.pair_vdw_b[p] = NaN;
      ctx.pair_vdw_C[p] = NaN;
      ctx.pair_vdw_D[p] = NaN;
    }
  }

  return ctx;
}
