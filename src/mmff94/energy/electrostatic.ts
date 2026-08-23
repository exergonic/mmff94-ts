/**
 * Electrostatic (Coulombic) energy.
 *
 * Halgren1996, part III, eq. (6):
 *
 *   E_elec = 332.0716 · q_i · q_j / (D · (R_ij + S))
 *
 * where:
 *   q_i, q_j = partial charges (e⁻, from the BCI model — see
 *              src/mmff94/charges.ts; the caller attaches them via
 *              assign_bci_charges(), or the term computes them on
 *              the fly)
 *   R_ij     = interatomic distance (Å)
 *   S        = 0.05 Å — the ELECTROSTATIC BUFFERING CONSTANT, added
 *              to every distance to temper close contacts
 *   D        = dielectric constant — 1.0 in vacuo, the MMFF94 default
 *              (the distance-dependent D = r form is an alternative
 *              solvent model, not used here)
 *
 * The factor 332.0716 converts from e²/Å to kcal/mol.
 *
 * Only pairs separated by THREE or more bonds are evaluated — 1-2
 * (bonded) and 1-3 (one atom between) pairs are excluded, exactly as
 * in the van der Waals term (the suite's BatchMin energies confirm
 * it: ammonia's electrostatic energy is zero). The closest included
 * pairs, the 1-4 pairs (path i-j-k-l), are scaled by 0.75 because the
 * torsion term's parameterization already captures part of that
 * interaction. The scaling lives HERE rather than in total.ts because
 * the term functions return totals, not pair lists — total.ts cannot
 * rescale individual pairs.
 */

import type { TypedMolecule } from '../../types.js';
import { distance } from '../../utils/vector.js';
import { assign_bci_charges } from '../charges.js';
import { nonbonded_context_for } from '../nonbonded-context.js';

const S = 0.05; // the electrostatic buffering constant (Å)

/**
 * Calculate the total electrostatic energy.
 */
export function calc_electrostatic_energy(molecule: TypedMolecule): number {
  // Partial charges from the BCI model — attached by the caller's
  // assign_bci_charges step (the value flowing down the pipeline), or
  // computed here if the term is called alone.
  const charges =
    molecule.partial_charges ?? assign_bci_charges(molecule).partial_charges!;

  const ctx = nonbonded_context_for(molecule);
  let total_energy = 0.0;

  for (let p = 0; p < ctx.n_pairs; p++) {
    const i = ctx.pair_i[p];
    const j = ctx.pair_j[p];
    const qi = charges[i];
    const qj = charges[j];
    if (qi === 0 || qj === 0) continue; // zero-charge pairs contribute nothing

    const r = distance(
      [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z],
      [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z],
    ) + S;

    let energy = (332.0716 * qi * qj) / r; // D = 1.0 (in vacuo)

    // 1-4 scaling: atoms exactly three bonds apart.
    if (ctx.pair_is_14[p]) energy *= 0.75;

    total_energy += energy;
  }

  return total_energy;
}

/**
 * Are i and j connected by a path of exactly three bonds (i-j-k-l)?
 * Shared with the gradient so the 1-4 scaling is applied to exactly
 * the same pairs in both.
 */
export function is_1_4_pair(i: number, j: number, adj: number[][]): boolean {
  // BFS from i, bounded at depth 3 — j must appear at exactly that
  // depth (not sooner: 1-2 and 1-3 pairs are excluded entirely, and a
  // pair with both a 1-3 and a 1-4 path is a 1-3 pair).
  const seen = new Set<number>([i]);
  const queue: [number, number][] = adj[i].map(nb => [nb, 1]);
  let head = 0;
  while (head < queue.length) {
    const [node, depth] = queue[head++];
    // BFS pops level by level, so the FIRST arrival of j is its
    // shortest path: a 1-2 or 1-3 pair is never 1-4, even when a
    // longer three-bond path also exists (shortest path wins).
    if (node === j) return depth === 3;
    if (depth >= 3 || seen.has(node)) continue;
    seen.add(node);
    for (const nb of adj[node]) {
      if (!seen.has(nb)) queue.push([nb, depth + 1]);
    }
  }
  return false;
}
