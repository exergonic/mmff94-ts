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
 *              compute_bci_charges(), or the term computes them on
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

import type { TypedMolecule } from '../../types';
import { distance, Vec3 } from '../../utils/vector';
import { compute_bci_charges } from '../charges';

const S = 0.05; // the electrostatic buffering constant (Å)

/**
 * Calculate the total electrostatic energy.
 */
export function calc_electrostatic_energy(molecule: TypedMolecule): number {
  // Partial charges from the BCI model — attached by the caller's
  // compute_bci_charges step (the value flowing down the pipeline), or
  // computed here if the term is called alone.
  const charges =
    molecule.partial_charges ?? compute_bci_charges(molecule).partial_charges!;

  // Adjacency: for the 1-2/1-3 exclusion and the 1-4 classification.
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  // 1-3 pair map: atoms sharing a common neighbor (angle pairs)
  const pairs_1_3: Set<number>[] = Array.from({ length: molecule.atoms.length }, () => new Set());
  for (let i = 0; i < molecule.atoms.length; i++) {
    for (const n1 of adj[i]) {
      for (const n2 of adj[n1]) {
        if (n2 !== i) pairs_1_3[i].add(n2);
      }
    }
  }

  const n = molecule.atoms.length;
  const positions: Vec3[] = molecule.atoms.map(a => [a.x, a.y, a.z]);
  let total_energy = 0.0;

  for (let i = 0; i < n; i++) {
    const qi = charges[i];
    if (qi === 0) continue;
    for (let j = i + 1; j < n; j++) {
      // Skip 1-2 (bonded) and 1-3 (share a common neighbor) pairs
      if (adj[i].includes(j)) continue;
      if (pairs_1_3[i].has(j)) continue;

      const qj = charges[j];
      if (qj === 0) continue;

      const r = distance(positions[i], positions[j]) + S;

      let energy = (332.0716 * qi * qj) / r; // D = 1.0 (in vacuo)

      // 1-4 scaling: atoms exactly three bonds apart.
      if (is_1_4_pair(i, j, adj)) energy *= 0.75;

      total_energy += energy;
    }
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
    if (node === j && depth === 3) return true;
    if (depth >= 3 || seen.has(node)) continue;
    seen.add(node);
    for (const nb of adj[node]) {
      if (!seen.has(nb)) queue.push([nb, depth + 1]);
    }
  }
  return false;
}
