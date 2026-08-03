/**
 * Gradient of the van der Waals energy.
 *
 * See energy/van-der-waals.ts for the energy — the buffered 14-7
 * potential, Halgren1996 eq. (8):
 *
 *   E_vdw = ε · f_rep(r) · f_att(r)
 *   f_rep = (a / (r + b))⁷        a = 1.07·R*, b = 0.07·R*
 *   f_att = C / (r⁷ + D) − 2      C = 1.12·R*⁷, D = 0.12·R*⁷
 *
 * whose derivative with respect to the interatomic distance is
 *
 *   f_rep' = −7·a⁷ / (r + b)⁸
 *   f_att' = −7·C·r⁶ / (r⁷ + D)²
 *
 *   dE/dr = ε · (f_rep'·f_att + f_rep·f_att')
 *
 * and dE/dx = dE/dr · dr/dx along the bond direction. The pair
 * parameters (R*, ε) come from vdw_pair_parameters() — the same
 * combination rules the energy term uses — and the pair enumeration
 * is identical: every i < j pair, excluding 1-2 and 1-3, with NO
 * 1-4 scaling (Halgren: "1,4-vdW interactions are not differentially
 * scaled in MMFF94").
 */

import type { TypedMolecule } from '../../types';
import { Vec3 } from '../../utils/vector';
import { vdw_pair_parameters, vdw_parameters_for } from '../energy/van-der-waals';
import { bond_length_derivatives } from './derivatives';

/**
 * Gradient of the vdW energy, dE/dx per atom.
 * Returns an array parallel to molecule.atoms[].
 */
export function calc_vdw_gradient(molecule: TypedMolecule): number[][] {
  const gradient: number[][] = molecule.atoms.map(() => [0, 0, 0]);

  // Build adjacency and 2-away sets for pair exclusions — same as the
  // energy term.
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  // Build 1-3 pair map: atoms sharing a common neighbor (angle pairs)
  const pairs_1_3: Set<number>[] = Array.from({ length: molecule.atoms.length }, () => new Set());
  for (let i = 0; i < molecule.atoms.length; i++) {
    for (const n1 of adj[i]) {
      for (const n2 of adj[n1]) {
        if (n2 !== i) pairs_1_3[i].add(n2);
      }
    }
  }

  for (let i = 0; i < molecule.atoms.length; i++) {
    const param_i = vdw_parameters_for(molecule, adj, i);
    if (!param_i) continue;

    const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];

    for (let j = i + 1; j < molecule.atoms.length; j++) {

      // Skip 1-2 (bonded) and 1-3 (share a common neighbor); 1-4 pairs
      // are NOT scaled in MMFF94.
      if (adj[i].includes(j)) continue;
      if (pairs_1_3[i].has(j)) continue;

      const param_j = vdw_parameters_for(molecule, adj, j);
      if (!param_j) continue;

      const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];

      const { R_ij, epsilon_ij } = vdw_pair_parameters(param_i, param_j);
      const r = Math.hypot(posI[0] - posJ[0], posI[1] - posJ[1], posI[2] - posJ[2]);

      // Buffered 14-7 expression and its r-derivative
      const a = 1.07 * R_ij;
      const b = 0.07 * R_ij;
      const C = 1.12 * Math.pow(R_ij, 7);
      const D = 0.12 * Math.pow(R_ij, 7);

      const r7 = Math.pow(r, 7);
      const f_rep = Math.pow(a / (r + b), 7);
      const f_att = C / (r7 + D) - 2;
      const f_rep_prime = -7.0 * Math.pow(a, 7) / Math.pow(r + b, 8);
      const f_att_prime = -7.0 * C * Math.pow(r, 6) / Math.pow(r7 + D, 2);

      const dE_dr = epsilon_ij * (f_rep_prime * f_att + f_rep * f_att_prime);

      const { d_dx_a, d_dx_b } = bond_length_derivatives(posI, posJ);
      for (let axis = 0; axis < 3; axis++) {
        gradient[i][axis] += dE_dr * d_dx_a[axis];
        gradient[j][axis] += dE_dr * d_dx_b[axis];
      }
    }
  }

  return gradient;
}
