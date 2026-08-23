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

import type { TypedMolecule } from '../../types.js';
import { Vec3 } from '../../utils/vector.js';
import { nonbonded_context_for } from '../nonbonded-context.js';
import { bond_length_derivatives } from './derivatives.js';

/**
 * Gradient of the vdW energy, dE/dx per atom.
 * Returns an array parallel to molecule.atoms[].
 */
export function calc_vdw_gradient(molecule: TypedMolecule): number[][] {
  const gradient: number[][] = molecule.atoms.map(() => [0, 0, 0]);

  // Pair list + combined parameters come from the cached context —
  // they are topology, not geometry (see nonbonded-context.ts).
  const ctx = nonbonded_context_for(molecule);

  for (let p = 0; p < ctx.n_pairs; p++) {
    const epsilon_ij = ctx.pair_epsilon_ij[p];
    if (isNaN(epsilon_ij)) continue; // an atom without vdW parameters

    const i = ctx.pair_i[p];
    const j = ctx.pair_j[p];
    const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];
    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];

    const r = Math.hypot(posI[0] - posJ[0], posI[1] - posJ[1], posI[2] - posJ[2]);

    // Buffered 14-7 expression and its r-derivative, with the fused
    // per-pair constants a/b/C/D from the context.
    const a = ctx.pair_vdw_a[p];
    const b = ctx.pair_vdw_b[p];
    const C = ctx.pair_vdw_C[p];
    const D = ctx.pair_vdw_D[p];

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

  return gradient;
}
