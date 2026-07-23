/**
 * Stretch-bend cross term energy.
 *
 * This is a CLASS II force field term — it couples bond stretching with
 * angle bending. Most force fields (UFF, GAFF, MM2/3) omit this term, but
 * MMFF94 includes it because bond lengths and angles are physically coupled:
 * when an H-C-H angle in methane closes, the C-H bonds shorten slightly.
 *
 *   E_sb = 2.51210 * [k_sb_IJK * (r_IJ - r_IJ0) + k_sb_KJI * (r_KJ - r_KJ0)]
 *                   * (theta - theta0)
 *
 * where:
 *   k_sb_IJK = stretch-bend force constant for the I-J bond side
 *   k_sb_KJI = stretch-bend force constant for the K-J bond side
 *   r_IJ     = length of bond i-j (one side of the angle)
 *   r_KJ     = length of bond k-j (the other side)
 *   r_IJ0    = equilibrium length of bond i-j (from BOND_PARAMS)
 *   r_KJ0    = equilibrium length of bond k-j (from BOND_PARAMS)
 *   theta    = current angle i-j-k
 *   theta0   = equilibrium angle i-j-k (from ANGLE_PARAMS)
 *   2.51210  = unit conversion factor
 *
 * The two k_sb values are indexed by (type_i, type_j, type_k). Many entries
 * have k_sb_IJK == k_sb_KJI (symmetric angles), but asymmetric environments
 * (e.g., C-C-O) can have different values.
 */

import type { TypedMolecule } from '../../types';
import { BOND_PARAMS, ANGLE_PARAMS, STRETCH_BEND_PARAMS, lookup_param } from '../parameters';
import { distance, angle_in_radians, Vec3 } from '../../utils/vector';

/**
 * Calculate the total stretch-bend cross term energy.
 */
export function calc_stretch_bend_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  // Build adjacency list
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  // Iterate over all possible central atoms j
  for (let j = 0; j < molecule.atoms.length; j++) {
    const neighbors = adj[j];
    const tj = molecule.atom_types[j];
    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];

    // For each pair of neighbors (i, k) — these define angles i-j-k
    for (let i_idx = 0; i_idx < neighbors.length; i_idx++) {
      for (let k_idx = i_idx + 1; k_idx < neighbors.length; k_idx++) {
        const i = neighbors[i_idx];
        const k = neighbors[k_idx];

        const ti = molecule.atom_types[i];
        const tk = molecule.atom_types[k];

        const t_min = Math.min(ti, tk);
        const t_max = Math.max(ti, tk);

        // 1. Look up stretch-bend parameters
        const sb_params = lookup_param(STRETCH_BEND_PARAMS, [t_min, tj, t_max]);
        if (!sb_params) continue;

        // 2. Look up equilibrium bond lengths from bond parameters
        const bond_ij = lookup_param(BOND_PARAMS, [t_min, tj]);
        const bond_kj = lookup_param(BOND_PARAMS, [t_max, tj]);
        if (!bond_ij || !bond_kj) continue;

        // 3. Look up equilibrium angle from angle parameters
        const ang_params = lookup_param(ANGLE_PARAMS, [t_min, tj, t_max]);
        if (!ang_params) continue;

        // 4. Compute current geometry
        const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];
        const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

        const r_IJ = distance(posI, posJ);
        const r_KJ = distance(posK, posJ);
        const theta_rad = angle_in_radians(posJ, posI, posK);
        const theta_deg = theta_rad * (180.0 / Math.PI);

        // 5. Accumulate energy
        const dr_IJ = r_IJ - bond_ij.r0;
        const dr_KJ = r_KJ - bond_kj.r0;
        const d_theta = theta_deg - ang_params.theta0;

        total_energy += 2.51210 * (sb_params.k_sb_IJK * dr_IJ + sb_params.k_sb_KJI * dr_KJ) * d_theta;
      }
    }
  }

  return total_energy;
}
