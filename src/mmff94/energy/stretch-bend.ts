/**
 * Stretch-bend cross term energy.
 *
 * Halgren1996, eq. (5):
 *
 *   E_sb = 2.51210 · [k_sb_IJK · (r_IJ − r_IJ0)
 *                     + k_sb_KJI · (r_KJ − r_KJ0)]
 *                    · (θ − θ₀)
 *
 * This is a CLASS II force field term — it couples bond stretching with
 * angle bending. Most force fields (UFF, GAFF, MM2/3) omit this term, but
 * MMFF94 includes it because bond lengths and angles are physically coupled:
 * when an H-C-H angle in methane closes, the C-H bonds shorten slightly.
 *
 * Two separate k_sb values are used because asymmetric environments
 * (e.g., C-C-O vs O-C-C) have different coupling strengths for the two
 * sides of the angle. Many entries have k_sb_IJK = k_sb_KJI (symmetric
 * angles), but the parameter table stores both independently.
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

        // The table stores terminal types sorted (I ≤ K), so k_sb_IJK
        // belongs to the bond on the min-type side and k_sb_KJI to the
        // bond on the max-type side — whichever of i and k that is.
        const k_ij = ti <= tk ? sb_params.k_sb_IJK : sb_params.k_sb_KJI;
        const k_kj = ti <= tk ? sb_params.k_sb_KJI : sb_params.k_sb_IJK;

        // 2. Look up equilibrium bond lengths from bond parameters —
        //    each bond uses its own type pair (sorted), not the angle's
        //    sorted terminal types.
        const bond_ij = lookup_param(BOND_PARAMS, [Math.min(ti, tj), Math.max(ti, tj)]);
        const bond_kj = lookup_param(BOND_PARAMS, [Math.min(tk, tj), Math.max(tk, tj)]);
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

        total_energy += 2.51210 * (k_ij * dr_IJ + k_kj * dr_KJ) * d_theta;
      }
    }
  }

  return total_energy;
}
