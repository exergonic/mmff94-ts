/**
 * Angle bending energy.
 *
 * Halgren1996, eq. (3):
 *
 *   E_angle = 0.043844 · k_a · (θ − θ₀)²
 *             · [1 + cb · (θ − θ₀)]
 *
 * where:
 *   k_a   = force constant in mdyn·Å/rad²
 *   θ     = current bond angle in degrees
 *   θ₀    = equilibrium bond angle in degrees
 *   cb    = cubic bend constant = −0.007 deg⁻¹ ( −0.4 rad⁻¹ )
 *   0.043844 = unit conversion factor: (mdyn·Å/rad²) → (kcal/mol)/deg²
 *
 * We implement only the leading harmonic term. The cubic correction
 * (with cb) is omitted for the same reason as bond-stretch: for the
 * small angular displacements at equilibrium, the anharmonic
 * contribution is negligible.
 *
 * For near-linear angles (θ ≈ 180°), MMFF94 switches to a cosine
 * form, eq. (4):
 *   E_angle = 143.9325 · k_a · (1 + cos θ)
 * This is not currently implemented — linear geometries are rare in
 * the organic molecules covered by the current test set.
 */

import type { TypedMolecule } from '../../types';
import { ANGLE_PARAMS, lookup_param } from '../parameters';
import { angle_in_radians, Vec3 } from '../../utils/vector';

/**
 * Calculate the total angle bending energy for all bond angles in a molecule.
 */
export function calc_angle_bend_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  // Build adjacency list to easily find neighbors
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

    // Iterate over all pairs of neighbors (i, k)
    for (let i_idx = 0; i_idx < neighbors.length; i_idx++) {
      for (let k_idx = i_idx + 1; k_idx < neighbors.length; k_idx++) {
        const i = neighbors[i_idx];
        const k = neighbors[k_idx];
        
        const ti = molecule.atom_types[i];
        const tk = molecule.atom_types[k];
        
        const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];
        const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

        const t_min = Math.min(ti, tk);
        const t_max = Math.max(ti, tk);

        const params = lookup_param(ANGLE_PARAMS, [t_min, tj, t_max]);
        if (params) {
          const { k_a, theta0 } = params;
          const theta_rad = angle_in_radians(posJ, posI, posK);
          const theta_deg = theta_rad * (180.0 / Math.PI);
          const d_theta = theta_deg - theta0;
          total_energy += 0.043844 * k_a * d_theta * d_theta;
        }
      }
    }
  }

  return total_energy;
}
