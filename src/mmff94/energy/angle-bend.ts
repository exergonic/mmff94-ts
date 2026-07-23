/**
 * Angle bending energy.
 *
 * MMFF94 uses a harmonic potential for angle bending:
 *
 *   E_angle = 0.043844 * k_a * (θ − θ₀)²
 *
 * where:
 *   k_a   = force constant in mdyn·Å/rad²
 *   θ     = current bond angle in degrees
 *   θ₀    = equilibrium bond angle in degrees
 *   0.043844 = unit conversion factor: (mdyn·Å/rad²) → (kcal/mol)/deg²
 *
 * The angle is defined by three consecutive atoms: i − j − k,
 * where j is the central atom. The equilibrium angle θ₀ depends
 * on the types of all three atoms.
 *
 * Parameters are indexed by (type_i, type_j, type_k). The wildcard
 * lookup works the same as for bond stretching: exact match first,
 * then wildcard at terminal positions (i or k), then wildcard at
 * the central position (j).
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
