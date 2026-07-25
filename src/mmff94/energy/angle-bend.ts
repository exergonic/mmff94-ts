/**
 * Angle bending energy.
 *
 * Halgren1996, eq. (3):
 *
 *   E_angle = 0.043844 · (k_a / 2) · Δθ² · (1 + cb · Δθ)
 *
 * where:
 *   k_a   = force constant in mdyn·Å/rad²
 *   Δθ    = θ − θ₀, the deviation from the reference angle (degrees)
 *   θ₀    = equilibrium bond angle in degrees
 *   cb    = cubic bend constant = −0.007 deg⁻¹ ( −0.4 rad⁻¹ )
 *   0.043844 = unit conversion factor: (mdyn·Å/rad²) → (kcal/mol)/deg²
 *
 * The ½ factor is part of Halgren's definition (same as bond stretch).
 * The parameter table stores the full k_a; the ½ is applied here.
 *
 * For near-linear angles (θ₀ > 150°), eq. (3) is replaced by eq. (4):
 *   E_angle = 143.9325 · k_a · (1 + cos θ)
 * which avoids the singularity of the cubic expansion at θ = 180°.
 * This is checked via theta0, since linear reference angles define
 * the use of the cosine form regardless of the actual angle.
 */

import type { TypedMolecule } from '../../types';
import { ANGLE_PARAMS, lookup_param } from '../parameters';
import { angle_in_radians, Vec3 } from '../../utils/vector';

const CB = -0.007; // cubic bend constant, deg⁻¹

/**
 * Calculate the total angle bending energy for all bond angles in a molecule.
 */
export function calc_angle_bend_energy(molecule: TypedMolecule): number {
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
        if (!params) continue;

        const { k_a, theta0 } = params;
        const theta_rad = angle_in_radians(posJ, posI, posK);
        const theta_deg = theta_rad * (180.0 / Math.PI);

        if (theta0 > 150.0) {
          // Halgren1996 eq. (4): cosine form for near-linear angles
          total_energy += 143.9325 * k_a * (1.0 + Math.cos(theta_rad));
        } else {
          // Halgren1996 eq. (3): harmonic + cubic expansion
          const delta_theta = theta_deg - theta0;
          const half_k_a = 0.5 * k_a;
          const harmonic = 0.043844 * half_k_a * delta_theta * delta_theta;
          const anharmonic = 1.0 + CB * delta_theta;
          total_energy += harmonic * anharmonic;
        }
      }
    }
  }

  return total_energy;
}
