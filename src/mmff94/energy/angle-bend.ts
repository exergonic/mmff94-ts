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
 * For linear centers (the lin flag in mmffprop.par, e.g. sp carbon),
 * eq. (3) is replaced by eq. (4):
 *   E_angle = 143.9325 · k_a · (1 + cos θ)
 * which avoids the singularity of the cubic expansion at θ = 180°.
 * The linear flag comes from the atom-type properties, not from θ₀.
 */

import type { TypedMolecule } from '../../types';
import { angle_in_radians, Vec3 } from '../../utils/vector';
import { make_class_context, angle_parameters } from '../parameters/parameter-classes';

// Cubic bend constant. The paper gives cb = −0.007 deg⁻¹ "(or, more
// precisely, −0.4 rad⁻¹)" — BatchMin uses the precise radian value
// converted to degrees. The rounding matters for large deviations:
// GESNIB's near-linear C(37)–C(37)–C(22) angles (Δ ≈ 47°) differ by
// 0.035 kcal/mol each between the two forms.
const CB = -0.4 * (Math.PI / 180); // −0.0069813… deg⁻¹

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
  const ctx = make_class_context(molecule, adj);

  // Iterate over all possible central atoms j
  for (let j = 0; j < molecule.atoms.length; j++) {
    const neighbors = adj[j];
    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];

    // Iterate over all pairs of neighbors (i, k)
    for (let i_idx = 0; i_idx < neighbors.length; i_idx++) {
      for (let k_idx = i_idx + 1; k_idx < neighbors.length; k_idx++) {
        const i = neighbors[i_idx];
        const k = neighbors[k_idx];

        const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];
        const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

        // Class-aware resolution: BTij/ring classes select the entry
        // (small rings have their own θ₀ ≈ 90°/60°); misses fall to the
        // empirical rules. Linear centers (the lin flag, e.g. sp C) use
        // the cosine form regardless of the lookup result.
        const { k_a, theta0, linear } = angle_parameters(ctx, i, j, k);
        const theta_rad = angle_in_radians(posJ, posI, posK);

        if (linear) {
          // Halgren1996 eq. (4): cosine form for linear centers
          total_energy += 143.9325 * k_a * (1.0 + Math.cos(theta_rad));
        } else {
          // Halgren1996 eq. (3): harmonic + cubic expansion
          const theta_deg = theta_rad * (180.0 / Math.PI);
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
