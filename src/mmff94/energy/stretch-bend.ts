/**
 * Stretch-bend cross term energy.
 *
 * This is a CLASS II force field term — it couples bond stretching with
 * angle bending. Most force fields (UFF, GAFF, MM2/3) omit this term, but
 * MMFF94 includes it because bond lengths and angles are physically coupled:
 * when an H-C-H angle in methane closes, the C-H bonds shorten slightly.
 *
 *   E_sb = 2.51210 * k_sb * [(r₁ − r₁₀) + (r₂ − r₂₀)] * (θ − θ₀)
 *
 * where:
 *   k_sb  = stretch-bend force constant
 *   r₁    = length of bond i−j (one side of the angle)
 *   r₂    = length of bond j−k (the other side)
 *   r₁₀   = equilibrium length of bond i−j
 *   r₂₀   = equilibrium length of bond j−k
 *   θ     = current angle i−j−k
 *   θ₀    = equilibrium angle i−j−k
 *   2.51210 = unit conversion factor
 *
 * The three parameters (k_sb, k_b, r₀) for the two bonds and the one
 * angle all come from the same angle parameter entry (type_i, type_j, type_k).
 *
 * Note: the equilibrium bond lengths r₁₀ and r₂₀ come from the BOND
 * parameters for the respective bond types, NOT from the stretch-bend table.
 * The stretch-bend table only provides k_sb.
 */

import type { TypedMolecule } from '../../types';

/**
 * Calculate the total stretch-bend cross term energy.
 */
export function calc_stretch_bend_energy(molecule: TypedMolecule): number {
  // TODO: implement stretch-bend energy.
  //
  // This term requires both bond stretching and angle bending parameters.
  // It must be calculated AFTER bond and angle parameters are loaded.
  //
  // For each angle i−j−k:
  //   1. Look up k_sb from STRETCH_BEND_PARAMS by (type_i, type_j, type_k).
  //   2. Look up (k_b, r₁₀) for bond i−j from BOND_PARAMS.
  //   3. Look up (k_b, r₂₀) for bond j−k from BOND_PARAMS.
  //   4. Look up θ₀ from ANGLE_PARAMS for angle i−j−k.
  //   5. Compute r₁ = |pos_i − pos_j|, r₂ = |pos_j − pos_k|, θ = angle(i,j,k).
  //   6. Accumulate: E += 2.51210 * k_sb * [(r₁ − r₁₀) + (r₂ − r₂₀)] * (θ − θ₀)
  //
  // Return the sum in kcal/mol.
  return 0.0;
}
