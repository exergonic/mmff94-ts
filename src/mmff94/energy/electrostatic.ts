/**
 * Electrostatic (Coulombic) energy.
 *
 * MMFF94 uses Coulomb's law with partial charges derived from the
 * bond charge increment (BCI) model:
 *
 *   E_elec = 332.0716 * q_i * q_j / (ε * r)
 *
 * where:
 *   q_i, q_j = partial charges on atoms i and j (e⁻, from BCI)
 *   r        = distance between atoms (Å)
 *   ε        = dielectric constant
 *   332.0716 = unit conversion factor: (e⁻)²/Å → kcal/mol
 *
 * The dielectric constant ε depends on the environment:
 *   ε = 1       for in-vacuo calculations
 *   ε = r       for distance-dependent dielectric (MMFF94 default)
 *   ε = 4.0     for protein/interior calculations
 *
 * The BCI model assigns a charge increment to each bond type.
 * Each atom's partial charge is the sum of the increments from
 * all bonds it participates in. This means partial charges depend
 * on connectivity and atom types, not on fixed per-atom tables.
 *
 * The BCI values are stored in src/mmff94/parameters/bci.ts and
 * the charge computation lives in src/mmff94/atom-types.ts.
 *
 * 1-4 SCALING: For atoms exactly three bonds apart, the electrostatic
 * energy is multiplied by 0.75. Applied in total.ts.
 */

import type { TypedMolecule } from '../../types';

/**
 * Calculate the total electrostatic energy between all non-bonded atom pairs.
 *
 * Requires that molecule.partial_charges is non-null
 * (call compute_bci_charges() first).
 *
 * Excludes 1-2 and 1-3 pairs. 1-4 pairs are calculated with full
 * electrostatic energy; the 0.75 scaling factor is applied in total.ts.
 */
export function calc_electrostatic_energy(molecule: TypedMolecule): number {
  // TODO: implement electrostatic energy.
  //
  // This function must run AFTER compute_bci_charges() fills in
  // molecule.partial_charges[].
  //
  // For each pair of atoms (i, j) where i < j:
  //   1. Skip if 1-2 or 1-3 pair.
  //   2. Compute r = distance.
  //   3. Look up partial charges q_i, q_j.
  //   4. Choose ε (1.0 for in-vacuo, r for distance-dependent).
  //   5. Accumulate: E += 332.0716 * q_i * q_j / (ε * r).
  //   6. (1-4 scaling applied in total.ts.)
  //
  // Return the sum in kcal/mol (before 1-4 scaling).
  return 0.0;
}
