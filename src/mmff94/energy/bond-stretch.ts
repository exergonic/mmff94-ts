/**
 * Bond stretching energy.
 *
 * MMFF94 uses a harmonic potential for bond stretching:
 *
 *   E_bond = 143.88 * k_b * (r − r₀)²
 *
 * where:
 *   k_b  = force constant in millidynes per Ångström (mdyn/Å)
 *   r    = current bond length in Å
 *   r₀   = equilibrium bond length in Å
 *   143.88 = unit conversion factor: (mdyn/Å) → (kcal/mol)/Å²
 *
 * A harmonic approximation is used instead of a full Morse potential
 * because for most organic molecules at room temperature, bonds do
 * not stretch far enough from equilibrium to feel the anharmonicity.
 * MMFF94 is parametrized for equilibrium geometries and conformational
 * energies, not bond dissociation.
 *
 * The parameters (k_b, r₀) are indexed by the MMFF94 types of the two
 * bonded atoms. If no exact (type_i, type_j) match is found, the lookup
 * falls through to a wildcard entry if one exists.
 */

import type { TypedMolecule } from '../../types';

// Placeholder: import { BOND_PARAMS } from '../parameters/bond';

/**
 * Calculate the total bond stretching energy for all bonds in a molecule.
 */
export function calc_bond_stretch_energy(molecule: TypedMolecule): number {
  // TODO: implement bond stretching energy calculation.
  //
  // For each bond in molecule.bonds:
  //   1. Get the two atoms and their current positions.
  //   2. Compute r = distance between them.
  //   3. Look up (k_b, r₀) from BOND_PARAMS by (type_i, type_j).
  //      - Try (type_i, type_j) in both orders.
  //      - Fall back to wildcard '*' if no exact match.
  //   4. Accumulate: E += 143.88 * k_b * (r - r₀)²
  //
  // Return the sum in kcal/mol.
  return 0.0;
}
