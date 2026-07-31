/**
 * Out-of-plane bending energy.
 *
 * MMFF94 uses a dedicated out-of-plane bending term for sp²-hybridized
 * (trigonal planar) centers, rather than the "improper torsion" approach
 * used by UFF and GAFF. These are NOT the same thing:
 *
 *   An improper torsion (UFF/GAFF) treats the oop deformation as a
 *   dihedral rotation, which couples it with the true torsion term.
 *
 *   MMFF94's oop term is a PURE out-of-plane bending — it measures
 *   the distance of the central atom from the plane defined by its
 *   three substituents, independent of any dihedral rotation.
 *
 *   E_oop = 0.043844 * k_oop/2 * χ²
 *
 * where:
 *   k_oop  = out-of-plane force constant
 *   χ      = out-of-plane deformation angle (degrees)
 *           (the angle between the vector from the central atom to
 *            any one substituent and the plane defined by the other
 *            two substituents and the central atom)
 *   0.043844 = unit conversion factor: same as angle bending
 *
 * This term applies to:
 *   - Carbonyl carbon (C=O)
 *   - Olefinic carbon (C=C, both sp² carbons)
 *   - Aromatic carbon (in benzene rings)
 *   - Trigonal planar nitrogen (amide N, pyridine N)
 *   - Trigonal planar oxygen (carbonyl-like? — check parameter table)
 *
 * The parameter k_oop is indexed by the type of the CENTRAL atom only.
 */

import type { TypedMolecule } from '../../types';

/**
 * Calculate the total out-of-plane bending energy.
 *
 * Evaluated for every sp² center where the central atom has exactly
 * three bonded neighbors.
 */
export function calc_oop_energy(molecule: TypedMolecule): number {
  // TODO: implement out-of-plane bending energy.
  //
  // For each atom j (the central atom) with exactly 3 bonded neighbors
  // (i, k, l):
  //   1. Check if atom j's type is a known sp² center:
  //      - Look at the atom type, or check by element and coordination.
  //   2. Compute χ: the oop deformation angle.
  //      - Let (i, k, l) be the three neighbors of j.
  //      - Compute the plane through (i, k, l).
  //      - Compute the angle between the normal of that plane and the
  //        vector from the plane's centroid to j.
  //      - Alternatively: the Wilson out-of-plane angle.
  //   3. Look up k_oop from OOP_PARAMS by atom type of j.
  //   4. Accumulate: E += 0.043844 * k_oop/2 * χ².
  //
  // Return the sum in kcal/mol.
  return 0.0;
}
