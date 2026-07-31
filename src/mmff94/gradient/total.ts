/**
 * Gradient of the total MMFF94 energy with respect to atomic positions.
 *
 * Returns dE/dx, dE/dy, dE/dz for every atom as a 2D array:
 *   gradient[i] = [dE/dx_i, dE/dy_i, dE/dz_i]
 *
 * The gradient is the negative of the force on each atom:
 *   F_i = −∇_i E
 *
 * STATUS: NOT YET IMPLEMENTED — returns a zero gradient for every atom.
 * The plan (AGENTS.md Phase 5) is one gradient file per energy term,
 * mirroring the energy/ layout, with finite-difference cross-checks in
 * tests/gradient.test.ts.
 */

import type { TypedMolecule } from '../../types';

/**
 * Compute the full gradient of the MMFF94 energy.
 *
 * Returns an array parallel to molecule.atoms[]:
 *   result[i] = [dE/dx_i, dE/dy_i, dE/dz_i]
 * with units of kcal/mol/Å.
 */
export function calc_gradient(molecule: TypedMolecule): number[][] {
  // TODO: implement full MMFF94 gradient.
  //
  // For each energy term:
  //   1. Call the corresponding gradient function (to be created).
  //   2. Sum the contributions into the per-atom gradient array.
  //
  // Return: gradient[i] = [dE/dx_i, dE/dy_i, dE/dz_i].
  return molecule.atoms.map(() => [0, 0, 0]);
}
