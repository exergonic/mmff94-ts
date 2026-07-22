/**
 * Gradient of the total MMFF94 energy with respect to atomic positions.
 *
 * Returns dE/dx, dE/dy, dE/dz for every atom as a 2D array:
 *   gradient[i] = [dE/dx_i, dE/dy_i, dE/dz_i]
 *
 * The gradient is the negative of the force on each atom:
 *   F_i = −∇_i E
 *
 * Each energy term contributes its own analytical derivative.
 * The gradient files in this directory mirror the energy/ layout:
 * one file per term, each exporting a function that returns the
 * gradient contribution for that term.
 *
 * Analytical gradients are derived from the same functional forms
 * documented in the corresponding energy/ files. They are
 * cross-checked against finite-difference calculations in the
 * test suite (tests/gradient.test.ts).
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
  // The gradient functions follow the same pattern as the energy functions:
  // one function per term, same signature, returning the gradient contribution.
  //
  // Return: gradient[i] = [dE/dx_i, dE/dy_i, dE/dz_i].
  return molecule.atoms.map(() => [0, 0, 0]);
}
