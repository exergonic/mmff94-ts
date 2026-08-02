/**
 * Gradient of the total MMFF94 energy with respect to atomic positions.
 *
 * Returns dE/dx, dE/dy, dE/dz for every atom as a 2D array:
 *   gradient[i] = [dE/dx_i, dE/dy_i, dE/dz_i]
 * with units of kcal/mol/Å.
 *
 * The gradient is the negative of the force on each atom:
 *   F_i = −∇_i E
 *
 * Each term's gradient lives in its own file (mirroring the energy/
 * layout) and is validated against finite differences in
 * tests/gradient.test.ts. The 1-4 electrostatic scaling (×0.75) is
 * applied inside the electrostatic term's gradient, exactly as the
 * energy term applies it — the term functions are the only place
 * that scaling can live, since they return full arrays, not pairs.
 */
import type { TypedMolecule } from '../../types';
import { calc_bond_stretch_gradient } from './bond-stretch';
import { calc_angle_bend_gradient } from './angle-bend';
import { calc_stretch_bend_gradient } from './stretch-bend';
import { calc_torsion_gradient } from './torsion';
import { calc_vdw_gradient } from './van-der-waals';
import { calc_electrostatic_gradient } from './electrostatic';
import { calc_oop_gradient } from './out-of-plane';

/**
 * Compute the full gradient of the MMFF94 energy.
 *
 * Returns an array parallel to molecule.atoms[]:
 *   result[i] = [dE/dx_i, dE/dy_i, dE/dz_i]
 * with units of kcal/mol/Å.
 */
export function calc_gradient(molecule: TypedMolecule): number[][] {
  const gradient: number[][] = molecule.atoms.map(() => [0, 0, 0]);

  const terms = [
    calc_bond_stretch_gradient(molecule),
    calc_angle_bend_gradient(molecule),
    calc_stretch_bend_gradient(molecule),
    calc_torsion_gradient(molecule),
    calc_vdw_gradient(molecule),
    calc_electrostatic_gradient(molecule),
    calc_oop_gradient(molecule),
  ];

  for (const term of terms) {
    for (let a = 0; a < molecule.atoms.length; a++) {
      gradient[a][0] += term[a][0];
      gradient[a][1] += term[a][1];
      gradient[a][2] += term[a][2];
    }
  }

  return gradient;
}
