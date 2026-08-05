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
import type { Molecule } from '../../types.js';
import { prepare_molecule } from '../prepare.js';
import { calc_bond_stretch_gradient } from './bond-stretch.js';
import { calc_angle_bend_gradient } from './angle-bend.js';
import { calc_stretch_bend_gradient } from './stretch-bend.js';
import { calc_torsion_gradient } from './torsion.js';
import { calc_vdw_gradient } from './van-der-waals.js';
import { calc_electrostatic_gradient } from './electrostatic.js';
import { calc_oop_gradient } from './out-of-plane.js';

/**
 * Compute the full gradient of the MMFF94 energy.
 *
 * Accepts a bare Molecule straight from parse_sdf() — atom typing and
 * BCI charges are assigned on demand (prepare_molecule). An
 * already-prepared TypedMolecule passes through untouched.
 *
 * Returns an array parallel to molecule.atoms[]:
 *   result[i] = [dE/dx_i, dE/dy_i, dE/dz_i]
 * with units of kcal/mol/Å.
 */
export function calc_gradient(molecule: Molecule): number[][] {
  const prepared = prepare_molecule(molecule);
  const gradient: number[][] = prepared.atoms.map(() => [0, 0, 0]);

  const terms = [
    calc_bond_stretch_gradient(prepared),
    calc_angle_bend_gradient(prepared),
    calc_stretch_bend_gradient(prepared),
    calc_torsion_gradient(prepared),
    calc_vdw_gradient(prepared),
    calc_electrostatic_gradient(prepared),
    calc_oop_gradient(prepared),
  ];

  for (const term of terms) {
    for (let a = 0; a < prepared.atoms.length; a++) {
      gradient[a][0] += term[a][0];
      gradient[a][1] += term[a][1];
      gradient[a][2] += term[a][2];
    }
  }

  return gradient;
}
