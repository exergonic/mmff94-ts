/**
 * Gradient of the out-of-plane bending energy.
 *
 * See energy/out-of-plane.ts for the energy — Halgren1996 eq. (6):
 *
 *   E_oop = OOP_UNIT · (k_oop/2) · χ²
 *
 * with χ in degrees, so for each of the three Wilson angles at a
 * tri-coordinate center j:
 *
 *   dE/dx = OOP_UNIT · k_oop · χ_deg · dχ_deg/dx
 *
 * The geometric factor dχ/dx comes from derivatives.ts, which
 * differentiates the same unit-normal construction as
 * wilson_oop_angle() (unit vectors (i−j), (k−j), (l−j); normal
 * û_ji × û_jk; sin χ = n̂·û_jl) — same sign convention, so the
 * gradient agrees with the energy's own χ by construction. The
 * center's k_oop comes from the shared oop_force_constant() (sorted
 * substituent multiset, then the per-central-type wildcard).
 *
 * Note the k_oop sign carries through: amine N's explicit zero
 * contributes nothing, and amide N's NEGATIVE constant correctly
 * produces forces that push the nitrogen toward its pyramidal
 * equilibrium — the gradient of χ² is 2χ·χ′, and the negative k
 * flips the sign of the whole term.
 */

import type { TypedMolecule } from '../../types.js';
import { Vec3, wilson_oop_angle } from '../../utils/vector.js';
import { oop_force_constant } from '../energy/out-of-plane.js';
import { oop_angle_derivatives, RAD_PER_DEG } from './derivatives.js';

const OOP_UNIT = 143.9325 * (Math.PI / 180) ** 2; // exact form of the published 0.043844 (see angle-bend.ts)

/**
 * Gradient of the out-of-plane bending energy, dE/dx per atom.
 * Returns an array parallel to molecule.atoms[].
 */
export function calc_oop_gradient(molecule: TypedMolecule): number[][] {
  const gradient: number[][] = molecule.atoms.map(() => [0, 0, 0]);

  // Build adjacency list — same as the energy term
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  for (let j = 0; j < molecule.atoms.length; j++) {
    const neighbors = adj[j];
    if (neighbors.length !== 3) continue;

    const [a, c, d] = neighbors;
    const k_oop = oop_force_constant(molecule, j, a, c, d);
    if (k_oop === undefined) continue;

    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];
    const posA: Vec3 = [molecule.atoms[a].x, molecule.atoms[a].y, molecule.atoms[a].z];
    const posC: Vec3 = [molecule.atoms[c].x, molecule.atoms[c].y, molecule.atoms[c].z];
    const posD: Vec3 = [molecule.atoms[d].x, molecule.atoms[d].y, molecule.atoms[d].z];

    // The three Wilson angles at j — the same turns as the energy
    // term: each substituent takes a turn as the out-of-plane atom,
    // with the plane through j and the other two.
    const angles: { chi: number; deriv: ReturnType<typeof oop_angle_derivatives>; atoms: [number, number, number, number] }[] = [
      {
        chi: wilson_oop_angle(posD, posJ, posC, posA),
        deriv: oop_angle_derivatives(posD, posJ, posC, posA),
        atoms: [d, j, c, a], // i, j, k, l in the helper's convention
      },
      {
        chi: wilson_oop_angle(posA, posJ, posD, posC),
        deriv: oop_angle_derivatives(posA, posJ, posD, posC),
        atoms: [a, j, d, c],
      },
      {
        chi: wilson_oop_angle(posA, posJ, posC, posD),
        deriv: oop_angle_derivatives(posA, posJ, posC, posD),
        atoms: [a, j, c, d],
      },
    ];

    for (const { chi, deriv, atoms } of angles) {
      // dE/dχ_deg (kcal/mol per degree), then the chain rule:
      // dE/dx = dE/dχ_deg · dχ_deg/dx, with dχ_deg/dx = dχ_rad/dx / (π/180)
      const dE_dchi = OOP_UNIT * k_oop * chi;
      const [atom_i, atom_j, atom_k, atom_l] = atoms;
      const contribs = [
        [atom_i, deriv.d_dx_i],
        [atom_j, deriv.d_dx_j],
        [atom_k, deriv.d_dx_k],
        [atom_l, deriv.d_dx_l],
      ] as const;
      for (const [atom, d_dx] of contribs) {
        for (let axis = 0; axis < 3; axis++) {
          gradient[atom][axis] += dE_dchi * d_dx[axis] / RAD_PER_DEG;
        }
      }
    }
  }

  return gradient;
}
