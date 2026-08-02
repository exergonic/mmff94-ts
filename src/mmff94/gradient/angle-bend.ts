/**
 * Gradient of the angle bending energy.
 *
 * See energy/angle-bend.ts for the energy — two forms:
 *
 *   Regular (eq. 3):  E = 0.043844 · (k_a/2) · Δθ² · (1 + cb·Δθ)
 *     with Δθ = θ_deg − θ₀ (degrees) and cb = −0.007 deg⁻¹:
 *     dE/dθ_deg = 0.043844 · k_a · Δθ · (1 + cb·Δθ)
 *               + 0.043844 · (k_a/2) · Δθ² · cb
 *
 *   Linear centers (eq. 4):  E = 143.9325 · k_a · (1 + cos θ)
 *     dE/dθ = −143.9325 · k_a · sin θ
 *
 * The energy uses θ in DEGREES for eq. (3), so dE/dθ_deg must be
 * converted by (180/π) before multiplying with the geometric
 * dθ_rad/dx from derivatives.ts (radians per Ångström).
 *
 * Parameter resolution and angle enumeration are identical to the
 * energy term: every pair of neighbors (i, k) of every central atom j,
 * with the same class-aware angle_parameters() lookup — including the
 * lin flag that selects the cosine form.
 */

import type { TypedMolecule } from '../../types';
import { Vec3, angle_in_radians } from '../../utils/vector';
import { make_class_context, angle_parameters } from '../parameters/parameter-classes';
import { angle_derivatives, RAD_PER_DEG } from './derivatives';

const ANGLE_UNIT = 0.043844; // (mdyn·Å/rad²) → (kcal/mol)/deg²
const LINEAR_UNIT = 143.9325; // eq. (4) unit — same as bond stretch
const CB = -0.007; // cubic bend constant, deg⁻¹

/**
 * Gradient of the angle bending energy, dE/dx per atom.
 * Returns an array parallel to molecule.atoms[].
 */
export function calc_angle_bend_gradient(molecule: TypedMolecule): number[][] {
  const gradient: number[][] = molecule.atoms.map(() => [0, 0, 0]);

  // Build adjacency list — same as the energy term
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }
  const ctx = make_class_context(molecule, adj);

  for (let j = 0; j < molecule.atoms.length; j++) {
    const neighbors = adj[j];
    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];

    for (let i_idx = 0; i_idx < neighbors.length; i_idx++) {
      for (let k_idx = i_idx + 1; k_idx < neighbors.length; k_idx++) {
        const i = neighbors[i_idx];
        const k = neighbors[k_idx];

        const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];
        const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

        const { k_a, theta0, linear } = angle_parameters(ctx, i, j, k);

        // dE/dθ (kcal/mol per radian of the angle)
        let dE_dtheta: number;
        if (linear) {
          // eq. (4) derivative; θ in radians here
          const theta_rad = angle_in_radians(posJ, posI, posK);
          dE_dtheta = -LINEAR_UNIT * k_a * Math.sin(theta_rad);
        } else {
          const theta_deg = angle_in_radians(posJ, posI, posK) / RAD_PER_DEG;
          const delta_theta = theta_deg - theta0;
          dE_dtheta =
            (ANGLE_UNIT * k_a * delta_theta * (1.0 + CB * delta_theta) +
              ANGLE_UNIT * (k_a / 2.0) * delta_theta * delta_theta * CB) /
            RAD_PER_DEG; // dθ_deg → dθ_rad
        }

        const { d_dx_i, d_dx_j, d_dx_k } = angle_derivatives(posJ, posI, posK);
        for (let a = 0; a < 3; a++) {
          gradient[i][a] += dE_dtheta * d_dx_i[a];
          gradient[j][a] += dE_dtheta * d_dx_j[a];
          gradient[k][a] += dE_dtheta * d_dx_k[a];
        }
      }
    }
  }

  return gradient;
}
