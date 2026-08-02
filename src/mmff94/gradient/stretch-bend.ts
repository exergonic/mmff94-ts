/**
 * Gradient of the stretch-bend cross term energy.
 *
 * See energy/stretch-bend.ts for the energy — Halgren1996 eq. (5):
 *
 *   E_sb = 2.51210 · (k_ij·Δr_IJ + k_kj·Δr_KJ) · Δθ
 *
 * with Δr in Å, Δθ in degrees. The chain rule touches all three
 * internal coordinates:
 *
 *   dE/dx = 2.51210 · [ k_ij · (dΔr_IJ/dx) · Δθ
 *                       + k_kj · (dΔr_KJ/dx) · Δθ
 *                       + (k_ij·Δr_IJ + k_kj·Δr_KJ) · (dΔθ/dx) ]
 *
 * where dΔθ/dx is the angle derivative converted from radians to
 * degrees (the constants above are degree-based). The parameters come
 * from stretch_bend_angle_terms() — the SAME resolution the energy
 * term uses (strbnd class → table → element-row defaults), so the two
 * can never disagree about which angles exist or with which constants.
 */

import type { TypedMolecule } from '../../types';
import { Vec3, angle_in_radians } from '../../utils/vector';
import { make_class_context } from '../parameters/parameter-classes';
import {
  stretch_bend_angle_terms,
} from '../energy/stretch-bend';
import {
  bond_length_derivatives,
  angle_derivatives,
  RAD_PER_DEG,
} from './derivatives';

const SB_UNIT = 2.51210; // eq. (5) unit conversion

/**
 * Gradient of the stretch-bend energy, dE/dx per atom.
 * Returns an array parallel to molecule.atoms[].
 */
export function calc_stretch_bend_gradient(molecule: TypedMolecule): number[][] {
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

        const terms = stretch_bend_angle_terms(ctx, molecule, i, j, k);
        if (!terms || terms.linear) continue;

        const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];
        const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

        const r_IJ = Math.hypot(posI[0] - posJ[0], posI[1] - posJ[1], posI[2] - posJ[2]);
        const r_KJ = Math.hypot(posK[0] - posJ[0], posK[1] - posJ[1], posK[2] - posJ[2]);
        const theta_deg = angle_in_radians(posJ, posI, posK) / RAD_PER_DEG;

        const dr_IJ = r_IJ - terms.r0_ij;
        const dr_KJ = r_KJ - terms.r0_kj;
        const d_theta = theta_deg - terms.theta0;

        // dE/dx from the product rule — see the header. dE/dΔθ (per
        // degree of the angle) multiplies the degree-based angle
        // derivative; the k·Δr terms multiply the bond derivatives.
        const dE_ddelta_theta = SB_UNIT * (terms.k_ij * dr_IJ + terms.k_kj * dr_KJ);

        const { d_dx_a: dIJ_dx_i, d_dx_b: dIJ_dx_j } = bond_length_derivatives(posI, posJ);
        const { d_dx_a: dKJ_dx_k, d_dx_b: dKJ_dx_j } = bond_length_derivatives(posK, posJ);
        const { d_dx_i: dtheta_dx_i, d_dx_j: dtheta_dx_j, d_dx_k: dtheta_dx_k } =
          angle_derivatives(posJ, posI, posK);
        // degrees per Å
        const dtheta_deg_dx_i = dtheta_dx_i.map(v => v / RAD_PER_DEG) as Vec3;
        const dtheta_deg_dx_j = dtheta_dx_j.map(v => v / RAD_PER_DEG) as Vec3;
        const dtheta_deg_dx_k = dtheta_dx_k.map(v => v / RAD_PER_DEG) as Vec3;

        for (let a = 0; a < 3; a++) {
          // dE/dx = 2.51210·[k_ij·dΔr_IJ/dx·Δθ + k_kj·dΔr_KJ/dx·Δθ
          //                + (k_ij·Δr_IJ + k_kj·Δr_KJ)·dΔθ/dx]
          gradient[i][a] +=
            SB_UNIT * terms.k_ij * dIJ_dx_i[a] * d_theta +
            dE_ddelta_theta * dtheta_deg_dx_i[a];
          gradient[j][a] +=
            SB_UNIT * terms.k_ij * dIJ_dx_j[a] * d_theta +
            SB_UNIT * terms.k_kj * dKJ_dx_j[a] * d_theta +
            dE_ddelta_theta * dtheta_deg_dx_j[a];
          gradient[k][a] +=
            SB_UNIT * terms.k_kj * dKJ_dx_k[a] * d_theta +
            dE_ddelta_theta * dtheta_deg_dx_k[a];
        }
      }
    }
  }

  return gradient;
}
