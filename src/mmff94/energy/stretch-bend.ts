/**
 * Stretch-bend cross term energy.
 *
 * Halgren1996, eq. (5):
 *
 *   E_sb = 2.51210 · [k_sb_IJK · (r_IJ − r_IJ0)
 *                     + k_sb_KJI · (r_KJ − r_KJ0)]
 *                    · (θ − θ₀)
 *
 * This is a CLASS II force field term — it couples bond stretching with
 * angle bending. Most force fields (UFF, GAFF, MM2/3) omit this term, but
 * MMFF94 includes it because bond lengths and angles are physically coupled:
 * when an H-C-H angle in methane closes, the C-H bonds shorten slightly.
 *
 * Two separate k_sb values are used because asymmetric environments
 * (e.g., C-C-O vs O-C-C) have different coupling strengths for the two
 * sides of the angle. Many entries have k_sb_IJK = k_sb_KJI (symmetric
 * angles), but the parameter table stores both independently.
 */

import type { TypedMolecule } from '../../types';
import {
  STRETCH_BEND_PARAMS,
  DEFAULT_STRETCH_BEND,
  ELEMENT_ROW,
  lookup_param,
  type StretchBendParams,
} from '../parameters';
import { distance, angle_in_radians, Vec3 } from '../../utils/vector';
import { make_class_context, type ClassContext, strbnd_type, bond_parameters, angle_parameters } from '../parameters/parameter-classes';

/**
 * All parameters a stretch-bend angle i-j-k needs: the two coupling
 * constants (k_ij for the i-j bond side, k_kj for the k-j side), the
 * two equilibrium bond lengths, the equilibrium angle, and the linear
 * flag.
 *
 * Shared by the energy term and its gradient so both resolve the
 * same parameters by construction — the resolution chain (strbnd
 * class → table → element-row defaults) is the part that carries
 * the chemistry, and it must not drift between the two.
 */
export interface StretchBendAngleTerms {
  k_ij: number;
  k_kj: number;
  r0_ij: number;
  r0_kj: number;
  theta0: number;
  linear: boolean;
}

/**
 * Resolve the stretch-bend parameters for the angle i-j-k.
 * Returns undefined when no parameters exist at all (no table entry,
 * no element-row default, or a missing bond parameter).
 */
export function stretch_bend_angle_terms(
  ctx: ClassContext,
  molecule: TypedMolecule,
  i: number,
  j: number,
  k: number,
): StretchBendAngleTerms | undefined {
  const tj = molecule.atom_types[j];
  const ti = molecule.atom_types[i];
  const tk = molecule.atom_types[k];

  const t_min = Math.min(ti, tk);
  const t_max = Math.max(ti, tk);

  // The stretch-bend class is a remap of the angle class
  // (GetStrBndType): for BT-flagged angles the class splits 1/2 by
  // which side carries the flag, so the ring/BT classes resolve to
  // different keys than the angle term's.
  const cls = strbnd_type(ctx, i, j, k);

  // 1. Look up stretch-bend parameters
  let sb_params: StretchBendParams | undefined;
  if (cls !== 0) {
    const keys = [
      `${cls}-${t_min}-${tj}-${t_max}`,
      `${cls}-0-${tj}-${t_max}`,
      `${cls}-${t_min}-${tj}-0`,
      `${cls}-0-${tj}-0`,
    ];
    for (const key of keys) {
      sb_params = STRETCH_BEND_PARAMS[key];
      if (sb_params) break;
    }
  } else {
    sb_params = lookup_param(STRETCH_BEND_PARAMS, [t_min, tj, t_max]);
  }
  // No stretch-bend entry: use the default F(I_J,K)/F(K_J,I)
  // from the element-row table (mmffdfsb.par) — BatchMin
  // evaluates every angle, and the defaults are the only values
  // for e.g. Si angles (a stretched Si–C bond can carry several
  // kcal/mol here, so skipping the angle outright is wrong).
  let k_ij: number;
  let k_kj: number;
  if (sb_params) {
    // The table stores terminal types sorted (I ≤ K), so k_sb_IJK
    // belongs to the bond on the min-type side and k_sb_KJI to the
    // bond on the max-type side — whichever of i and k that is.
    k_ij = ti <= tk ? sb_params.k_sb_IJK : sb_params.k_sb_KJI;
    k_kj = ti <= tk ? sb_params.k_sb_KJI : sb_params.k_sb_IJK;
  } else {
    const rowa = ELEMENT_ROW[molecule.atoms[i].element] ?? 0;
    const rowb = ELEMENT_ROW[molecule.atoms[j].element] ?? 0;
    const rowc = ELEMENT_ROW[molecule.atoms[k].element] ?? 0;
    const direct = DEFAULT_STRETCH_BEND[`${rowa}-${rowb}-${rowc}`];
    const F = direct ?? DEFAULT_STRETCH_BEND[`${rowc}-${rowb}-${rowa}`];
    if (!F) return undefined;
    // F[0] belongs to the I-J bond when the stored row order is
    // (rowa, rowb, rowc); for the reversed match it belongs to K-J.
    k_ij = direct ? F[0] : F[1];
    k_kj = direct ? F[1] : F[0];
  }

  // 2. Look up equilibrium bond lengths from bond parameters —
  //    each bond uses its own type pair (sorted), not the angle's
  //    sorted terminal types.
  const bond_ij = bond_parameters(ctx, i, j);
  const bond_kj = bond_parameters(ctx, j, k);
  if (!bond_ij || !bond_kj) return undefined;

  // 3. Equilibrium angle from the same class-aware resolution the
  //    angle term uses, so ring and BT-flagged angles share θ₀.
  const { theta0, linear } = angle_parameters(ctx, i, j, k);

  return { k_ij, k_kj, r0_ij: bond_ij.r0, r0_kj: bond_kj.r0, theta0, linear };
}

/**
 * Calculate the total stretch-bend cross term energy.
 */
export function calc_stretch_bend_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  // Build adjacency list
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }
  const ctx = make_class_context(molecule, adj);

  // Iterate over all possible central atoms j
  for (let j = 0; j < molecule.atoms.length; j++) {
    const neighbors = adj[j];
    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];

    // For each pair of neighbors (i, k) — these define angles i-j-k
    for (let i_idx = 0; i_idx < neighbors.length; i_idx++) {
      for (let k_idx = i_idx + 1; k_idx < neighbors.length; k_idx++) {
        const i = neighbors[i_idx];
        const k = neighbors[k_idx];

        // All parameter resolution lives in stretch_bend_angle_terms —
        // shared with the gradient so the two can never drift apart.
        const terms = stretch_bend_angle_terms(ctx, molecule, i, j, k);
        if (!terms || terms.linear) continue;

        // 4. Compute current geometry
        const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];
        const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

        const r_IJ = distance(posI, posJ);
        const r_KJ = distance(posK, posJ);
        const theta_rad = angle_in_radians(posJ, posI, posK);
        const theta_deg = theta_rad * (180.0 / Math.PI);

        // 5. Accumulate energy
        const dr_IJ = r_IJ - terms.r0_ij;
        const dr_KJ = r_KJ - terms.r0_kj;
        const d_theta = theta_deg - terms.theta0;

        total_energy += 2.51210 * (terms.k_ij * dr_IJ + terms.k_kj * dr_KJ) * d_theta;
      }
    }
  }

  return total_energy;
}
