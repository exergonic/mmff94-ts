/**
 * Van der Waals (non-bonded) energy.
 *
 * MMFF94 uses a "buffered 14-7" potential — NOT the Lennard-Jones 12-6
 * used by UFF and GAFF. The form is:
 *
 *   E_vdw = ε_ij · [ (1.07·R* / (r + 0.07·R*))⁷ ·
 *                    (1.12·R*⁷ / (r⁷ + 0.12·R*⁷) − 2) ]
 *
 * where R* and ε are derived from the per-atom parameters A_i, α_i, N_i,
 * G_i, and DA via Waldman-Hagler combination rules with 1/r⁶ damping.
 *
 * The buffer terms (0.07·R* and 0.12·R*⁷) eliminate the singularity at
 * r = 0 that plagues the standard Lennard-Jones potential, giving a finite
 * repulsive wall.
 *
 * At r = R*, the expression simplifies to E = −ε — the well depth.
 *
 * PER-ATOM REDUCED RADIUS:
 *   R_i = A_i · α_i^0.25
 *
 * COMBINATION RULES (non-hydrogen-bond pairs):
 *   g = (R_i − R_j) / (R_i + R_j)
 *   R_ij = 0.5 · (R_i + R_j) · [1 + 0.2 · (1 − exp(−12 · g²))]
 *   ε_ij  = 181.16 · G_i · G_j · α_i · α_j /
 *           [α_i / √(N_i) + α_j / √(N_j)] / R_ij⁶
 *
 * The Waldman-Hagler combination (with B = 0.2, Beta = 12) accounts for
 * the different sizes of heteronuclear pairs. The 1/R_ij⁶ damping of the
 * well depth is part of the MMFF94 vdW model.
 *
 * For hydrogen-bond donor-acceptor pairs, R_ij is a simple arithmetic
 * mean and is further scaled by 0.8; ε is halved.
 *
 * 1-4 SCALING: vdW is multiplied by 0.5 for atoms exactly three bonds
 * apart. This is applied in total.ts, not here.
 */

import type { TypedMolecule } from '../../types';
import { VDW_PARAMS } from '../parameters';
import { distance, Vec3 } from '../../utils/vector';

/**
 * Calculate the total van der Waals energy between all non-bonded atom pairs.
 *
 * Excludes 1-2 (bonded) and 1-3 (angle) pairs. 1-4 pairs are included
 * at full strength (0.5 scaling applied externally in total.ts).
 */
export function calc_vdw_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  // Build adjacency and 2-away sets for pair exclusions
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  // Build 1-3 pair map: atoms sharing a common neighbor (angle pairs)
  const pairs_1_3: Set<number>[] = Array.from({ length: molecule.atoms.length }, () => new Set());
  for (let i = 0; i < molecule.atoms.length; i++) {
    for (const n1 of adj[i]) {
      for (const n2 of adj[n1]) {
        if (n2 !== i) pairs_1_3[i].add(n2);
      }
    }
  }

  for (let i = 0; i < molecule.atoms.length; i++) {
    const param_i = VDW_PARAMS[molecule.atom_types[i]];
    if (!param_i) continue;

    // Per-atom reduced radius
    const R_i = param_i.A_i * Math.pow(param_i.alpha_i, 0.25);
    const sqrt_alpha_over_N_i = Math.sqrt(param_i.alpha_i / param_i.N_i);

    const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];

    for (let j = i + 1; j < molecule.atoms.length; j++) {

      // Skip 1-2 (bonded) and 1-3 (share a common neighbor)
      if (adj[i].includes(j)) continue;
      if (pairs_1_3[i].has(j)) continue;

      const param_j = VDW_PARAMS[molecule.atom_types[j]];
      if (!param_j) continue;

      const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];
      const r = distance(posI, posJ);

      // Per-atom reduced radius for j
      const R_j = param_j.A_i * Math.pow(param_j.alpha_i, 0.25);
      const sqrt_alpha_over_N_j = Math.sqrt(param_j.alpha_i / param_j.N_i);

      const isDonor = param_i.DA === 1 || param_j.DA === 1;
      const isAcceptor = param_i.DA === 2 || param_j.DA === 2;

      // Combined radius and well depth
      let R_ij: number;
      let epsilon_ij: number;

      if (isDonor && isAcceptor) {
        // Hydrogen bond donor-acceptor: arithmetic mean, epsilon halved,
        // R_ij further scaled by 0.8
        const R_ij_unscaled = 0.5 * (R_i + R_j);
        const R_ij6_unscaled = Math.pow(R_ij_unscaled, 6);
        epsilon_ij = 0.5 * (181.16 * param_i.G_i * param_j.G_i * param_i.alpha_i * param_j.alpha_i) /
                     (sqrt_alpha_over_N_i + sqrt_alpha_over_N_j) / R_ij6_unscaled;
        R_ij = 0.8 * R_ij_unscaled;
      } else {
        // Non-hydrogen-bond: Waldman-Hagler combination
        const asymmetry = (R_i - R_j) / (R_i + R_j);
        const asymmetry2 = asymmetry * asymmetry;
        R_ij = 0.5 * (R_i + R_j) * (1.0 + 0.2 * (1.0 - Math.exp(-12.0 * asymmetry2)));
        const R_ij6 = Math.pow(R_ij, 6);
        epsilon_ij = (181.16 * param_i.G_i * param_j.G_i * param_i.alpha_i * param_j.alpha_i) /
                     (sqrt_alpha_over_N_i + sqrt_alpha_over_N_j) / R_ij6;
      }

      // Buffered 14-7 expression
      const R_ij6 = Math.pow(R_ij, 6);
      const R_ij7 = R_ij6 * R_ij;
      const buff_r_plus = r + 0.07 * R_ij;
      const buff_r7_plus = Math.pow(r, 7) + 0.12 * R_ij7;

      const repulsive_term = Math.pow(1.07 * R_ij / buff_r_plus, 7);
      const attractive_term = 1.12 * R_ij7 / buff_r7_plus - 2;

      total_energy += epsilon_ij * repulsive_term * attractive_term;
    }
  }

  return total_energy;
}
