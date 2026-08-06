/**
 * Van der Waals (non-bonded) energy.
 *
 * MMFF94 uses a "buffered 14-7" potential — NOT the Lennard-Jones 12-6
 * used by UFF and GAFF. The form, shown in Halgren1996, eq. (8), is:
 *
 *   E_vdw = ε_ij · [ (1.07·R* / (r + 0.07·R*))⁷ ·
 *                    (1.12·R*⁷ / (r⁷ + 0.12·R*⁷) − 2) ]
 *
 * The buffer terms (0.07·R* and 0.12·R*⁷) eliminate the singularity at
 * r = 0 that plagues the standard Lennard-Jones potential, giving a finite
 * repulsive wall. The name "14-7" refers to what the exponents would be
 * if the buffer constants were set to zero.
 *
 * At r = R*, the expression simplifies to E = −ε — the well depth.
 *
 * Per-atom reduced radius, eq. (9):
 *   R_i = A_i · α_i^0.25
 *
 * Waldman-Hagler combination rules, eqs. (10-11):
 *   g = (R_i − R_j) / (R_i + R_j)
 *   R_ij = 0.5 · (R_i + R_j) · [1 + 0.2 · (1 − exp(−12 · g²))]
 *   (B = 0.2 and Beta = 12 are the global vdW parameters from the
 *    mmffvdw.par file header)
 *
 * Slater-Kirkwood well depth, eq. (12):
 *   ε_ij  = 181.16 · G_i · G_j · α_i · α_j /
 *           [√(α_i / N_i) + √(α_j / N_j)] / R_ij⁶
 *
 * The 1/R_ij⁶ damping is part of the MMFF94 vdW model. The Waldman-Hagler
 * combination accounts for the different sizes of heteronuclear pairs.
 *
 * For hydrogen-bond donor-acceptor pairs (DA flags), R_ij is a simple
 * arithmetic mean and is further scaled by 0.8; ε is halved.
 *
 * 1-4 SCALING: none. Halgren 1996 (p. 496): "1,4-vdW interactions are
 * not differentially scaled in MMFF94" — unlike MM2/MM3/GAFF, no 0.5
 * factor is applied at three-bond separation. (Only the electrostatic
 * term carries a 1-4 factor, 0.75, applied in total.ts.)
 */

import type { TypedMolecule } from '../../types.js';
import { VDW_PARAMS } from '../parameters/index.js';
import { distance, Vec3 } from '../../utils/vector.js';

/**
 * The pair parameters for one vdW interaction: combined radius R*
 * and well depth ε. The combination rules are the chemistry — the
 * Waldman-Hagler equations, the arithmetic-mean shortcut whenever a
 * hydrogen-bond donor is involved, and the ×0.8/halving for
 * donor-acceptor pairs.
 *
 * Shared with the gradient so the energy and its derivative always
 * combine the same pair.
 */
export interface VdwPairParams {
  R_ij: number;
  epsilon_ij: number;
}

/**
 * Per-atom vdW parameters for atom i, with the sulfinate-S=O bridge.
 *
 * The S=O oxygen of an anionic sulfinate (bonded to a type-73 sulfur)
 * is typed 7 by the reference typing rules but keyed 32 in every
 * parameter table — the same story as the angle bridge in angle.ts —
 * and the reference energies were computed with the 32 typing. The
 * 7-typed oxygen therefore uses type 32's vdW parameters (JALSOE,
 * SO18A reproduce BatchMin exactly with this mapping).
 *
 * Shared with the gradient so the energy and its derivative use the
 * same parameters.
 */
export function vdw_parameters_for(
  molecule: TypedMolecule,
  adj: number[][],
  i: number,
): { A_i: number; alpha_i: number; N_i: number; G_i: number; DA: number } | undefined {
  const t = molecule.atom_types[i];
  if (t === 7 && adj[i].some(nb => molecule.atom_types[nb] === 73)) {
    return VDW_PARAMS[32];
  }
  // Metal-hydrate oxidation states: the original program typed the
  // +2/+1 cations (FE+2, CU+1) with their own vdW rows (87/97, the
  // paper's values); OpenBabel's canonical typing collapses them onto
  // the +3/+2 rows (88/98). The .mmd's formal charge carries the
  // oxidation state — bridge to the +2/+1 rows so the hydrates use
  // the reference's own parameters (FE2PW3/CU1PW1; same pattern as
  // the sulfinate bridge above).
  if (t === 88 && molecule.atoms[i].formal_charge === 2) return VDW_PARAMS[87];
  if (t === 98 && molecule.atoms[i].formal_charge === 1) return VDW_PARAMS[97];
  return VDW_PARAMS[t];
}

/**
 * Combine the per-atom vdW parameters of i and j (see the header of
 * this file for the equations). The VDW_PARAMS entries carry the
 * per-atom reduced radius pieces (A_i, α_i, N_i, G_i, DA flag).
 */
export function vdw_pair_parameters(
  param_i: { A_i: number; alpha_i: number; N_i: number; G_i: number; DA: number },
  param_j: { A_i: number; alpha_i: number; N_i: number; G_i: number; DA: number },
): VdwPairParams {
  // Per-atom reduced radius (eq. 9) and the Slater-Kirkwood pieces
  const R_i = param_i.A_i * Math.pow(param_i.alpha_i, 0.25);
  const R_j = param_j.A_i * Math.pow(param_j.alpha_i, 0.25);
  const sqrt_alpha_over_N_i = Math.sqrt(param_i.alpha_i / param_i.N_i);
  const sqrt_alpha_over_N_j = Math.sqrt(param_j.alpha_i / param_j.N_i);

  const is_donor = param_i.DA === 1 || param_j.DA === 1;    // DA flag: 1 = H-bond donor
  const is_acceptor = param_i.DA === 2 || param_j.DA === 2; //          2 = H-bond acceptor

  // Combined radius and well depth. The arithmetic-mean combination
  // applies whenever EITHER atom is a donor (MMFF part V) — the
  // Waldman-Hagler rules are only for pairs with no donor at all
  // (e.g. C...C, O...O). ε is halved and R* scaled by 0.8 only for
  // donor-acceptor pairs; ε always uses the UNSCALED R*.
  let R_ij: number;
  let epsilon_ij: number;

  if (is_donor) {
    const R_ij_unscaled = 0.5 * (R_i + R_j);
    const R_ij6_unscaled = Math.pow(R_ij_unscaled, 6);
    if (is_acceptor) {
      // Hydrogen bond donor-acceptor: epsilon halved, R* × 0.8
      epsilon_ij = 0.5 * (181.16 * param_i.G_i * param_j.G_i * param_i.alpha_i * param_j.alpha_i) /
                   (sqrt_alpha_over_N_i + sqrt_alpha_over_N_j) / R_ij6_unscaled;
      R_ij = 0.8 * R_ij_unscaled;
    } else {
      // Donor present but no acceptor (e.g. H-N...H-C): arithmetic
      // mean, full epsilon
      epsilon_ij = (181.16 * param_i.G_i * param_j.G_i * param_i.alpha_i * param_j.alpha_i) /
                   (sqrt_alpha_over_N_i + sqrt_alpha_over_N_j) / R_ij6_unscaled;
      R_ij = R_ij_unscaled;
    }
  } else {
    // No donor: Waldman-Hagler combination
    const asymmetry = (R_i - R_j) / (R_i + R_j);
    const asymmetry2 = asymmetry * asymmetry;
    R_ij = 0.5 * (R_i + R_j) * (1.0 + 0.2 * (1.0 - Math.exp(-12.0 * asymmetry2)));
    const R_ij6 = Math.pow(R_ij, 6);
    epsilon_ij = (181.16 * param_i.G_i * param_j.G_i * param_i.alpha_i * param_j.alpha_i) /
                 (sqrt_alpha_over_N_i + sqrt_alpha_over_N_j) / R_ij6;
  }

  return { R_ij, epsilon_ij };
}

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
    const param_i = vdw_parameters_for(molecule, adj, i);
    if (!param_i) continue;

    const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];

    for (let j = i + 1; j < molecule.atoms.length; j++) {

      // Skip 1-2 (bonded) and 1-3 (share a common neighbor)
      if (adj[i].includes(j)) continue;
      if (pairs_1_3[i].has(j)) continue;

      const param_j = vdw_parameters_for(molecule, adj, j);
      if (!param_j) continue;

      const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];
      const r = distance(posI, posJ);

      // Combined radius and well depth — the combination rules live
      // in vdw_pair_parameters(), shared with the gradient.
      const { R_ij, epsilon_ij } = vdw_pair_parameters(param_i, param_j);

      // Buffered 14-7 expression, Halgren1996 eq. (8)
      const r7 = Math.pow(r, 7);
      const R_ij7 = Math.pow(R_ij, 7);

      const repulsive_term = Math.pow(1.07 * R_ij / (r + 0.07 * R_ij), 7);
      const attractive_term = 1.12 * R_ij7 / (r7 + 0.12 * R_ij7) - 2;

      total_energy += epsilon_ij * repulsive_term * attractive_term;
    }
  }

  return total_energy;
}
