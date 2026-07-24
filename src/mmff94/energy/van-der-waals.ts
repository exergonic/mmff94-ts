/**
 * Van der Waals (non-bonded) energy.
 *
 * MMFF94 uses a "buffered 14-7" potential — NOT the Lennard-Jones 12-6
 * used by UFF and GAFF. The form is:
 *
 *   E_vdw = ε_ij · [ (1.07·R* / (r + 0.07·R*))⁷ ·
 *                    (1.12·R*⁷ / (r⁷ + 0.12·R*⁷) − 2) ]
 *
 * The buffer terms (0.07·R* and 0.12·R*⁷) eliminate the singularity at
 * r = 0 that plagues the standard Lennard-Jones potential, giving a finite
 * repulsive wall. This makes MMFF94 more numerically stable during
 * optimization when atoms may approach very closely.
 *
 * At r = R* (the equilibrium distance), the expression simplifies to
 * E = −ε — the well depth — because both buffer fractions become exactly 1.
 *
 * COMBINATION RULES:
 *   R*_ij = 0.5 · (R*_i + R*_j)                     (arithmetic mean)
 *   ε_ij  = 181.16 · G_i · G_j · α_i · α_j /
 *           [ α_i / √(N_i) + α_j / √(N_j) ]         (Slater-Kirkwood)
 *
 * Slater-Kirkwood is used instead of the geometric-mean combination
 * found in simpler force fields because geometric-mean systematically
 * overestimates well depths for heteronuclear pairs (Halgren 1996, §III.C).
 *
 * 1-4 SCALING: For atoms exactly three bonds apart, the vdW energy
 * is multiplied by 0.5. This is applied in total.ts, not here.
 * This function computes the FULL unscaled energy for every pair.
 */

import type { TypedMolecule } from '../../types';
import { VDW_PARAMS } from '../parameters';
import { distance, Vec3 } from '../../utils/vector';

/**
 * Calculate the total van der Waals energy between all non-bonded atom pairs.
 *
 * Excludes 1-2 (bonded) and 1-3 (angle) pairs. 1-4 pairs are included
 * at full strength (scaling applied externally in total.ts).
 *
 * Direct O(n²) pair loop — fine for molecules under ~100 atoms.
 * Swap to a neighbor-list when benchmarking shows the need.
 */
export function calc_vdw_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  // Build adjacency: for each atom, which other atoms is it bonded to?
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  // For each atom i, collect all atoms j that are 1-3 (share a common neighbor)
  // so we can skip them in the pair loop.
  const adj2: Set<number>[] = Array.from({ length: molecule.atoms.length }, () => new Set());
  for (let i = 0; i < molecule.atoms.length; i++) {
    for (const n1 of adj[i]) {
      for (const n2 of adj[n1]) {
        if (n2 !== i) adj2[i].add(n2);
      }
    }
  }

  for (let i = 0; i < molecule.atoms.length; i++) {
    const ti = molecule.atom_types[i];
    const params_i = VDW_PARAMS[ti];
    if (!params_i) continue;

    const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];

    for (let j = i + 1; j < molecule.atoms.length; j++) {

      // Skip 1-2 pairs (bonded)
      if (adj[i].includes(j)) continue;

      // Skip 1-3 pairs (share a common neighbor)
      if (adj2[i].has(j)) continue;

      const tj = molecule.atom_types[j];
      const params_j = VDW_PARAMS[tj];
      if (!params_j) continue;

      const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];
      const r = distance(posI, posJ);

      // Combination rules
      const R_star_ij = 0.5 * (params_i.R_star + params_j.R_star);

      // Slater-Kirkwood combination for well depth
      const alpha_over_sqrt_N_i = params_i.alpha_i / Math.sqrt(params_i.N_i);
      const alpha_over_sqrt_N_j = params_j.alpha_i / Math.sqrt(params_j.N_i);
      const epsilon_ij = 181.16 * params_i.G_i * params_j.G_i *
                         params_i.alpha_i * params_j.alpha_i /
                         (alpha_over_sqrt_N_i + alpha_over_sqrt_N_j);

      // Buffered 14-7 expression
      const buff_r_plus = r + 0.07 * R_star_ij;
      const term1 = Math.pow(1.07 * R_star_ij / buff_r_plus, 7);

      const r7 = Math.pow(r, 7);
      const R_star_ij_7 = Math.pow(R_star_ij, 7);
      const buff_r7_plus = r7 + 0.12 * R_star_ij_7;
      const term2 = 1.12 * R_star_ij_7 / buff_r7_plus - 2;

      total_energy += epsilon_ij * term1 * term2;
    }
  }

  return total_energy;
}
