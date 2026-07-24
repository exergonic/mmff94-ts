/**
 * Torsion (dihedral) energy.
 *
 * Halgren1996, eq. (7):
 *
 *   E_tors = ½ · [V₁ · (1 + cos τ) + V₂ · (1 − cos 2τ) + V₃ · (1 + cos 3τ)]
 *
 * Written more compactly using our stored phases γ_n:
 *
 *   E_tors = Σ (V_n / 2) · [1 + cos(n · τ − γ_n)]   for n = 1, 2, 3
 *
 * where:
 *   V_n    = barrier height for the n-th term (kcal/mol)
 *   n      = periodicity (1, 2, or 3)
 *   τ      = current dihedral angle (degrees)
 *   γ_n    = phase shift for the n-th term (degrees)
 *
 * The two forms are equivalent when:
 *   γ₁ = 0°    →  1 + cos τ        (n=1)
 *   γ₂ = 180°  →  1 − cos 2τ       (n=2, since cos(2τ−180°) = −cos 2τ)
 *   γ₃ = 0°    →  1 + cos 3τ       (n=3)
 *
 * Convention: τ = 0° when i−j and k−l bonds are eclipsed (cis).
 * τ = 180° when staggered (trans). Sign follows IUPAC right-hand rule.
 *
 * Only single bonds (bond_order === 1) are evaluated for torsional
 * strain. Double and triple bonds are kept planar by angle bend and
 * out-of-plane terms, not by torsion.
 *
 * The parameter lookup tries type order (i, j, k, l) first, then the
 * reverse (l, k, j, i), to catch symmetric parameter definitions.
 */

import type { TypedMolecule } from '../../types';
import { TORSION_PARAMS, lookup_param } from '../parameters';
import { dihedral_angle, Vec3 } from '../../utils/vector';

/**
 * Calculate the total torsional (dihedral) energy.
 */
export function calc_torsion_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  // Build adjacency list for neighbor lookups
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  // Iterate over all bonds; only single bonds contribute torsion
  for (const bond of molecule.bonds) {
    if (bond.bond_order !== 1) continue;

    const j = bond.atom1;
    const k = bond.atom2;

    const tj = molecule.atom_types[j];
    const tk = molecule.atom_types[k];

    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];
    const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

    // Neighbors of j (excluding k) and neighbors of k (excluding j)
    const i_neighbors = adj[j].filter(n => n !== k);
    const l_neighbors = adj[k].filter(n => n !== j);

    if (i_neighbors.length === 0 || l_neighbors.length === 0) continue;

    for (const i of i_neighbors) {
      const ti = molecule.atom_types[i];
      const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];

      for (const l of l_neighbors) {
        const tl = molecule.atom_types[l];
        const posL: Vec3 = [molecule.atoms[l].x, molecule.atoms[l].y, molecule.atoms[l].z];

        // Look up torsion parameters for ordered types (i, j, k, l)
        let params = lookup_param(TORSION_PARAMS, [ti, tj, tk, tl]);

        // If not found, try reverse order (l, k, j, i)
        if (!params) {
          params = lookup_param(TORSION_PARAMS, [tl, tk, tj, ti]);
        }

        if (!params) continue;

        // Compute dihedral angle in degrees
        const tau_rad = dihedral_angle(posI, posJ, posK, posL);
        const tau_deg = tau_rad * (180.0 / Math.PI);

        // Evaluate each Fourier term
        for (const term of params.terms) {
          const angle_rad = (term.periodicity * tau_deg - term.gamma) * (Math.PI / 180.0);
          total_energy += (term.V / 2.0) * (1.0 + Math.cos(angle_rad));
        }
      }
    }
  }

  return total_energy;
}
