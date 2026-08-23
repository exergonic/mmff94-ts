/**
 * Gradient of the torsional (dihedral) energy.
 *
 * See energy/torsion.ts for the energy — Halgren1996 eq. (7):
 *
 *   E_tors = ½·V₁·(1 + cos τ) + ½·V₂·(1 − cos 2τ) + ½·V₃·(1 + cos 3τ)
 *
 * whose derivative with respect to the dihedral angle is
 *
 *   dE/dτ = −½·V₁·sin τ + V₂·sin 2τ − 1.5·V₃·sin 3τ
 *
 * The geometric factor dτ/dx comes from derivatives.ts, which
 * differentiates the same normalized-cross-product construction as
 * dihedral_angle() — same sign convention, same handedness, so the
 * gradient agrees with the energy's own angle definition by
 * construction.
 *
 * Dihedral enumeration and parameter resolution are identical to the
 * energy term: every i on j and l on k around every central bond,
 * with the shared torsion_terms() helper (TTijkl class → step-down →
 * part-IV empirical rules).
 */

import type { TypedMolecule } from '../../types.js';
import { Vec3, dihedral_angle } from '../../utils/vector.js';
import { class_context_for } from '../parameters/parameter-classes.js';
import { torsion_terms } from '../energy/torsion.js';
import { dihedral_derivatives } from './derivatives.js';

/**
 * Gradient of the torsion energy, dE/dx per atom.
 * Returns an array parallel to molecule.atoms[].
 */
export function calc_torsion_gradient(molecule: TypedMolecule): number[][] {
  const gradient: number[][] = molecule.atoms.map(() => [0, 0, 0]);

  // Build adjacency list — same as the energy term
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }
  const ctx = class_context_for(molecule, adj);

  // Each bond is the central j-k of its dihedrals i-j-k-l
  for (const bond of molecule.bonds) {
    const j = bond.atom1;
    const k = bond.atom2;

    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];
    const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

    // Substituents of j (excluding k) and of k (excluding j)
    const i_neighbors = adj[j].filter(n => n !== k);
    const l_neighbors = adj[k].filter(n => n !== j);

    if (i_neighbors.length === 0 || l_neighbors.length === 0) continue;

    for (const i of i_neighbors) {
      const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];

      for (const l of l_neighbors) {
        // A 3-membered ring closes on itself: i and l are the same
        // atom — not a dihedral, and the energy term skips it too
        // (BatchMin's d == a rule; see energy/torsion.ts).
        if (l === i) continue;

        const posL: Vec3 = [molecule.atoms[l].x, molecule.atoms[l].y, molecule.atoms[l].z];

        const terms = torsion_terms(ctx, molecule, i, j, k, l);
        if (!terms) continue;

        const tau_rad = dihedral_angle(posI, posJ, posK, posL);

        // dE/dτ — the derivative of the three-term Fourier series
        const dE_dtau =
          -0.5 * terms.v1 * Math.sin(tau_rad) +
          terms.v2 * Math.sin(2.0 * tau_rad) -
          1.5 * terms.v3 * Math.sin(3.0 * tau_rad);

        const { d_dx_i, d_dx_j, d_dx_k, d_dx_l } = dihedral_derivatives(posI, posJ, posK, posL);
        for (let a = 0; a < 3; a++) {
          gradient[i][a] += dE_dtau * d_dx_i[a];
          gradient[j][a] += dE_dtau * d_dx_j[a];
          gradient[k][a] += dE_dtau * d_dx_k[a];
          gradient[l][a] += dE_dtau * d_dx_l[a];
        }
      }
    }
  }

  return gradient;
}
