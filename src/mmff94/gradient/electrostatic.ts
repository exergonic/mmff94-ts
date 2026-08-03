/**
 * Gradient of the electrostatic (Coulombic) energy.
 *
 * See energy/electrostatic.ts for the energy — Halgren1996 part III,
 * eq. (6):
 *
 *   E_elec = 332.0716 · q_i · q_j / (r + S)      S = 0.05 Å
 *
 * whose derivative with respect to the interatomic distance is
 *
 *   dE/dr = −332.0716 · q_i · q_j / (r + S)²
 *
 * and dE/dx = dE/dr · dr/dx along the bond direction.
 *
 * The two non-trivial consistency requirements are shared with the
 * energy term by construction:
 *
 *   - The pair list: every i < j pair excluding 1-2 and 1-3, with the
 *     1-4 pairs scaled by 0.75 — via the exported is_1_4_pair() (the
 *     same BFS the energy uses). A pair the energy scales at 0.75
 *     contributes 0.75× its unscaled gradient here.
 *   - The charges: the same BCI partial charges (computed on demand
 *     if the caller did not call assign_bci_charges()).
 */

import type { TypedMolecule } from '../../types';
import { Vec3 } from '../../utils/vector';
import { assign_bci_charges } from '../charges';
import { is_1_4_pair } from '../energy/electrostatic';
import { bond_length_derivatives } from './derivatives';

const S = 0.05; // the electrostatic buffering constant (Å)
const ELEC_UNIT = 332.0716; // e²/Å → kcal/mol
const SCALE_1_4 = 0.75; // MMFF94's 1-4 electrostatic scaling

/**
 * Gradient of the electrostatic energy, dE/dx per atom.
 * Returns an array parallel to molecule.atoms[].
 */
export function calc_electrostatic_gradient(molecule: TypedMolecule): number[][] {
  const gradient: number[][] = molecule.atoms.map(() => [0, 0, 0]);

  // Partial charges from the BCI model — same as the energy term.
  const charges =
    molecule.partial_charges ?? assign_bci_charges(molecule).partial_charges!;

  // Adjacency: for the 1-2/1-3 exclusion and the 1-4 classification.
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  // 1-3 pair map: atoms sharing a common neighbor (angle pairs)
  const pairs_1_3: Set<number>[] = Array.from({ length: molecule.atoms.length }, () => new Set());
  for (let i = 0; i < molecule.atoms.length; i++) {
    for (const n1 of adj[i]) {
      for (const n2 of adj[n1]) {
        if (n2 !== i) pairs_1_3[i].add(n2);
      }
    }
  }

  const n = molecule.atoms.length;
  const positions: Vec3[] = molecule.atoms.map(a => [a.x, a.y, a.z]);

  for (let i = 0; i < n; i++) {
    const qi = charges[i];
    if (qi === 0) continue;
    for (let j = i + 1; j < n; j++) {
      // Skip 1-2 (bonded) and 1-3 (share a common neighbor) pairs
      if (adj[i].includes(j)) continue;
      if (pairs_1_3[i].has(j)) continue;

      const qj = charges[j];
      if (qj === 0) continue;

      const r = Math.hypot(
        positions[i][0] - positions[j][0],
        positions[i][1] - positions[j][1],
        positions[i][2] - positions[j][2],
      );
      const r_buffered = r + S;

      // dE/dr — derivative of eq. (6), with the 1-4 scaling applied
      // to the SAME pairs the energy term scales.
      let dE_dr = -ELEC_UNIT * qi * qj / (r_buffered * r_buffered);
      if (is_1_4_pair(i, j, adj)) dE_dr *= SCALE_1_4;

      const { d_dx_a, d_dx_b } = bond_length_derivatives(positions[i], positions[j]);
      for (let axis = 0; axis < 3; axis++) {
        gradient[i][axis] += dE_dr * d_dx_a[axis];
        gradient[j][axis] += dE_dr * d_dx_b[axis];
      }
    }
  }

  return gradient;
}
