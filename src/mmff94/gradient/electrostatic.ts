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

import type { TypedMolecule } from '../../types.js';
import { Vec3 } from '../../utils/vector.js';
import { assign_bci_charges } from '../charges.js';
import { nonbonded_context_for } from '../nonbonded-context.js';
import { bond_length_derivatives } from './derivatives.js';

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

  // Pair list + 1-4 classification come from the cached context —
  // they are topology, not geometry (see nonbonded-context.ts). The
  // context's exactly-depth-3 flag is the same set is_1_4_pair
  // returns (same BFS, same shortest-path-wins semantics); that
  // function stays exported in energy/electrostatic.ts for callers
  // holding only an adjacency list.
  const ctx = nonbonded_context_for(molecule);

  for (let p = 0; p < ctx.n_pairs; p++) {
    const i = ctx.pair_i[p];
    const j = ctx.pair_j[p];
    const qi = charges[i];
    const qj = charges[j];
    if (qi === 0 || qj === 0) continue;

    const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];
    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];

    const r = Math.hypot(posI[0] - posJ[0], posI[1] - posJ[1], posI[2] - posJ[2]);
    const r_buffered = r + S;

    // dE/dr — derivative of eq. (6), with the 1-4 scaling applied
    // to the SAME pairs the energy term scales.
    let dE_dr = -ELEC_UNIT * qi * qj / (r_buffered * r_buffered);
    if (ctx.pair_is_14[p]) dE_dr *= SCALE_1_4;

    const { d_dx_a, d_dx_b } = bond_length_derivatives(posI, posJ);
    for (let axis = 0; axis < 3; axis++) {
      gradient[i][axis] += dE_dr * d_dx_a[axis];
      gradient[j][axis] += dE_dr * d_dx_b[axis];
    }
  }

  return gradient;
}
