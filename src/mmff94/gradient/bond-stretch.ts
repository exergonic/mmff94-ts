/**
 * Gradient of the bond stretching energy.
 *
 * Halgren1996, eq. (2) — see energy/bond-stretch.ts for the energy:
 *
 *   E_bond = 143.9325 · (k_b/2) · Δr² · [1 + cs·Δr + 7/12·cs²·Δr²]
 *
 * with Δr = r − r₀ and cs = −2 Å⁻¹. The chain rule gives
 *
 *   dE/dr = 143.9325 · k_b · Δr · [1 + cs·Δr + 7/12·cs²·Δr²]
 *         + 143.9325 · (k_b/2) · Δr² · [cs + 7/6·cs²·Δr]
 *
 * and the geometric factor dr/dx_i = (pos_i − pos_j)/r is the unit
 * vector along the bond (shared via derivatives.ts, which mirrors the
 * energy term's distance() path).
 *
 * The parameter lookup must be identical to the energy term's —
 * including the class-1 (BTij) branch for conjugated single bonds —
 * or a bond the energy counts would contribute a different force here.
 */

import type { TypedMolecule } from '../../types.js';
import { Vec3 } from '../../utils/vector.js';
import { bond_parameters, class_context_for } from '../parameters/parameter-classes.js';
import { empirical_bond_parameters } from '../parameters/empirical.js';
import { bond_length_derivatives } from './derivatives.js';

const BOND_UNIT = 143.9325; // (mdyn/Å) → (kcal/mol)/Å²
const CS = -2.0; // cubic stretch constant, Å⁻¹

/**
 * Gradient of the bond stretching energy, dE/dx per atom.
 * Returns an array parallel to molecule.atoms[]:
 *   result[i] = [dE/dx_i, dE/dy_i, dE/dz_i]  (kcal/mol/Å)
 */
export function calc_bond_stretch_gradient(molecule: TypedMolecule): number[][] {
  const gradient: number[][] = molecule.atoms.map(() => [0, 0, 0]);

  // Adjacency for the BTij class queries — same as the energy term.
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }
  const ctx = class_context_for(molecule, adj);

  for (const bond of molecule.bonds) {
    const a1 = bond.atom1;
    const a2 = bond.atom2;
    const pos1: Vec3 = [molecule.atoms[a1].x, molecule.atoms[a1].y, molecule.atoms[a1].z];
    const pos2: Vec3 = [molecule.atoms[a2].x, molecule.atoms[a2].y, molecule.atoms[a2].z];

    let params = bond_parameters(ctx, a1, a2);
    if (!params) {
      // Mirror the energy term's miss path: the part V empirical rules
      // generate parameters for a bond with no stored row (e.g. the
      // sp2-C-P bond of vinyl phosphine — Halgren's Table III has no
      // 0-2-26 row). Without this the bond's gradient silently
      // vanishes while its energy is nonzero — the optimizer sees no
      // restoring force on the bond and collapses it (found via
      // dogfooding: P-C ended at 0.957 Å vs the 1.83 equilibrium).
      params = empirical_bond_parameters(molecule.atoms[a1], molecule.atoms[a2]);
      if (!params) continue;
    }

    // dE/dr — the derivative of eq. (2) w.r.t. r (see header)
    const { d_dx_a, d_dx_b } = bond_length_derivatives(pos1, pos2);
    const r = Math.hypot(pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2] - pos2[2]);
    const dr = r - params.r0;
    const k_b = params.k_b;
    const anharmonic = 1.0 + CS * dr + (7.0 / 12.0) * CS * CS * dr * dr;
    const d_anharmonic = CS + (7.0 / 6.0) * CS * CS * dr;
    const dE_dr = BOND_UNIT * k_b * dr * anharmonic + BOND_UNIT * (k_b / 2.0) * dr * dr * d_anharmonic;

    for (let a = 0; a < 3; a++) {
      gradient[a1][a] += dE_dr * d_dx_a[a];
      gradient[a2][a] += dE_dr * d_dx_b[a];
    }
  }

  return gradient;
}
