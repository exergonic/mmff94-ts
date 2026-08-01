/**
 * Total MMFF94 energy.
 *
 * Sums all seven energy terms.
 *
 * 1-4 SCALING applies to atoms that are exactly three bonds apart:
 *   - Electrostatic:      multiply by 0.75
 *   - Van der Waals:      NOT scaled — Halgren 1996 (p. 496): "1,4-vdW
 *                         interactions are not differentially scaled in
 *                         MMFF94". The ×0.5 common in MM2/GAFF is a
 *                         different force field's convention.
 *
 * The electrostatic factor is part of the MMFF94 specification and
 * compensates for the fact that 1-4 charges are partially captured by
 * the torsion term's parameterization.
 *
 * STATUS: the scaling is NOT yet applied — all terms are summed
 * unscaled (see TODO below). When it lands, it will be applied HERE
 * rather than in the individual term functions, so that each term
 * function is simple and testable in isolation: the individual term
 * functions calculate the FULL (unscaled) energy for every pair, and
 * total.ts decides which pairs get scaled.
 */

import type { TypedMolecule, EnergyComponents } from '../../types';
import { calc_bond_stretch_energy } from './bond-stretch';
import { calc_angle_bend_energy } from './angle-bend';
import { calc_stretch_bend_energy } from './stretch-bend';
import { calc_torsion_energy } from './torsion';
import { calc_vdw_energy } from './van-der-waals';
import { calc_electrostatic_energy } from './electrostatic';
import { calc_oop_energy } from './out-of-plane';

/**
 * Compute the full MMFF94 energy for a typed molecule.
 *
 * Returns an EnergyComponents object with every term broken out
 * separately plus the total.
 */
export function calc_energy(molecule: TypedMolecule): EnergyComponents {
  const bond_stretch  = calc_bond_stretch_energy(molecule);
  const angle_bend    = calc_angle_bend_energy(molecule);
  const stretch_bend  = calc_stretch_bend_energy(molecule);
  const torsion       = calc_torsion_energy(molecule);
  const van_der_waals = calc_vdw_energy(molecule);
  const electrostatic = calc_electrostatic_energy(molecule);
  const out_of_plane  = calc_oop_energy(molecule);

  // TODO: apply 1-4 scaling to electrostatic only (×0.75). vdW is NOT
  // scaled at 1-4 (Halgren 1996, p. 496 — see the header above).
  // This requires identifying 1-4 atom pairs (atoms exactly 3 bonds apart)
  // and multiplying the appropriate pairwise contributions by 0.75.
  // TODO: stretch-bend interactions are omitted when Halgren's equation 4
  // is used to calculate the angle bending energy. In angle-bend.ts, we use
  // equation 4 when the reference angle theta0 > 150 degrees. Keep track of
  // which angles use equation 4 and which use equation 3, OR use the same
  // condition for both equations. If the angle uses equation 4, simply omit
  // stretch-bend interactions for it.

  const total = bond_stretch + angle_bend + stretch_bend +
                torsion + van_der_waals + electrostatic + out_of_plane;

  return {
    total,
    bond_stretch,
    angle_bend,
    stretch_bend,
    torsion,
    van_der_waals,
    electrostatic,
    out_of_plane,
  };
}
