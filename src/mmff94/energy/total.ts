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
 * The scaling is applied inside the electrostatic term — it is the
 * only scaled term, and the term functions return totals (not pair
 * lists), so total.ts cannot rescale individual pairs. Each term
 * stays testable in isolation.
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

  // 1-4 scaling is applied inside the electrostatic term (×0.75 for
  // atoms exactly three bonds apart) — it is the only scaled term,
  // and the scaling cannot live here because the term functions
  // return totals, not pair lists. vdW is NOT scaled at 1-4
  // (Halgren 1996, p. 496 — see the header above).

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
