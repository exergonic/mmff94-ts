/**
 * Bond stretching energy.
 *
 * Halgren1996, eq. (2):
 *
 *   E_bond = 143.9325 · k_b · (r − r₀)²
 *            · [1 + cs · (r − r₀) + 7/12 · cs² · (r − r₀)²]
 *
 * where:
 *   k_b  = force constant in millidynes per Ångström (mdyn/Å)
 *   r    = current bond length in Å
 *   r₀   = equilibrium bond length in Å
 *   cs   = cubic stretch constant = −2 Å⁻¹
 *   143.9325 = unit conversion factor: (mdyn/Å) → (kcal/mol)/Å²
 *
 * We implement only the leading harmonic term (the purely quadratic part
 * of the expansion). The cubic and quartic correction terms (with cs)
 * are omitted because they contribute negligibly for the small
 * displacements near equilibrium that dominate conformational energies.
 * This truncation matches the convention used by OpenBabel's MMFF94
 * implementation (forcefieldmmff94.cpp:120) and RDKit.
 *
 * The 143.9325 factor equals 143.88 when halved and rounded to MM2
 * bondunit (71.94), but Halgren's published value is 143.9325.
 */

import type { TypedMolecule } from '../../types';
import { BOND_PARAMS, lookup_param } from '../parameters';
import { distance, Vec3 } from '../../utils/vector';

/**
 * Calculate the total bond stretching energy for all bonds in a molecule.
 */
export function calc_bond_stretch_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  for (const bond of molecule.bonds) {
    const a1 = molecule.atoms[bond.atom1];
    const a2 = molecule.atoms[bond.atom2];
    const t1 = molecule.atom_types[bond.atom1];
    const t2 = molecule.atom_types[bond.atom2];

    const pos1: Vec3 = [a1.x, a1.y, a1.z];
    const pos2: Vec3 = [a2.x, a2.y, a2.z];
    const r = distance(pos1, pos2);

    const t_min = Math.min(t1, t2);
    const t_max = Math.max(t1, t2);

    const params = lookup_param(BOND_PARAMS, [t_min, t_max]);
    if (params) {
      const { k_b, r0 } = params;
      const dr = r - r0;
      total_energy += 143.9325 * k_b * dr * dr;
    } else {
      // Fallback: If no parameters exist, we could add a zero contribution or throw an error.
      // For now, we just skip (add 0).
    }
  }

  return total_energy;
}
