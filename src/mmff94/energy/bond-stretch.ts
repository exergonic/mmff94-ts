/**
 * Bond stretching energy.
 *
 * MMFF94 uses a harmonic potential for bond stretching:
 *
 *   E_bond = 143.88 * k_b * (r − r₀)²
 *
 * where:
 *   k_b  = force constant in millidynes per Ångström (mdyn/Å)
 *   r    = current bond length in Å
 *   r₀   = equilibrium bond length in Å
 *   143.88 = unit conversion factor: (mdyn/Å) → (kcal/mol)/Å²
 *
 * A harmonic approximation is used instead of a full Morse potential
 * because for most organic molecules at room temperature, bonds do
 * not stretch far enough from equilibrium to feel the anharmonicity.
 * MMFF94 is parametrized for equilibrium geometries and conformational
 * energies, not bond dissociation.
 *
 * The parameters (k_b, r₀) are indexed by the MMFF94 types of the two
 * bonded atoms. If no exact (type_i, type_j) match is found, the lookup
 * falls through to a wildcard entry if one exists.
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
      total_energy += 143.88 * k_b * dr * dr;
    } else {
      // Fallback: If no parameters exist, we could add a zero contribution or throw an error.
      // For now, we just skip (add 0).
    }
  }

  return total_energy;
}
