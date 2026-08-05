/**
 * Bond stretching energy.
 *
 * Halgren1996, eq. (2):
 *
 *   E_bond = 143.9325 · k_b/2 · (r − r₀)²
 *            · [1 + cs · (r − r₀) + 7/12 · cs² · (r − r₀)²]
 *
 * where:
 *   k_b  = force constant in millidynes per Ångström (mdyn/Å)
 *   r    = current bond length in Å
 *   r₀   = equilibrium bond length in Å
 *   cs   = cubic stretch constant = −2 Å⁻¹
 *   143.9325 = unit conversion factor: (mdyn/Å) → (kcal/mol)/Å²
 *
 */

import type { TypedMolecule } from '../../types';
import { distance, Vec3 } from '../../utils/vector';
import { make_class_context, bond_parameters } from '../parameters/parameter-classes';
import { empirical_bond_parameters } from '../parameters/empirical';

/**
 * Calculate the total bond stretching energy for all bonds in a molecule.
 */
export function calc_bond_stretch_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  // Adjacency for the BTij class queries (conjugated single bonds).
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }
  const ctx = make_class_context(molecule, adj);

  for (const bond of molecule.bonds) {
    const a1 = molecule.atoms[bond.atom1];
    const a2 = molecule.atoms[bond.atom2];

    const pos1: Vec3 = [a1.x, a1.y, a1.z];
    const pos2: Vec3 = [a2.x, a2.y, a2.z];
    const r = distance(pos1, pos2);

    // Class-aware lookup: a BTij=1 bond (conjugated single bond) uses
    // the class-1 entry — '0-2-2' is the C=C double-bond parameter and
    // would badly overestimate a diene's central single bond.
    let params = bond_parameters(ctx, bond.atom1, bond.atom2);
    if (!params) {
      // Part V empirical-rule generation (eqs. 18-19): the designed
      // fallback for a bond with no stored row — the suite's only case
      // is OHMW1's hydroxide O–H. The generation is itself validated:
      // OHMW1's stretch now matches the reference to 6e-5 (the two
      // measured constants live in empirical.ts).
      params = empirical_bond_parameters(molecule.atoms[bond.atom1], molecule.atoms[bond.atom2]);
      if (!params) continue;
    }

    const { k_b, r0 } = params;
    const dr = r - r0;
    const cs = -2.0;

    // Halgren1996 eq. (2): E = 143.9325 · (k_b/2) · Δr² · [1 + cs·Δr + 7/12·cs²·Δr²]
    const half_k_b = 0.5 * k_b;
    const harmonic = 143.9325 * half_k_b * dr * dr;
    const anharmonic = 1.0 + cs * dr + (7.0 / 12.0) * cs * cs * dr * dr;
    total_energy += harmonic * anharmonic;
  }

  return total_energy;
}
