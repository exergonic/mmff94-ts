/**
 * Out-of-plane bending energy.
 *
 * Halgren1996, eq. (6):
 *
 *   E_oop = 0.043844 · (k_oop / 2) · χ²
 *
 * where:
 *   k_oop  = out-of-plane force constant (mdyn·Å/rad²)
 *   χ      = Wilson out-of-plane angle (degrees) between the bond j→l
 *            and the plane through (i, j, k) — see wilson_oop_angle()
 *   0.043844 = unit conversion factor: same as angle bending
 *
 * MMFF94 uses this dedicated term for TRICOORDINATE centers, not the
 * "improper torsion" approach used by UFF and GAFF. Those are NOT the
 * same thing: an improper torsion treats the deformation as a dihedral
 * rotation, which couples it with the true torsion term, while the oop
 * term measures how far one bond sticks out of the plane of the other
 * two — a pure bending coordinate.
 *
 * The term is NOT restricted to planar (sp²) centers. The three angles
 * at a given center all share one force constant, and the angle-bending
 * term uses reference values that may average less than 120°; together
 * these let MMFF94 represent pyramidal equilibrium geometries and fit
 * inversion barriers (Halgren 1996). What actually decides whether a
 * center contributes is the parameter table:
 *
 *   - Planar sp² centers (carbonyl C, vinylic C, aromatic C) have
 *     substantial positive k_oop that keeps them planar.
 *   - Amine N (type 8) has k_oop = 0: pyramidalization is handled
 *     entirely by the angle-bend reference angles.
 *   - Amide N (type 10) has NEGATIVE k_oop, which actively favors the
 *     nonplanar geometry — this is why MMFF94 (unlike MMFF94s) gives
 *     pyramidal delocalized trigonal nitrogens.
 *   - Tetravalent atoms (e.g. alkane C) have no oop parameters at all.
 *
 * The force constant is keyed by the central atom type AND the three
 * substituent types (key i-j-k-l, central j in the second position).
 * The three substituent types are interchangeable — the same k applies
 * to all three angles at the center — so the lookup matches the sorted
 * multiset of substituent types, falling back to the per-central-type
 * wildcard entry ("0-j-0-0") when no specific entry exists.
 */

import type { TypedMolecule } from '../../types';
import { OOP_PARAMS, lookup_param } from '../parameters';
import { wilson_oop_angle, Vec3 } from '../../utils/vector';

const OOP_UNIT = 0.043844; // (mdyn·Å/rad²)·deg² → kcal/mol, same as angle bending

/**
 * The out-of-plane force constant for the tri-coordinate center j
 * with substituents a, c, d.
 *
 * The three substituent types are interchangeable — the same k
 * applies to all three angles at the center — so the lookup matches
 * the sorted multiset of substituent types, falling back to the
 * per-central-type wildcard entry ("0-j-0-0") when no specific entry
 * exists. Shared with the gradient so both terms use the same k.
 */
export function oop_force_constant(
  molecule: TypedMolecule,
  j: number,
  a: number,
  c: number,
  d: number,
): number | undefined {
  const tj = molecule.atom_types[j];

  const sorted = [
    molecule.atom_types[a],
    molecule.atom_types[c],
    molecule.atom_types[d],
  ].sort((x, y) => x - y);
  let params = lookup_param(OOP_PARAMS, [sorted[0], tj, sorted[1], sorted[2]]);

  // No specific entry: fall back to the per-central-type wildcard
  // ("0-j-0-0"). This is what applies amine N's explicit zero and
  // amide N's negative constant to every substituent combination.
  if (!params) {
    params = lookup_param(OOP_PARAMS, [0, tj, 0, 0]);
  }
  return params?.k_oop;
}

/**
 * Calculate the total out-of-plane bending energy.
 *
 * Evaluated for every tri-coordinate center (exactly three bonded
 * neighbors). Each center contributes three Wilson angles — one for
 * each substituent taking its turn as the out-of-plane atom l, with
 * the plane defined by the other two — all sharing one k_oop.
 */
export function calc_oop_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  // Build adjacency list
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }

  for (let j = 0; j < molecule.atoms.length; j++) {
    const neighbors = adj[j];
    if (neighbors.length !== 3) continue;

    const [a, c, d] = neighbors;

    // The three substituents are interchangeable — the same k_oop
    // applies to all three Wilson angles at the center (resolved by
    // the shared oop_force_constant() helper, same for the gradient).
    const k_oop = oop_force_constant(molecule, j, a, c, d);
    if (k_oop === undefined) continue;

    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];
    const posA: Vec3 = [molecule.atoms[a].x, molecule.atoms[a].y, molecule.atoms[a].z];
    const posC: Vec3 = [molecule.atoms[c].x, molecule.atoms[c].y, molecule.atoms[c].z];
    const posD: Vec3 = [molecule.atoms[d].x, molecule.atoms[d].y, molecule.atoms[d].z];

    // The three Wilson angles at j: each substituent takes a turn as
    // the out-of-plane atom, with the plane through j and the other two.
    const chi_a = wilson_oop_angle(posD, posJ, posC, posA);  // a out of plane (d, c)
    const chi_c = wilson_oop_angle(posA, posJ, posD, posC);  // c out of plane (a, d)
    const chi_d = wilson_oop_angle(posA, posJ, posC, posD);  // d out of plane (a, c)

    total_energy +=
      OOP_UNIT * (k_oop / 2.0) * (chi_a * chi_a + chi_c * chi_c + chi_d * chi_d);
  }

  return total_energy;
}
