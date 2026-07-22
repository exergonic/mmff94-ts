/**
 * MMFF94 atom type assignment.
 *
 * This is the hardest single piece of the force field. Every atom in the
 * molecule must be assigned one of the ~350 MMFF94 atom types before
 * any energy term can look up its parameters. A wrong type cascades
 * into wrong parameters for EVERY energy term.
 *
 * The decision tree considers:
 *   - Element (C, N, O, H, S, P, halogen, etc.)
 *   - Coordination number (how many bonded neighbors, including implicit H)
 *   - Bond orders to neighbors (single, double, triple, aromatic)
 *   - What elements the neighbors are (e.g., C=O carbonyl vs. C=C alkene)
 *   - What types the neighbors have (e.g., carbonyl C vs. aromatic C)
 *   - Ring membership (is the atom itself in a ring?)
 *   - Formal charge
 *
 * The logic follows Halgren's 1996 description (J. Comput. Chem. 17, 520-552)
 * and is cross-checked against OpenBabel's atom typing.
 *
 * Each major branch in the decision tree is preceded by a comment
 * explaining WHY the chemistry requires that distinction.
 */

import type { Molecule, TypedMolecule } from '../types';

/**
 * Assign an MMFF94 atom type to every atom in the molecule.
 *
 * Returns a TypedMolecule with the same atoms and bonds plus
 * an atom_types array parallel to atoms[].
 *
 * If atom typing fails for any atom (e.g., unrecognized element
 * or coordination environment), that atom is assigned type 1
 * (generic sp³ carbon), which is a safe-but-approximate default.
 */
export function assign_atom_types(molecule: Molecule): TypedMolecule {
  // TODO: implement the full MMFF94 atom type decision tree.
  //
  // For now, this is a stub that assigns type 1 (generic sp³ C)
  // to every atom whose element is not H, and type 5 (generic H
  // bonded to sp³ C) to hydrogens.
  const atom_types: number[] = molecule.atoms.map(atom => {
    if (atom.element === 'H') return 5;
    if (atom.element === 'C') return 1;
    if (atom.element === 'N') return 8;
    if (atom.element === 'O') return 6;
    if (atom.element === 'S') return 15;
    if (atom.element === 'F') return 11;
    if (atom.element === 'Cl') return 12;
    if (atom.element === 'Br') return 13;
    if (atom.element === 'I') return 14;
    return 1; // fallback: generic sp³ C
  });

  return {
    ...molecule,
    atom_types,
  };
}

/**
 * Compute partial charges for every atom using the MMFF94 bond charge
 * increment (BCI) model.
 *
 * MMFF94 does NOT store per-atom partial charges in its parameter tables.
 * Instead, each BOND TYPE has a charge increment value. The partial charge
 * on an atom is the SUM of the BCI values of every bond it participates in.
 *
 * Formal charges (e.g., from an SDF or from a sketcher) override the BCI sum:
 * if an atom carries a formal charge, we use that directly instead.
 *
 * Mutates the TypedMolecule in place by setting partial_charges[].
 */
export function compute_bci_charges(molecule: TypedMolecule): void {
  // TODO: implement BCI charge calculation.
  //
  // Steps:
  // 1. For each bond, look up the BCI value by (type_i, type_j).
  // 2. Add +BCI to atom1's running total and -BCI to atom2's (or vice versa
  //    depending on the sign convention — must match OpenBabel's convention).
  // 3. Handle explicit formal charges: check if the molecule has formal charge
  //    info (from SDF or the sketcher) and override for those atoms.
  // 4. The sum of all partial charges should equal the net molecular charge.
  molecule.partial_charges = molecule.atoms.map(() => 0.0);
}
