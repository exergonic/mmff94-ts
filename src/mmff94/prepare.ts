/**
 * One-call preparation for the simple API path.
 *
 * The top-level functions (calc_energy, calc_gradient, the
 * optimizers) accept a bare Molecule straight from parse_sdf() and
 * normalize it here: assign atom types if absent, then BCI charges if
 * absent. An already-prepared TypedMolecule passes through untouched
 * (no re-typing, no re-charging) — the expert path pays nothing.
 */

import type { Molecule, TypedMolecule } from '../types';
import { assign_atom_types } from './assign-atom-types';
import { assign_bci_charges } from './charges';

export function prepare_molecule(molecule: Molecule): TypedMolecule {
  let prepared = molecule as TypedMolecule;
  if (!prepared.atom_types) prepared = assign_atom_types(prepared);
  if (!prepared.partial_charges) prepared = assign_bci_charges(prepared);
  return prepared;
}
