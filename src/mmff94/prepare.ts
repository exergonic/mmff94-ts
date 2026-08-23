/**
 * One-call preparation for the simple API path.
 *
 * The top-level functions (calc_energy, calc_gradient, the
 * optimizers) accept a bare Molecule straight from parse_sdf() and
 * normalize it here: assign atom types if absent, then BCI charges if
 * absent. An already-prepared TypedMolecule passes through untouched
 * (no re-typing, no re-charging) — the expert path pays nothing.
 */

import type { Molecule, TypedMolecule } from '../types.js';
import { assign_atom_types } from './assign-atom-types.js';
import { assign_bci_charges } from './charges.js';

export function prepare_molecule(molecule: Molecule): TypedMolecule {
  let prepared = molecule as TypedMolecule;
  // Length-checked guards: an empty or stale annotation array must not
  // skip typing/charging silently (a re-parsed or hand-edited Molecule
  // can carry leftovers). A right-length-but-wrong array is
  // undetectable here — the expert path owns that contract.
  const needs_typing =
    !prepared.atom_types || prepared.atom_types.length !== molecule.atoms.length;
  if (needs_typing) prepared = assign_atom_types(prepared);
  const needs_charging =
    !prepared.partial_charges ||
    prepared.partial_charges.length !== molecule.atoms.length;
  if (needs_charging) prepared = assign_bci_charges(prepared);
  return prepared;
}
