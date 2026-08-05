/**
 * mmff94 — the Merck Molecular Force Field in pure TypeScript.
 *
 * This barrel exports every public function a consumer needs to
 * assign atom types, compute energy, compute gradients, and optimize
 * geometry using the MMFF94 force field.
 *
 * Every function is pure — no hidden state, no global setup.
 */

export {
  assign_atom_types,
} from './assign-atom-types.js';

export {
  assign_bci_charges,
} from './charges.js';

export {
  calc_energy,
} from './energy/total.js';

export {
  calc_gradient,
} from './gradient/total.js';
