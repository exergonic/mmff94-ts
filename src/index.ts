/**
 * mmff94-ts — Merck Molecular Force Field in pure TypeScript.
 *
 * Usage:
 *
 *   import { parse_sdf, assign_atom_types, calc_energy } from 'mmff94-ts';
 *
 *   const mol = parse_sdf(sdf_text);
 *   const typed = assign_atom_types(mol);
 *   const energy = calc_energy(typed);
 *   console.log(`Total energy: ${energy.total} kcal/mol`);
 *
 * See examples/quickstart.ts for a complete walk-through.
 */

// Types — the data model everything else builds on
export type {
  Atom,
  Bond,
  Molecule,
  TypedMolecule,
  EnergyComponents,
  OptimizationResult,
} from './types';

// SDF parser — read molecules from standard file formats
export { parse_sdf } from './sdf';

// MMFF94 force field
export {
  assign_atom_types,
  assign_bci_charges,
  calc_energy,
  calc_gradient,
} from './mmff94';

// Geometry optimizers — L-BFGS is the primary minimizer; steepest
// descent is the fallback (robust but slow — Armijo line search).
export { optimize_lbfgs } from './optimize/l-bfgs';
export { optimize_steepest_descent } from './optimize/steepest-descent';

// Utility — for advanced users who want to build their own tools
export {
  distance,
  angle_in_radians,
  dihedral_angle,
  rotate_around_axis,
} from './utils/vector';
