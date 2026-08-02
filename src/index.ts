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
  compute_bci_charges,
  calc_energy,
  calc_gradient,
} from './mmff94';

// Geometry optimizers — L-BFGS is the primary minimizer. The
// steepest-descent fallback is not exported yet: it is still a stub
// (AGENTS.md Phase 6), and a placeholder that returns converged: false
// does not belong in the public API.
export {
  optimize_lbfgs,
} from './optimize/l-bfgs';

// Utility — for advanced users who want to build their own tools
export {
  distance,
  angle_in_radians,
  dihedral_angle,
  rotate_around_axis,
} from './utils/vector';
