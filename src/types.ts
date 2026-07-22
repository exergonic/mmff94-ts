/**
 * Portable data model for a chemical molecule.
 *
 * These types are the common currency of the entire library.
 * Every function accepts or returns one of these shapes.
 * No external dependency — just plain objects a chemist can read.
 */

/** A single atom with its 3D position in Ångströms. */
export interface Atom {
  index: number;
  element: string;     // Periodic-table symbol: 'C', 'N', 'O', 'H', 'S', etc.
  x: number;
  y: number;
  z: number;
}

/** A bond connecting two atoms by their index in the atoms array. */
export interface Bond {
  atom1: number;
  atom2: number;
  bond_order: number;  // 1 = single, 2 = double, 3 = triple
}

/** A molecule — the minimum needed to evaluate a force field. */
export interface Molecule {
  atoms: Atom[];
  bonds: Bond[];
  name?: string;       // Optional human-readable label
}

/**
 * A molecule that has been through atom typing.
 *
 * Every atom now has an MMFF94 type number, which is the key that
 * unlocks the right parameters for every energy term.
 * Partial charges are filled in by compute_bci_charges().
 */
export interface TypedMolecule extends Molecule {
  atom_types: number[];          // MMFF94 type index per atom, same order as atoms[]
  partial_charges?: number[];    // Partial charge (e⁻) per atom, from BCI model
}

/**
 * The seven energy components that make up the total MMFF94 energy.
 * All values are in kcal/mol.
 */
export interface EnergyComponents {
  total: number;
  bond_stretch: number;
  angle_bend: number;
  stretch_bend: number;       // Cross term — couples bond and angle
  torsion: number;
  van_der_waals: number;
  electrostatic: number;
  out_of_plane: number;
}

/** Result of geometry optimization. */
export interface OptimizationResult {
  molecule: Molecule;          // Optimized geometry, atoms in new positions
  energy: EnergyComponents;     // Final energy after optimization
  iterations: number;           // Number of optimization cycles
  converged: boolean;           // True if max gradient fell below the threshold
  final_max_gradient: number;   // Largest force component (kcal/mol/Å) at finish
}
