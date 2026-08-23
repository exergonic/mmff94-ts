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
  /** Formal charge in electrons (0 = neutral). Read from the SDF/MOL
   * charge field when present; used by atom typing (charged type
   * variants) and by the BCI charge model's primary charges. */
  formal_charge?: number;
  /** Reference partial charge in electrons, when the structure file
   * carries one (the validation suite's .mmd pchg column). Not used by
   * any energy computation — the library computes BCI charges itself;
   * this field exists so the suite parser can expose the reference for
   * comparison. */
  partial_charge?: number;
  /** Atom label from the structure file (the suite's .mmd name field,
   * e.g. "O6", "FE1") — how the reference files address individual
   * atoms (the formal-charge file keys on it). */
  label?: string;
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
 * Partial charges are attached by assign_bci_charges() (which
 * returns the charged molecule); the energy terms fall back to
 * computing them on demand when the field is absent.
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
  /** Optimized geometry — typed and charged (the optimizer works on a
   *  prepared copy, so the result carries atom_types and
   *  partial_charges for follow-up per-term work at the minimum). */
  molecule: TypedMolecule;
  energy: EnergyComponents;     // Final energy after optimization
  iterations: number;           // Number of optimization cycles
  converged: boolean;          // True if a convergence gate tripped (see the optimizer options)
  final_max_gradient: number;   // Largest force component (kcal/mol/Å) at finish
  /** RMS of the force components (kcal/mol/Å) at finish — the TINKER-style
   *  convergence signal; NaN when the molecule has no atoms. */
  final_rms_gradient?: number;
}
