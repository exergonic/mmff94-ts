/**
 * mmff94-ts quickstart — complete pipeline from SDF to minimized energy.
 *
 * Run with: npx tsx examples/quickstart.ts
 * (Requires tsx: npm install -g tsx or npx tsx)
 */

import { parse_sdf, assign_atom_types, compute_bci_charges, calc_energy } from '../src/index';

// A minimal ethane molecule in V2000 MOL format.
// This is the simplest test case: 2 carbons, 6 hydrogens, one C–C bond.
const ETHANE_SDF = `
  -ISIS-  0123456789

  8  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0900    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8900   -0.6300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3400    0.6300    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4300    0.6300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8800   -0.8900    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8800    0.0000   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
M  END
`;

function main() {
  // Step 1: parse the SDF string into our Molecule data model
  const molecule = parse_sdf(ETHANE_SDF);
  console.log(`Parsed ${molecule.atoms.length} atoms and ${molecule.bonds.length} bonds`);

  // Step 2: assign MMFF94 atom types to every atom
  const typed = assign_atom_types(molecule);
  console.log(`Atom types: ${typed.atom_types.join(', ')}`);

  // Step 3: compute partial charges from bond charge increments
  compute_bci_charges(typed);

  // Step 4: calculate the full MMFF94 energy
  const energy = calc_energy(typed);
  console.log(`\nMMFF94 Energy (kcal/mol):`);
  console.log(`  Bond stretch:     ${energy.bond_stretch.toFixed(4)}`);
  console.log(`  Angle bend:       ${energy.angle_bend.toFixed(4)}`);
  console.log(`  Stretch-bend:     ${energy.stretch_bend.toFixed(4)}`);
  console.log(`  Torsion:          ${energy.torsion.toFixed(4)}`);
  console.log(`  Van der Waals:    ${energy.van_der_waals.toFixed(4)}`);
  console.log(`  Electrostatic:    ${energy.electrostatic.toFixed(4)}`);
  console.log(`  Out-of-plane:     ${energy.out_of_plane.toFixed(4)}`);
  console.log(`  ─────────────────────────`);
  console.log(`  TOTAL:            ${energy.total.toFixed(4)}`);
}

main();
