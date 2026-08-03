/**
 * mmff94-ts quickstart — the full pipeline from SDF to per-term MMFF94 energy,
 * gradient, and L-BFGS geometry optimization.
 *
 * Run with: npx tsx examples/quickstart.ts
 * (Requires tsx: npm install -g tsx or npx tsx)
 */

import {
  parse_sdf,
  assign_atom_types,
  compute_bci_charges,
  calc_energy,
  calc_gradient,
  optimize_lbfgs,
} from '../src/index';

// A minimal ethane molecule in V2000 MOL format: 2 carbons, 6 hydrogens,
// one C–C bond. The geometry is the STAGGERED (gauche 60°) global minimum,
// MMFF94-optimized — the torsion term is still nonzero here, because
// MMFF94's H-C-C-H Fourier terms (V1/V2) don't vanish at 60°.
const ETHANE_SDF = `ethane
  produced by Avogadro 07242607333D; MMFF94 optimized

  8  7  0  0  0  0  0  0  0  0999 V2000
    0.7560    0.0000   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7560   -0.0000   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1404   -0.5122    0.8871 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1404   -0.5122   -0.8871 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1404    1.0244    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1404    0.5122    0.8871 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1404    0.5122   -0.8871 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1404   -1.0244    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
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

  // Step 3: compute partial charges from bond charge increments.
  // compute_bci_charges returns a COPY of the molecule with the
  // charges attached — the value that flows into the energy terms
  // (which also compute the charges on demand if given a bare typed
  // molecule).
  const charged = compute_bci_charges(typed);

  // Step 4: calculate the full MMFF94 energy
  const energy = calc_energy(charged);
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

  // Step 5: analytical gradient — dE/dx_i per atom (kcal/mol/Å)
  const gradient = calc_gradient(typed);
  const max_g = Math.max(...gradient.flat().map(Math.abs));
  console.log(`\nGradient: max|dE/dx| = ${max_g.toFixed(4)} kcal/mol/Å`);

  // Step 6: L-BFGS geometry optimization. The SDF geometry is already
  // MMFF94-optimized, so nudge one hydrogen first — the optimizer must
  // walk it back to the staggered minimum.
  typed.atoms[2].x += 0.3;
  const e_before = calc_energy(typed).total;
  const result = optimize_lbfgs(
    typed,
    m => ({ energy: calc_energy(m), gradient: calc_gradient(m) }),
    { gradient_tolerance: 0.05 },
  );
  console.log(`\nOptimization (L-BFGS, ${result.iterations} iterations):`);
  console.log(`  Energy: ${e_before.toFixed(4)} → ${result.energy.total.toFixed(4)} kcal/mol`);
  console.log(`  Converged: ${result.converged} (max|dE/dx| = ${result.final_max_gradient.toExponential(1)})`);
}

main();
