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
  assign_bci_charges,
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
  // ── SIMPLE PATH: a parsed molecule is all you need ──────────────
  const molecule = parse_sdf(ETHANE_SDF);
  console.log(`Parsed ${molecule.atoms.length} atoms and ${molecule.bonds.length} bonds`);

  // Energy with the full per-term breakdown — atom typing and BCI
  // charges happen on demand. The simple path is never a subset: it
  // returns exactly what the rich path returns.
  const energy = calc_energy(molecule);
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

  // Geometry optimization — one call, no preparation, no callback.
  // The SDF geometry is already MMFF94-optimized, so nudge one
  // hydrogen first: the optimizer must walk it back to the staggered
  // minimum. The result carries the full per-term energy AND the
  // typed/charged molecule at the minimum.
  molecule.atoms[2].x += 0.3;
  const result = optimize_lbfgs(molecule, { gradient_tolerance: 0.05 });
  console.log(`\nOptimization (L-BFGS, ${result.iterations} iterations):`);
  console.log(`  Energy: ${energy.total.toFixed(4)} → ${result.energy.total.toFixed(4)} kcal/mol`);
  console.log(`  Converged: ${result.converged} (max|dE/dx| = ${result.final_max_gradient.toExponential(1)})`);

  // ── RICH PATH: every layer explicitly, like OpenBabel ───────────
  // The same molecule, step by step — for inspecting the types, the
  // charges, and for handing a custom oracle to the optimizer.
  const typed = assign_atom_types(molecule);
  console.log(`\nAtom types: ${typed.atom_types.join(', ')}`);

  const charged = assign_bci_charges(typed);
  console.log(
    `Partial charges: ${charged.partial_charges!.map(q => q.toFixed(3)).join(', ')}`,
  );

  const gradient = calc_gradient(charged);
  console.log(`Gradient: max|dE/dx| = ${Math.max(...gradient.flat().map(Math.abs)).toFixed(4)} kcal/mol/Å`);

  // Custom oracle: the rich path keeps full control of the optimizer.
  const custom = optimize_lbfgs(
    charged,
    m => ({ energy: calc_energy(m), gradient: calc_gradient(m) }),
    { gradient_tolerance: 0.05 },
  );
  console.log(`\nOptimization with a custom oracle: ${custom.iterations} iterations, E = ${custom.energy.total.toFixed(4)} kcal/mol`);

  // The simple path's result molecule is fully prepared — the rich
  // path is reachable from it without re-doing anything.
  console.log(`Types at the minimum: ${result.molecule.atom_types.join(', ')}`);
}

main();
