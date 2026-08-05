// CSWAT probe: methanethiolate·water complex — the S⁻–H(water) bond
// (types 72-31) has NO parameter row in any implementation, so all
// three (mmff94-ts, Tinker, OpenBabel) generate it via the part V
// empirical rules (eqs. 18-19). This probe prints the model's
// generated parameters and the energy components for cross-implementation
// comparison. The typing is hand-set (identical to Tinker's classes).
import { make_class_context, bond_parameters } from '../../src/mmff94/parameters/parameter-classes';
import { empirical_bond_parameters } from '../../src/mmff94/parameters/empirical';
import { calc_energy } from '../../src/mmff94/energy/total';
import { assign_bci_charges } from '../../src/mmff94/charges';
import type { TypedMolecule } from '../../src/types';

const atoms = [
  { index: 0, element: 'C', x: 0.0, y: 0.0, z: 0.0 },
  { index: 1, element: 'H', x: 0.77, y: 1.02, z: -0.35 },
  { index: 2, element: 'H', x: -0.77, y: 0.75, z: 0.35 },
  { index: 3, element: 'H', x: 0.0, y: -0.35, z: 1.04 },
  { index: 4, element: 'S', x: 1.82, y: 0.0, z: 0.0 },
  { index: 5, element: 'O', x: 4.2, y: 0.6, z: 0.0 },
  { index: 6, element: 'H', x: 3.17, y: 0.33, z: 0.0 },
  { index: 7, element: 'H', x: 4.536, y: 1.502, z: 0.0 },
];
const bonds = [
  { atom1: 0, atom2: 1, bond_order: 1 },
  { atom1: 0, atom2: 2, bond_order: 1 },
  { atom1: 0, atom2: 3, bond_order: 1 },
  { atom1: 0, atom2: 4, bond_order: 1 }, // C–S⁻
  { atom1: 4, atom2: 6, bond_order: 1 }, // S⁻–H(water) — THE empirical bond
  { atom1: 5, atom2: 6, bond_order: 1 },
  { atom1: 5, atom2: 7, bond_order: 1 },
];
const mol: TypedMolecule = {
  name: 'CSWAT',
  atoms,
  bonds,
  atom_types: [1, 5, 5, 5, 72, 70, 31, 31],
};

const adj: number[][] = atoms.map(() => []);
for (const b of bonds) {
  adj[b.atom1].push(b.atom2);
  adj[b.atom2].push(b.atom1);
}
const ctx = make_class_context(mol, adj);

// The C–S⁻ bond (par row 0-1-72 exists) and the S⁻–H bond (par-less):
console.log('C-S-  par:', JSON.stringify(bond_parameters(ctx, 0, 4)));
const sH = bonds[4];
const emp = empirical_bond_parameters(mol.atoms[4], mol.atoms[6])!;
console.log('S-H   empirical r0 =', emp.r0.toFixed(5), ' k_b =', emp.k_b.toFixed(4));
console.log('S-H   r(S-H in geom) =', Math.hypot(3.17 - 1.82, 0.33, 0).toFixed(5));

const e = calc_energy(mol);
console.log('E     stretch =', e.bond_stretch.toFixed(6), ' angle =', e.angle_bend.toFixed(6));
console.log('      strbnd =', e.stretch_bend.toFixed(6), ' torsion =', e.torsion.toFixed(6),
  ' oop =', e.out_of_plane.toFixed(6));
console.log('      vdw =', e.van_der_waals.toFixed(6), ' elec =', e.electrostatic.toFixed(6),
  ' total =', e.total.toFixed(6));
const charged = assign_bci_charges(mol);
console.log('q(S) =', charged.partial_charges?.[4].toFixed(4),
  ' q(H7) =', charged.partial_charges?.[6].toFixed(4),
  ' q(O6) =', charged.partial_charges?.[5].toFixed(4));
