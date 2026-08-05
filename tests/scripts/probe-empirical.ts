// OHMW1 empirical bond-generation probe: the eq. (18)-(19) parameters
// for the hydroxide O-H (the suite's only par-less bond) and the
// stretch contribution they produce.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { calc_energy } from '../../src/mmff94/energy/total';
import { make_class_context } from '../../src/mmff94/parameters/parameter-classes';
import { empirical_bond_parameters } from '../../src/mmff94/parameters/empirical';
import { distance, Vec3 } from '../../src/utils/vector';

const name = process.argv[2] ?? 'OHMW1';
const text = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');
const mol = parse_mmd(text).find(m => m.name === name)!;
const typed = assign_atom_types(mol);

const adj: number[][] = Array.from({ length: mol.atoms.length }, () => []);
for (const bond of mol.bonds) {
  adj[bond.atom1].push(bond.atom2);
  adj[bond.atom2].push(bond.atom1);
}
const ctx = make_class_context(typed, adj);

console.log('per-atom: idx elem type');
for (let a = 0; a < mol.atoms.length; a++) {
  console.log(`  ${a} ${mol.atoms[a].element} t${typed.atom_types[a]}`);
}

const comp = calc_energy(typed);
console.log(`\nmodel stretch (current, O-H skipped): ${comp.bond_stretch.toFixed(6)}`);
console.log('reference BOND STRETCHING (bmin): 0.5756 (need to check)');

// For every bond lacking a stored parameter, print the empirical values.
for (const bond of mol.bonds) {
  const emp = empirical_bond_parameters(ctx, typed, bond);
  if (!emp) continue;
  const a1 = typed.atoms[bond.atom1];
  const a2 = typed.atoms[bond.atom2];
  const r = distance(
    [a1.x, a1.y, a1.z] as Vec3,
    [a2.x, a2.y, a2.z] as Vec3,
  );
  const dr = r - emp.r0;
  const cs = -2.0;
  const harmonic = 143.9325 * 0.5 * emp.k_b * dr * dr;
  const anharmonic = 1.0 + cs * dr + (7.0 / 12.0) * cs * cs * dr * dr;
  const e = harmonic * anharmonic;
  console.log(
    `\nbond ${bond.atom1}-${bond.atom2} (${a1.element}${typed.atom_types[bond.atom1]}-${a2.element}${typed.atom_types[bond.atom2]}) ` +
      `r=${r.toFixed(4)}  r0_emp=${emp.r0.toFixed(4)}  k_emp=${emp.k_b.toFixed(3)}  E=${e.toExponential(3)}`,
  );
  console.log(`  new stretch total: ${(comp.bond_stretch + e).toFixed(6)}`);
}
