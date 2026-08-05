// TAJSUS BT-flag / torsion-class ground truth: which atoms are 78-81,
// what the properties say, and how the 80-81 bond gets class 1.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';
import { make_class_context, bond_type_flag, torsion_class } from '../../src/mmff94/parameters/parameter-classes';
import { ATOM_TYPE_PROPERTIES } from '../../src/mmff94/parameters';

const name = process.argv[2] ?? 'TAJSUS';
const text = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');
const mol = parse_mmd(text).find(m => m.name === name)!;
const typed = assign_atom_types(mol);

const adj: number[][] = Array.from({ length: mol.atoms.length }, () => []);
for (const bond of mol.bonds) {
  adj[bond.atom1].push(bond.atom2);
  adj[bond.atom2].push(bond.atom1);
}
const ctx = make_class_context(typed, adj);

console.log('per-atom: idx elem type arom sbmb  nbrs(order:type)');
for (let a = 0; a < mol.atoms.length; a++) {
  const t = typed.atom_types[a];
  const p = ATOM_TYPE_PROPERTIES[t];
  const nbrs = adj[a]
    .map(n => {
      const b = mol.bonds.find(b => (b.atom1 === a && b.atom2 === n) || (b.atom1 === n && b.atom2 === a))!;
      return `${n}(${b.bond_order}:${typed.atom_types[n]})`;
    })
    .join(' ');
  console.log(
    `${String(a).padStart(2)} ${mol.atoms[a].element.padEnd(2)} t${String(t).padStart(2)} ` +
      `${p?.arom ? 'arom ' : '     '}${p?.sbmb ? 'sbmb ' : '     '} ${nbrs}`,
  );
}

console.log('\nBT flags for the N-chain bonds:');
for (const bond of mol.bonds) {
  const t1 = typed.atom_types[bond.atom1];
  const t2 = typed.atom_types[bond.atom2];
  const interesting = [78, 79, 80, 81, 9].includes(t1) && [78, 79, 80, 81, 9].includes(t2);
  if (interesting) {
    console.log(
      `bond ${bond.atom1}-${bond.atom2} (t${t1}-t${t2}) order=${bond.bond_order} ` +
        `BT=${bond_type_flag(ctx, bond.atom1, bond.atom2)}`,
    );
  }
}

console.log('\nclasses for the 80-81 dihedrals:');
for (const bond of mol.bonds) {
  if (!(typed.atom_types[bond.atom1] === 80 && typed.atom_types[bond.atom2] === 81) &&
      !(typed.atom_types[bond.atom1] === 81 && typed.atom_types[bond.atom2] === 80)) continue;
  const j = bond.atom1;
  const k = bond.atom2;
  for (const i of adj[j].filter(n => n !== k)) {
    for (const l of adj[k].filter(n => n !== j)) {
      if (i === l) continue;
      console.log(
        `dihedral ${i}-${j}-${k}-${l} t[${typed.atom_types[i]},${typed.atom_types[j]},${typed.atom_types[k]},${typed.atom_types[l]}] ` +
          `class=${torsion_class(ctx, i, j, k, l)}`,
      );
    }
  }
}

console.log('\nproperty rows for the involved types:');
for (const t of [9, 32, 35, 41, 78, 79, 80, 81]) {
  console.log(`t${t}: ${JSON.stringify(ATOM_TYPE_PROPERTIES[t])}`);
}
