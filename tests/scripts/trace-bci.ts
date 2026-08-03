// Trace BCI contributions for one atom of a molecule.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';
import { BCI_PARAMS, BCI_DEFAULT_PARAMS } from '../../src/mmff94/parameters';
import { make_class_context, bond_type_flag } from '../../src/mmff94/parameters/parameter-classes';

const name = process.argv[2] ?? 'CUBTUO';
const atomIdx = parseInt(process.argv[3] ?? '1', 10);
const text = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');
const mol = parse_mmd(text).find(m => m.name === name)!;
const typed = assign_atom_types(mol);

const adj: number[][] = Array.from({ length: mol.atoms.length }, () => []);
for (const bond of mol.bonds) {
  adj[bond.atom1].push(bond.atom2);
  adj[bond.atom2].push(bond.atom1);
}
const ctx = make_class_context(typed, adj);

let total = 0;
for (const bond of mol.bonds) {
  if (bond.atom1 !== atomIdx && bond.atom2 !== atomIdx) continue;
  const other = bond.atom1 === atomIdx ? bond.atom2 : bond.atom1;
  const ti = typed.atom_types[atomIdx];
  const tj = typed.atom_types[other];
  const t_min = Math.min(ti, tj);
  const t_max = Math.max(ti, tj);
  const cls = bond_type_flag(ctx, bond.atom1, bond.atom2);
  const entry = BCI_PARAMS[`${cls}-${t_min}-${t_max}`];
  let bci: number;
  if (entry) bci = entry.bci;
  else {
    const pa = BCI_DEFAULT_PARAMS[t_min]?.pbci ?? 0;
    const pb = BCI_DEFAULT_PARAMS[t_max]?.pbci ?? 0;
    bci = pa - pb;
  }
  const sign = typed.atom_types[atomIdx] === t_min ? -1 : 1;
  total += sign * bci;
  console.log(
    `bond ${atomIdx}-${other} (${mol.atoms[other].element}${bond.bond_order}) types ${ti}-${tj} ` +
    `cls ${cls} bci ${bci.toFixed(4)} contribution ${(sign * bci).toFixed(4)} ` +
    `${entry ? 'explicit' : 'default'}`,
  );
}
console.log(`total BCI for atom ${atomIdx}: ${total.toFixed(4)}`);
