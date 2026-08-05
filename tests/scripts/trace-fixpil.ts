// Trace FIXPIL's N81 charge.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types';
import { assign_bci_charges } from '../../src/mmff94/charges';

const mols = parse_mmd(readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8'));
const mol = mols.find(m => m.name === 'FIXPIL')!;
const typed = assign_atom_types(mol);
const ch = assign_bci_charges(typed);
typed.atoms.forEach((a, i) => {
  if (typed.atom_types[i] === 81) {
    const nbrs = mol.bonds.filter(b => b.atom1 === i || b.atom2 === i)
      .map(b => (b.atom1 === i ? b.atom2 : b.atom1));
    console.log(`atom ${i} type ${typed.atom_types[i]} fchg ${a.formal_charge} charge ${ch.partial_charges![i].toFixed(3)} nbrs ${nbrs.join(',')} nbrTypes ${nbrs.map(n => typed.atom_types[n]).join(',')}`);
  }
});
console.log('n atoms', typed.atoms.length);
