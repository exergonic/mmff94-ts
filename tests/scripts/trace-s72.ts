// Full per-atom charge decomposition for the S72-anion family.
import { readFileSync } from 'fs';
import { parse_mmd } from '../../src/utils/mmd-parser';
import { assign_atom_types } from '../../src/mmff94/atom-types';
import { assign_bci_charges } from '../../src/mmff94/charges';
import { BCI_PARAMS, BCI_DEFAULT_PARAMS } from '../../src/mmff94/parameters';
import { make_class_context, bond_type_flag } from '../../src/mmff94/parameters/parameter-classes';

const mmdText = readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8');

function reference_charges(name: string, nAtoms: number): number[] {
  const pchg = new Array(nAtoms).fill(0);
  let inMol = false;
  for (const line of mmdText.split('\n')) {
    const head = line.match(/^\s*\d+\s+\[(\w+),/);
    if (head) { inMol = head[1] === name; continue; }
    if (!inMol) continue;
    const p = line.trim().split(/\s+/);
    if (p.length >= 20) {
      const serial = parseInt(p[p.length - 1], 10);
      if (!isNaN(serial)) pchg[serial - 1] = parseFloat(p[p.length - 4]);
    }
  }
  return pchg;
}

for (const name of process.argv.slice(2)) {
  console.log(`\n===== ${name}`);
  const mol = parse_mmd(mmdText).find(m => m.name === name)!;
  const typed = assign_atom_types(mol);
  const ch = assign_bci_charges(typed);
  const ref = reference_charges(name, mol.atoms.length);
  const adj: number[][] = Array.from({ length: mol.atoms.length }, () => []);
  for (const bond of mol.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }
  const ctx = make_class_context(typed, adj);
  // per-atom BCI sum and q0-term, replicating charges.ts
  const bciSum = new Array(mol.atoms.length).fill(0);
  for (const bond of mol.bonds) {
    const ti = typed.atom_types[bond.atom1], tj = typed.atom_types[bond.atom2];
    const t_min = Math.min(ti, tj), t_max = Math.max(ti, tj);
    const cls = bond_type_flag(ctx, bond.atom1, bond.atom2);
    const entry = BCI_PARAMS[`${cls}-${t_min}-${t_max}`];
    let bci: number;
    if (entry) bci = entry.bci;
    else {
      bci = (BCI_DEFAULT_PARAMS[t_min]?.pbci ?? 0) - (BCI_DEFAULT_PARAMS[t_max]?.pbci ?? 0);
    }
    const flow = (a1: number, a2: number, s: number) => { bciSum[a1] += s; bciSum[a2] -= s; };
    if (entry) {
      if (ti === t_min) flow(bond.atom1, bond.atom2, -bci);
      else flow(bond.atom1, bond.atom2, bci);
    } else if (ti === t_min) flow(bond.atom1, bond.atom2, bci);
    else flow(bond.atom1, bond.atom2, -bci);
  }
  const q0Term = mol.atoms.map((_, i) => ch.partial_charges![i] - bciSum[i]);
  for (let i = 0; i < mol.atoms.length; i++) {
    const t = typed.atom_types[i];
    const dev = ch.partial_charges![i] - ref[i];
    if (Math.abs(dev) > 0.001 || t === 72 || t === 41 || t === 25) {
      console.log(
        `  a${i} ${mol.atoms[i].element} t${t} bci=${bciSum[i].toFixed(3)} q0term=${q0Term[i].toFixed(3)} ours=${ch.partial_charges![i].toFixed(3)} ref=${ref[i].toFixed(3)} dev=${dev.toFixed(3)} nbrs=${adj[i].map(n => `${mol.atoms[n].element}${typed.atom_types[n]}`).join(',')}`,
      );
    }
  }
}
