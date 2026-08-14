// Probe: skipped (both-pilp) dihedrals in the ERULE fragments — what would
// each empirical reading contribute at the reference geometry?
import { readFileSync } from 'fs';
import { join } from 'path';
import { parse_mmd } from '../../src/utils/mmd-parser.js';
import { assign_atom_types } from '../../src/mmff94/assign-atom-types.js';
import { calc_energy } from '../../src/mmff94/energy/total.js';
import { parse_bmin_log } from '../../tests/scripts/bmin-log.js';
import { make_class_context, torsion_class, lookup_torsion, get_bond_order, is_aromatic_bond } from '../../src/mmff94/parameters/parameter-classes.js';
import { ATOM_TYPE_PROPERTIES } from '../../src/mmff94/parameters/index.js';
import { empirical_torsion } from '../../src/mmff94/parameters/empirical.js';
import { dihedral_angle } from '../../src/utils/vector.js';

const suiteDir = join(process.cwd(), 'tests', 'fixtures', 'validation-suite');
const molecules = parse_mmd(readFileSync(join(suiteDir, 'MMFF94.mmd'), 'utf-8'));
const refEnergies = parse_bmin_log(readFileSync(join(suiteDir, 'MMFF94_bmin.log'), 'utf-8'));

for (const code of ['ERULE_01', 'ERULE_02', 'ERULE_04', 'ERULE_08']) {
  const mol = molecules.find(m => m.name === code)!;
  const typed = assign_atom_types(mol);
  const got = calc_energy(typed);
  const ref = refEnergies.get(code)!;
  console.log(`\n${code}: our torsion=${got.torsion.toFixed(5)} ref=${ref.torsion.toFixed(5)} delta=${(got.torsion - ref.torsion).toFixed(5)}`);

  const adj: number[][] = Array.from({ length: typed.atoms.length }, () => []);
  for (const bond of typed.bonds) { adj[bond.atom1].push(bond.atom2); adj[bond.atom2].push(bond.atom1); }
  const ctx = make_class_context(typed, adj);

  for (const bond of typed.bonds) {
    const j = bond.atom1, k = bond.atom2;
    const i_neighbors = adj[j].filter(n => n !== k);
    const l_neighbors = adj[k].filter(n => n !== j);
    if (i_neighbors.length === 0 || l_neighbors.length === 0) continue;
    for (const i of i_neighbors) for (const l of l_neighbors) {
      if (l === i) continue;
      const cls = torsion_class(ctx, i, j, k, l);
      const table = lookup_torsion(cls, typed.atom_types[i], typed.atom_types[j], typed.atom_types[k], typed.atom_types[l]);
      if (table) continue; // only the skipped ones
      const pj = ATOM_TYPE_PROPERTIES[typed.atom_types[j]];
      const pk = ATOM_TYPE_PROPERTIES[typed.atom_types[k]];
      const ej = typed.atoms[j].element, ek = typed.atoms[k].element;
      const emp = empirical_torsion(pj, pk, ej, ek, get_bond_order(ctx, j, k), is_aromatic_bond(ctx, j, k));
      if (!emp.skip) continue;
      const tau = dihedral_angle(
        [typed.atoms[i].x, typed.atoms[i].y, typed.atoms[i].z],
        [typed.atoms[j].x, typed.atoms[j].y, typed.atoms[j].z],
        [typed.atoms[k].x, typed.atoms[k].y, typed.atoms[k].z],
        [typed.atoms[l].x, typed.atoms[l].y, typed.atoms[l].z]);
      const deg = tau * 180 / Math.PI;
      // candidate readings at this geometry:
      const cand = (V: number, n: number) => 0.5 * V * (1 + Math.cos(n * tau));
      // rule (h) V3: sqrt(Vb*Vc)/Nbc with Nbc=(crd-1)(crd-1)
      const vb = { C: 2.12, N: 1.5, O: 0.2, Si: 1.22, P: 2.4, S: 0.48 }[ej] ?? 0;
      const vc = { C: 2.12, N: 1.5, O: 0.2, Si: 1.22, P: 2.4, S: 0.48 }[ek] ?? 0;
      const ub = { C: 2.0, N: 2.0, O: 2.0, Si: 1.25, P: 1.25, S: 1.25 }[ej] ?? 0;
      const uc = { C: 2.0, N: 2.0, O: 2.0, Si: 1.25, P: 1.25, S: 1.25 }[ek] ?? 0;
      const nbc = ((pj?.crd ?? 2) - 1) * ((pk?.crd ?? 2) - 1);
      const V3_h = Math.sqrt(vb * vc) / nbc;
      const V2_c = 6 * 0.4 * Math.sqrt(ub * uc);           // rule (c) π=0.4 (order-2 reading)
      const V2_g = 6 * 0.15 * Math.sqrt(ub * uc);          // rule (g) case (4) π=0.15
      console.log(`  SKIPPED ${i}-${j}-${k}-${l} types ${typed.atom_types[i]}-${typed.atom_types[j]}-${typed.atom_types[k]}-${typed.atom_types[l]} tau=${deg.toFixed(1)}° ` +
        `rule(h) V3=${V3_h.toFixed(4)} E=${cand(V3_h, 3).toFixed(5)} | rule(c) V2=${V2_c.toFixed(4)} E=${cand(V2_c, 2).toFixed(5)} | rule(g) V2=${V2_g.toFixed(4)} E=${cand(V2_g, 2).toFixed(5)}`);
    }
  }
}
