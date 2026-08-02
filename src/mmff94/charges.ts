/**
 * The MMFF94 partial-charge model: bond charge increments (BCI).
 *
 * MMFF94 does NOT store per-atom partial charges in its parameter
 * tables. Instead, each BOND TYPE pair has a charge increment (the
 * bci value, mmffchg.par), and the partial charge on an atom is the
 * sum of the increments of every bond it participates in:
 *
 *   q_i = Σ_j bci(i, j)
 *
 * The par stores each pair with its types in ascending order, and the
 * increment is SUBTRACTED from the smaller-type atom and ADDED to the
 * larger-type atom (Halgren part III). So in ammonia, the N(8)–H(23)
 * bond carries bci = +0.36: the N collects −0.36 from each of its
 * three H's (−1.08 total) and each H collects +0.36.
 *
 * Bonds whose type pair has no entry use the per-atom default values
 * (mmffpbci.par): the pair contributes P_i − P_j, the difference of
 * the two atoms' defaults — the same flow from the smaller type.
 *
 * The lookup is class-scoped like every other term: a conjugated
 * single bond (BTij = 1) uses the class-1 entry when one exists, and
 * falls to the per-atom defaults when it does not (it does NOT fall
 * to the class-0 entry — the reference is strict about this).
 *
 * Formal charges override the BCI sum in the full model (the fcadj
 * factors in mmffpbci.par). Our parser does not read formal charges
 * yet, so the override is not applied — the fixtures and the
 * neutral suite molecules are unaffected.
 */

import type { TypedMolecule } from '../types';
import { BCI_PARAMS, BCI_DEFAULT_PARAMS } from './parameters';
import { make_class_context, bond_type_flag } from './parameters/parameter-classes';

/**
 * Compute the partial charges for every atom of a typed molecule and
 * store them on `molecule.partial_charges` (the electrostatics term
 * and the future gradient both consume them).
 */
export function compute_bci_charges(molecule: TypedMolecule): void {
  // Adjacency + the shared class context (the BTij flag selects the
  // bci class for conjugated single bonds).
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }
  const ctx = make_class_context(molecule, adj);

  const charges = molecule.atoms.map(() => 0.0);
  for (const bond of molecule.bonds) {
    const ti = molecule.atom_types[bond.atom1];
    const tj = molecule.atom_types[bond.atom2];
    const t_min = Math.min(ti, tj);
    const t_max = Math.max(ti, tj);
    const cls = bond_type_flag(ctx, bond.atom1, bond.atom2);
    const entry = BCI_PARAMS[`${cls}-${t_min}-${t_max}`];

    let bci: number;
    if (entry) {
      bci = entry.bci;
    } else {
      // Unparametrized pair: the difference of the per-atom defaults.
      const pa = BCI_DEFAULT_PARAMS[t_min]?.pbci ?? 0;
      const pb = BCI_DEFAULT_PARAMS[t_max]?.pbci ?? 0;
      bci = pa - pb;
    }

    // The increment flows from the smaller type to the larger one.
    if (ti === t_min) {
      charges[bond.atom1] -= bci;
      charges[bond.atom2] += bci;
    } else {
      charges[bond.atom1] += bci;
      charges[bond.atom2] -= bci;
    }
  }

  molecule.partial_charges = charges;
}
