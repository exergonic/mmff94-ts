/**
 * Torsion (dihedral) energy.
 *
 * Halgren1996, eq. (7):
 *
 *   E_tors = ½ · V₁ · (1 + cos τ)
 *          + ½ · V₂ · (1 − cos 2τ)
 *          + ½ · V₃ · (1 + cos 3τ)
 *
 * where:
 *   V_n    = barrier height for the n-th term (kcal/mol)
 *   τ      = current dihedral angle (degrees)
 *
 * Convention: τ = 0° when i−j and k−l bonds are eclipsed (cis).
 * τ = 180° when anti (trans); staggered ethane is ±60° (gauche).
 * Sign follows IUPAC right-hand rule.
 *
 * Every bond is evaluated — an alkene's C=C torsion is real
 * (V₂ ≈ 12 kcal/mol, holding the alkene planar); only torsions with
 * no substituent on either central atom (or an H-centered central
 * bond) are skipped by construction.
 *
 * The parameters come from the class-scoped step-down lookup
 * (parameter-classes.ts: TTijkl selects the class, then the
 * asymmetric EqLvl3/EqLvl5 chain runs in the order-canonical
 * direction), falling to the part-IV empirical rules when the chain
 * misses entirely.
 */

import type { TypedMolecule } from '../../types.js';
import { dihedral_angle, Vec3 } from '../../utils/vector.js';
import { make_class_context, type ClassContext, torsion_class, lookup_torsion, get_bond_order, is_aromatic_bond } from '../parameters/parameter-classes.js';
import { ATOM_TYPE_PROPERTIES } from '../parameters/index.js';
import { empirical_torsion } from '../parameters/empirical.js';

/**
 * The three Fourier barrier heights for the dihedral i-j-k-l, resolved
 * exactly as the energy term resolves them: torsion class (TTijkl),
 * class-scoped step-down lookup, then the empirical rules (part V).
 *
 * Shared with the gradient so the two can never resolve different
 * parameters for the same dihedral.
 */
export interface TorsionTerms {
  v1: number;
  v2: number;
  v3: number;
}

/**
 * Resolve the torsion parameters for i-j-k-l.
 * Returns undefined when the dihedral is skipped (the empirical rule
 * says so — e.g. no substituent on either central atom).
 */
export function torsion_terms(
  ctx: ClassContext,
  molecule: TypedMolecule,
  i: number,
  j: number,
  k: number,
  l: number,
): TorsionTerms | undefined {
  const cls = torsion_class(ctx, i, j, k, l);
  const ti = molecule.atom_types[i];
  const tj = molecule.atom_types[j];
  const tk = molecule.atom_types[k];
  const tl = molecule.atom_types[l];
  const params = lookup_torsion(cls, ti, tj, tk, tl);
  if (params) {
    return { v1: params.v1, v2: params.v2, v3: params.v3 };
  }
  // The empirical rules take the j/k properties and elements plus the
  // graph facts they need (the j-k order and aromaticity) — the
  // pure-formula home in empirical.ts never touches the ClassContext.
  const emp = empirical_torsion(
    ATOM_TYPE_PROPERTIES[tj],
    ATOM_TYPE_PROPERTIES[tk],
    molecule.atoms[j].element,
    molecule.atoms[k].element,
    get_bond_order(ctx, j, k),
    is_aromatic_bond(ctx, j, k),
  );
  if (emp.skip) return undefined;
  return { v1: emp.v1, v2: emp.v2, v3: emp.v3 };
}

/**
 * Calculate the total torsional (dihedral) energy.
 */
export function calc_torsion_energy(molecule: TypedMolecule): number {
  let total_energy = 0.0;

  // Build adjacency list for neighbor lookups
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }
  const ctx = make_class_context(molecule, adj);

  // Each bond is the central j-k of its dihedrals i-j-k-l
  for (const bond of molecule.bonds) {
    const j = bond.atom1;
    const k = bond.atom2;

    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];
    const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

    // Substituents of j (excluding k) and of k (excluding j)
    const i_neighbors = adj[j].filter(n => n !== k);
    const l_neighbors = adj[k].filter(n => n !== j);

    if (i_neighbors.length === 0 || l_neighbors.length === 0) continue;

    for (const i of i_neighbors) {
      const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];

      for (const l of l_neighbors) {
        // A 3-membered ring closes on itself: i and l are the same
        // atom, so i-j-k-l is a triangle, not a dihedral. BatchMin
        // skips these (OBMol::FindTorsions: d == a); counting them
        // pins a spurious V3 term at τ = 0 (FUVDOP's triazine ring).
        if (l === i) continue;

        const posL: Vec3 = [molecule.atoms[l].x, molecule.atoms[l].y, molecule.atoms[l].z];

        // Class-scoped step-down lookup (TTijkl), with the part-IV
        // empirical rules as the final fallback — resolved by the
        // shared torsion_terms() helper (same for the gradient).
        const terms = torsion_terms(ctx, molecule, i, j, k, l);
        if (!terms) continue;

        // Halgren1996 eq. (7): the three Fourier terms, phases γ =
        // 0°/180°/0° folded in (cos(2τ−180°) = −cos 2τ)
        const tau_rad = dihedral_angle(posI, posJ, posK, posL);
        total_energy +=
          0.5 * terms.v1 * (1.0 + Math.cos(tau_rad)) +
          0.5 * terms.v2 * (1.0 - Math.cos(2.0 * tau_rad)) +
          0.5 * terms.v3 * (1.0 + Math.cos(3.0 * tau_rad));
      }
    }
  }

  return total_energy;
}
