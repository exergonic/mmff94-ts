/**
 * Torsion (dihedral) energy.
 *
 * Halgren1996, eq. (7):
 *
 *   E_tors = ½ · [V₁ · (1 + cos τ) + V₂ · (1 − cos 2τ) + V₃ · (1 + cos 3τ)]
 *
 * Written more compactly using our stored phases γ_n:
 *
 *   E_tors = Σ (V_n / 2) · [1 + cos(n · τ − γ_n)]   for n = 1, 2, 3
 *
 * where:
 *   V_n    = barrier height for the n-th term (kcal/mol)
 *   n      = periodicity (1, 2, or 3)
 *   τ      = current dihedral angle (degrees)
 *   γ_n    = phase shift for the n-th term (degrees)
 *
 * The two forms are equivalent when:
 *   γ₁ = 0°    →  1 + cos τ        (n=1)
 *   γ₂ = 180°  →  1 − cos 2τ       (n=2, since cos(2τ−180°) = −cos 2τ)
 *   γ₃ = 0°    →  1 + cos 3τ       (n=3)
 *
 * Convention: τ = 0° when i−j and k−l bonds are eclipsed (cis).
 * τ = 180° when anti (trans); staggered ethane is ±60° (gauche).
 * Sign follows IUPAC right-hand rule.
 *
 * Only single bonds (bond_order === 1) are evaluated for torsional
 * strain. Double and triple bonds are kept planar by angle bend and
 * out-of-plane terms, not by torsion.
 *
 * The parameter lookup tries type order (i, j, k, l) first, then the
 * reverse (l, k, j, i), to catch symmetric parameter definitions.
 */

import type { TypedMolecule } from '../../types';
import { TORSION_PARAMS, ATOM_PROPERTIES, type TorsionParams } from '../parameters';
import { dihedral_angle, Vec3 } from '../../utils/vector';
import { make_class_context, bond_type_flag, type ClassContext } from './bond-type';

// Torsion class TTijkl (MMFF part IV, p. 609):
//   1 = central bond BT-flagged (e.g. the central single bond of a
//       conjugated diene — V2 ≈ 1.8, not the alkene's 12)
//   2 = central bond single but a terminal bond BT-flagged
//   4 = all four atoms in the same 4-membered ring
//   5 = all four in a non-aromatic 5-ring with at least one sp3 C
//   0 = everything else (alkanes, alkenes, aromatic rings)
function torsion_class(ctx: ClassContext, i: number, j: number, k: number, l: number): number {
  const { mol } = ctx;
  const bt_ab = bond_type_flag(ctx, i, j);
  const bt_bc = bond_type_flag(ctx, j, k);
  const bt_cd = bond_type_flag(ctx, k, l);

  if (bt_bc === 1) return 1;

  // 4-ring: the closing bond i-l makes i-j-k-l a 4-cycle.
  if (ctx.adj[i].includes(l)) return 4;

  const central = mol.bonds.find(
    b => (b.atom1 === j && b.atom2 === k) || (b.atom1 === k && b.atom2 === j),
  );
  if (central && central.bond_order === 1) {
    if (bt_ab || bt_cd) return 2;
  }

  // 5-ring: a common neighbor x of i and l closes the cycle i-j-k-l-x.
  const types = [i, j, k, l].map(a => mol.atom_types[a]);
  if (types.includes(1) && types.every(t => !ATOM_PROPERTIES[t]?.arom)) {
    const x = ctx.adj[i].find(n => n !== j && n !== k && ctx.adj[l].includes(n));
    if (x !== undefined) return 5;
  }

  return 0;
}

// Torsion step-down chain (part I, p. 513): exact, then the asymmetric
// (EqLvl3, EqLvl5) terminal reductions. The par file stores each entry
// in ONE canonical direction; the order index decides which direction
// that is, and only that direction's chain is consulted — the other
// direction would hit wildcard defaults (e.g. '0-0-1-1-0') before the
// exact reversed entry.
function lookup_torsion(
  cls: number,
  ti: number,
  tj: number,
  tk: number,
  tl: number,
): TorsionParams | undefined {
  const lvl3 = (t: number) => ATOM_PROPERTIES[t]?.lvl3 ?? t;
  const lvl5 = (t: number) => ATOM_PROPERTIES[t]?.lvl5 ?? t;
  const chain = (a: number, b: number, c: number, d: number): TorsionParams | undefined => {
    const keys = [
      `${cls}-${a}-${b}-${c}-${d}`,
      `${cls}-${lvl3(a)}-${b}-${c}-${lvl5(d)}`,
      `${cls}-${lvl5(a)}-${b}-${c}-${lvl3(d)}`,
      `${cls}-${lvl5(a)}-${b}-${c}-${lvl5(d)}`,
    ];
    for (const key of keys) {
      const p = TORSION_PARAMS[key];
      if (p) return p;
    }
    return undefined;
  };

  const M = 136;
  const order =
    tk * M ** 3 + tj * M ** 2 + tl * M + ti - (tj * M ** 3 + tk * M ** 2 + ti * M + tl);
  return order >= 0 ? chain(ti, tj, tk, tl) : chain(tl, tk, tj, ti);
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

  // Iterate over all bonds (single and multiple — an alkene's C=C
  // torsion is real, V2 ≈ 12; only H-centered bonds are skipped)
  for (const bond of molecule.bonds) {
    const j = bond.atom1;
    const k = bond.atom2;

    const tj = molecule.atom_types[j];
    const tk = molecule.atom_types[k];

    const posJ: Vec3 = [molecule.atoms[j].x, molecule.atoms[j].y, molecule.atoms[j].z];
    const posK: Vec3 = [molecule.atoms[k].x, molecule.atoms[k].y, molecule.atoms[k].z];

    // Neighbors of j (excluding k) and neighbors of k (excluding j)
    const i_neighbors = adj[j].filter(n => n !== k);
    const l_neighbors = adj[k].filter(n => n !== j);

    if (i_neighbors.length === 0 || l_neighbors.length === 0) continue;

    for (const i of i_neighbors) {
      const ti = molecule.atom_types[i];
      const posI: Vec3 = [molecule.atoms[i].x, molecule.atoms[i].y, molecule.atoms[i].z];

      for (const l of l_neighbors) {
        const tl = molecule.atom_types[l];
        const posL: Vec3 = [molecule.atoms[l].x, molecule.atoms[l].y, molecule.atoms[l].z];

        // Class-scoped step-down lookup: TTijkl selects the class
        // (part IV p. 609 — conjugated central bond, terminal BT flags,
        // 4/5-rings), then the asymmetric EqLvl3/EqLvl5 chain runs in
        // both directions (part I p. 513). Exact types in either
        // direction always win — a wildcard default like '*-1-1-*'
        // must never preempt the exact reversed entry.
        const cls = torsion_class(ctx, i, j, k, l);
        const params = lookup_torsion(cls, ti, tj, tk, tl);

        if (!params) continue;

        // Compute dihedral angle in degrees
        const tau_rad = dihedral_angle(posI, posJ, posK, posL);
        const tau_deg = tau_rad * (180.0 / Math.PI);

        // Evaluate each Fourier term
        for (const term of params.terms) {
          const angle_rad = (term.periodicity * tau_deg - term.gamma) * (Math.PI / 180.0);
          total_energy += (term.V / 2.0) * (1.0 + Math.cos(angle_rad));
        }
      }
    }
  }

  return total_energy;
}
