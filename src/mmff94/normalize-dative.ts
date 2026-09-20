/**
 * Dative drawings of hypervalent sulfur.
 *
 * C[S+2]([O-])([O-])C is a common way to write dimethyl sulfone without
 * explicit double bonds, and it is the same molecule as CS(=O)(=O)C. MMFF94
 * has types for the second form only: its S(IV) and S(VI) classes are
 * charge-neutral, and no MMFF94 type carries a positive charge on sulfur.
 * Typed exactly as drawn, that sulfur falls to the thiol/sulfide class and
 * the S=O parameters are never reached (measured: class 15 and a −1 e net
 * charge, 70.15 kcal/mol where the hypervalent drawing gives class 18 and
 * 4.90 — see docs/mmff94-compliance.md §5 item 9).
 *
 * normalize_dative_sulfur() rewrites those bonds before any typing decision
 * reads a bond order or a formal charge: each S(+n)–O(−) single bond becomes
 * an S=O double bond and both charges move one step toward neutral. Every
 * other atom, bond and coordinate is passed through untouched, and a
 * molecule that already carries the hypervalent drawing comes back
 * unchanged (the same object, not a copy).
 *
 * It is deliberately narrow, on one principle: normalize only where MMFF94
 * has no type for the charge-separated form but does have one for the
 * multiply-bonded form.
 *
 *   - Sulfur–oxygen pairs qualify (this file).
 *   - N+–O− N-oxides do not: MMFF94 has dedicated N-oxide types, and its
 *     reference entries use the charge-separated drawing.
 *   - P+–C− ylides do not: the ylide types are exactly that charge-separated
 *     form, and the WITTIG suite entry is drawn that way.
 *
 * Normalizing those would destroy correct typing, so they are left as drawn.
 */
import type { Molecule } from '../types.js';

const NEUTRAL = 0;

export function normalize_dative_sulfur(molecule: Molecule): Molecule {
  let changed = false;
  const atoms = molecule.atoms.map((a) => ({ ...a }));
  const bonds = molecule.bonds.map((b) => ({ ...b }));

  for (const bond of bonds) {
    if (bond.bond_order !== 1) continue;
    for (const [s, o] of [
      [bond.atom1, bond.atom2],
      [bond.atom2, bond.atom1],
    ] as const) {
      const sa = atoms[s];
      const oa = atoms[o];
      if (sa.element !== 'S' || oa.element !== 'O') continue;
      const sq = sa.formal_charge ?? NEUTRAL;
      const oq = oa.formal_charge ?? NEUTRAL;
      // The pair must be a drawn dative bond: sulfur positive, oxygen
      // negative. Anything else (a thiolate, a sulfate ester oxygen) is a
      // genuine charge-separated species and stays as drawn.
      if (sq < 1 || oq > -1) continue;
      bond.bond_order = 2;
      sa.formal_charge = sq - 1;
      oa.formal_charge = oq + 1;
      changed = true;
    }
  }

  return changed ? { ...molecule, atoms, bonds } : molecule;
}
