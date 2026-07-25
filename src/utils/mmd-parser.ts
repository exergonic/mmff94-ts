/**
 * Parse an MMFF94 Validation Suite .mmd (MacroModel/BatchMin dat) file
 * into TypedMolecule objects with pre-assigned MMFF94 types and charges.
 *
 * Layout per molecule:
 *   Header:  N    [CODE,N,N,S,k] FF=mmff EM=default ...
 *   Atoms:   mmff_type  nbr1 bo1  nbr2 bo2  ...  x y z  label idx fchg pchg name serial
 *
 * Each atom line:
 *   1      MMFF94 type (integer)
 *   2-13   6 pairs of (neighbor_idx, bond_order), "0 0" = unused
 *   14-16  x, y, z coordinates
 *   17     label (e.g. "1XA")
 *   18     charge index (always 0)
 *   19     formal charge (float)
 *   20     partial charge (float)
 *   21+    atom name fragments + trailing serial number
 */

import type { Atom, Bond, Molecule } from '../types';

/**
 * Parse a complete .mmd file into an array of Molecules.
 *
 * The .mmd file stores coordinates, connectivity, formal charges, and
 * partial charges from the OPTIMOL reference computation. The first
 * numeric field per atom line is an internal BatchMin/OPTIMOL index,
 * NOT the MMFF94 atom type — so this function returns plain Molecules
 * without pre-assigned types. Atom typing is done separately via
 * assign_atom_types().
 */
export function parse_mmd(mmd_text: string): Molecule[] {
  const molecules: Molecule[] = [];
  const lines = mmd_text.split('\n');

  let i = 0;
  while (i < lines.length) {
    const line = lines[i];

    // Header: whitespace, atom count, then [CODE, ...
    const headerMatch = line.match(/^\s*(\d+)\s+\[(\w+),/);
    if (!headerMatch) { i++; continue; }

    const numAtoms = parseInt(headerMatch[1]);
    const code = headerMatch[2];
    i++;

    const bondMap = new Map<string, number>();

    const atoms: Atom[] = [];

    for (let a = 0; a < numAtoms; a++) {
      if (i >= lines.length) break;
      const parts = lines[i++].trim().split(/\s+/, 30);
      if (parts.length < 16) continue;

      // Skip the first field (internal BatchMin type index, not MMFF94 type)

      // Read (neighbor, bond_order) pairs up to 6 pairs
      let cursor = 1;
      for (let pair = 0; pair < 6; pair++) {
        if (cursor >= parts.length - 9) break;
        const nbr_raw = parts[cursor];
        if (nbr_raw.includes('.')) break;
        const nbr = parseInt(nbr_raw);
        const order = parseInt(parts[cursor + 1]);
        if (nbr > 0 && order > 0) {
          const a1 = a;
          const a2 = nbr - 1;
          const key = a1 < a2 ? `${a1}-${a2}` : `${a2}-${a1}`;
          bondMap.set(key, order);
        }
        cursor += 2;
      }

      const x = parseFloat(parts[cursor++]);
      const y = parseFloat(parts[cursor++]);
      const z = parseFloat(parts[cursor++]);

      // Skip label, charge index
      cursor += 2;

      // Formal charge (not stored, just consume)
      cursor++;

      // Partial charge (not stored, just consume)
      cursor++;

      // Extract element from the atom name (last several fields).
      // The atom name typically starts with the molecule prefix followed
      // by the element symbol and index, e.g. "FORM C1" → "C".
      const nameFields = parts.slice(cursor);
      let element = '?';
      for (const field of nameFields) {
        const match = field.match(/^([A-Z][a-z]?)\d*$/);
        if (match) { element = match[1]; break; }
      }

      atoms.push({ index: a, element, x, y, z });
    }

    if (atoms.length === 0) continue;

    // Build bond list
    const bonds: Bond[] = [];
    for (const [key, order] of bondMap) {
      const [a1, a2] = key.split('-').map(Number);
      bonds.push({ atom1: a1, atom2: a2, bond_order: order });
    }

    molecules.push({ name: code, atoms, bonds });
  }

  return molecules;
}
