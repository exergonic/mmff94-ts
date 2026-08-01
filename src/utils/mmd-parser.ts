/**
 * Parse a MacroModel/BatchMin .mmd structure file into Molecule objects.
 *
 * Layout per molecule:
 *   Header:  N    [CODE,N,N,S,k] FF=mmff EM=default ...
 *   Atoms:   bmin_type  nbr1 bo1  nbr2 bo2  ...  x y z  label idx fchg pchg name serial
 *
 * Each atom line:
 *   1      internal BatchMin type index — NOT the canonical MMFF94 type, skipped
 *   2-13   6 pairs of (neighbor_idx, bond_order), "0 0" = unused
 *   14-16  x, y, z coordinates
 *   17     label (e.g. "1XA")
 *   18     charge index (always 0)
 *   19     formal charge (float)
 *   20     partial charge (float)
 *   21+    atom name fragments + trailing serial number
 */

import type { Molecule, Atom, Bond } from '../types';

// MacroModel atom type (the .mmd first column) → element symbol.
// Extracted from OpenBabel's data/types.txt MMD column — OpenBabel's
// mmd reader uses exactly this table, so elements agree with the
// typing reference by construction. Types not listed here (31, 51,
// 201+) are ones OpenBabel cannot translate either; their molecules
// never make it into the reference.
const MMD_ELEMENT: Record<number, string> = {
  1: 'C', 2: 'C', 3: 'C', 10: 'C', 11: 'C', 12: 'C',
  15: 'O', 16: 'O', 18: 'O', 20: 'O', 23: 'O',
  24: 'N', 25: 'N', 26: 'N', 29: 'N', 32: 'N', 35: 'N', 36: 'N',
  41: 'H', 42: 'H', 43: 'H', 44: 'H', 45: 'H', 48: 'H',
  49: 'S', 52: 'S',
  53: 'P', 54: 'B', 55: 'B', 56: 'F', 57: 'Cl', 58: 'Br', 59: 'I',
  60: 'Si',
};

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

      // Element from the MacroModel type index (first field): OpenBabel's
      // mmd reader maps it through data/types.txt (MMD → atomic number),
      // so this is the authoritative source — atom names carry residue
      // prefixes ("UNCH C1", "CE05 C1") and are only a fallback.
      const mmd_type = parseInt(parts[0], 10);
      let element = MMD_ELEMENT[mmd_type] ?? '?';

      if (element === '?') {
        // Fallback: derive from the atom name field ("C1", "Cl1"). The
        // strict pattern (one letter + optional lowercase) rejects
        // residue prefixes like "UNCH" or "CE05" so the scan reaches
        // the element field. Types outside the table (31, 51, 201+)
        // are ones OpenBabel cannot translate either — such molecules
        // never appear in the typing reference.
        const nameFields = parts.slice(cursor);
        for (const field of nameFields) {
          const match = field.match(/^([A-Z][a-z]?)\d*$/);
          if (match) { element = match[1]; break; }
        }
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
