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

import type { TypedMolecule } from '../types';
import { ATOM_TYPES } from '../mmff94/parameters';

/**
 * Parse a complete .mmd file into an array of TypedMolecules.
 * Each molecule carries pre-assigned atom_types and partial_charges
 * from the OPTIMOL reference computation.
 */
export function parse_mmd(mmd_text: string): TypedMolecule[] {
  const molecules: TypedMolecule[] = [];
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

    // Track bonds as "a1-a2" → bond_order, to deduplicate while keeping order
    const bondMap = new Map<string, number>();

    const atoms: TypedMolecule['atoms'] = [];
    const atom_types: number[] = [];
    const partial_charges: number[] = [];

    for (let a = 0; a < numAtoms; a++) {
      if (i >= lines.length) break;
      const parts = lines[i++].trim().split(/\s+/);
      if (parts.length < 16) continue;

      const mmff_type = parseInt(parts[0]);

      // Read (neighbor, bond_order) pairs up to 6 pairs
      // Stop when we hit a field containing a decimal (the x coordinate)
      let cursor = 1;
      for (let pair = 0; pair < 6; pair++) {
        if (cursor >= parts.length - 9) break; // need 9 fields for xyz+metadata
        const nbr_raw = parts[cursor];
        // Coordinates have decimal points; neighbor indices don't
        if (nbr_raw.includes('.')) break;
        const nbr = parseInt(nbr_raw);
        const order = parseInt(parts[cursor + 1]);
        if (nbr > 0 && order > 0) {
          const a1 = a; // 0-based atom index
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

      // Partial charge
      const pchg = parseFloat(parts[cursor++]);
      partial_charges.push(pchg);

      // Derive element from MMFF94 type
      const typeDef = ATOM_TYPES[mmff_type];
      const element = typeDef ? typeDef.element : '?';

      atoms.push({ index: a, element, x, y, z });
      atom_types.push(mmff_type);
    }

    if (atoms.length === 0) continue;

    // Build bond list
    const bonds: TypedMolecule['bonds'] = [];
    for (const [key, order] of bondMap) {
      const [a1, a2] = key.split('-').map(Number);
      bonds.push({ atom1: a1, atom2: a2, bond_order: order });
    }

    molecules.push({
      name: code,
      atoms,
      bonds,
      atom_types,
      partial_charges,
    });
  }

  return molecules;
}
