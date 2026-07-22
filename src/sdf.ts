/**
 * Parse an SDF / MOL V2000 block into our Molecule data model.
 *
 * This is a standalone parser — not borrowed from Valence or any other project.
 * It reads the standard V2000 format that OpenBabel, RDKit, and PubChem all emit.
 *
 * V2000 layout (fixed-width, 3-line header, then atom and bond blocks):
 *   Lines 1-3:   header (name, program, comment — we ignore most of this)
 *   Line 4:      counts line (number of atoms, number of bonds, ...)
 *   Atom block:  1 line per atom, columns 0-9 = x, 10-19 = y, 20-29 = z, 31-33 = element
 *   Bond block:  1 line per bond, columns 0-2 = atom1, 3-5 = atom2, 6-8 = bond order
 *   "M END"      terminates the connection table
 */

import type { Atom, Bond, Molecule } from './types';

/**
 * Parse a V2000 MOL block or SDF record into a Molecule.
 *
 * Returns an empty Molecule ({ atoms: [], bonds: [] }) if the input
 * cannot be parsed, so callers can fall through gracefully.
 */
export function parse_sdf(sdf_text: string): Molecule {
  const lines = sdf_text.split('\n');
  if (lines.length < 4) return { atoms: [], bonds: [] };

  // Line 4 (index 3) is the counts line.
  // Columns 0-2 = number of atoms (right-justified, 3 chars).
  // Columns 3-5 = number of bonds (right-justified, 3 chars).
  const counts_line = lines[3];

  let num_atoms = 0;
  let num_bonds = 0;
  try {
    num_atoms = parseInt(counts_line.substring(0, 3).trim(), 10) || 0;
    num_bonds = parseInt(counts_line.substring(3, 6).trim(), 10) || 0;
  } catch {
    return { atoms: [], bonds: [] };
  }

  // Sanity: a molecule with 0 atoms or > 999 atoms is probably a parse error.
  if (num_atoms === 0 || num_atoms > 999) return { atoms: [], bonds: [] };

  // Parse atom block: lines 4 through 4 + num_atoms - 1
  const atoms: Atom[] = [];
  const atom_start = 4;
  for (let i = 0; i < num_atoms; i++) {
    const line = lines[atom_start + i];
    if (!line || line.length < 34) continue;

    const x = parse_float(line, 0, 10);
    const y = parse_float(line, 10, 20);
    const z = parse_float(line, 20, 30);
    const element = line.substring(31, 34).trim();

    if (!element) continue;      // No element symbol — malformed line, skip
    if (element === 'H' || element === 'C' || element === 'N' || element === 'O' ||
        element === 'S' || element === 'P' || element === 'F' || element === 'Cl' ||
        element === 'Br' || element === 'I' || element === 'Si' || element === 'B' ||
        element === 'Fe' || element === 'Cu' || element === 'Zn' || element === 'Mn') {
      // Known element — proceed
    } else if (element.length === 1 && element >= 'A' && element <= 'Z') {
      // Single uppercase letter that isn't in our list — still a valid element
    } else if (element.length === 2 && element[0] >= 'A' && element[0] <= 'Z' &&
               element[1] >= 'a' && element[1] <= 'z') {
      // Two-character element symbol — still valid
    } else {
      // Not a recognized element symbol — assign as '?' (atoms that pass through this
      // parser are expected to be known, but we treat unrecognized ones as inert)
      continue;
    }

    atoms.push({
      index: i,
      element,
      x: isNaN(x) ? 0 : x,
      y: isNaN(y) ? 0 : y,
      z: isNaN(z) ? 0 : z,
    });
  }

  // Parse bond block: starts after the atom block.
  const bond_start = atom_start + num_atoms;
  const bonds: Bond[] = [];
  for (let i = 0; i < num_bonds; i++) {
    const line = lines[bond_start + i];
    if (!line || line.length < 9) continue;

    const a1 = parseInt(line.substring(0, 3).trim(), 10) - 1;  // V2000 is 1-indexed
    const a2 = parseInt(line.substring(3, 6).trim(), 10) - 1;
    const order = parseInt(line.substring(6, 9).trim(), 10) || 1;

    // Skip bonds that reference atoms outside our parsed range
    if (a1 < 0 || a1 >= atoms.length || a2 < 0 || a2 >= atoms.length) continue;

    bonds.push({ atom1: a1, atom2: a2, bond_order: order });
  }

  // First line of the SDF is typically the molecule name
  const name = lines[0]?.trim() || undefined;

  return { atoms, bonds, name };
}

/** Parse a fixed-width float field from a string. */
function parse_float(line: string, start: number, end: number): number {
  const fragment = line.substring(start, end).trim();
  if (fragment === '') return NaN;
  return parseFloat(fragment);
}
