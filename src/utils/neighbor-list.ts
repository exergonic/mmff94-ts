/**
 * Build a pair list for non-bonded interactions (vdW + electrostatic).
 *
 * For small molecules (< 200 atoms), a simple O(n²) double loop over
 * all atom pairs is fast enough. For larger molecules, we build a
 * neighbor list using a distance cutoff to skip pairs that are too
 * far apart to contribute (MMFF94's buffered 14-7 potential decays
 * rapidly beyond ~10 Å).
 *
 * In all cases, 1-2 (bonded) and 1-3 (angle) pairs are EXCLUDED.
 * 1-4 (torsion) pairs ARE included — they are calculated in full
 * here, with the 0.5/0.75 scaling applied by total.ts.
 */

import type { Molecule } from '../types';

export interface PairInfo {
  i: number;
  j: number;
  distance: number;
  is_1_4: boolean;   // true if atoms are exactly 3 bonds apart
}

/**
 * Generate the list of all non-bonded atom pairs that interact
 * via van der Waals and electrostatic terms.
 *
 * @param molecule The molecule (with at least atoms[] and bonds[]).
 * @param cutoff   Distance cutoff in Å. Pairs beyond this are skipped.
 * @returns Array of pair information.
 */
export function build_pair_list(
  molecule: Molecule,
  cutoff: number = 12.0
): PairInfo[] {
  // TODO: implement pair list generation.
  //
  // Steps:
  //   1. Build an adjacency list from molecule.bonds[].
  //   2. For each atom pair (i, j) where i < j:
  //      a. Skip if i == j.
  //      b. Skip if bonded (1-2): check adjacency list.
  //      c. Skip if 1-3 (share a common neighbor): walk adjacency.
  //      d. Classify as 1-4 or higher: check if exactly 3 bonds apart
  //         by BFS on the adjacency graph up to depth 3.
  //      e. Compute distance.
  //      f. Skip if distance > cutoff.
  //      g. Add to the pair list.
  //
  // Return the sorted (by i, then j) list of pairs.
  return [];
}

/**
 * For each atom, find the list of atoms that are exactly 3 bonds away
 * (1-4 pairs). Used for 1-4 scaling in total.ts.
 */
export function find_1_4_pairs(molecule: Molecule): Set<string> {
  // Use a BFS from each atom to depth 3. Atoms reached at exactly
  // depth 3 are 1-4 pairs.
  // Return a Set of "i-j" strings (with i < j) for O(1) lookup.
  return new Set();
}
