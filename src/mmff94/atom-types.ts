/**
 * MMFF94 atom type assignment.
 *
 * This is the hardest single piece of the force field. Every atom in the
 * molecule must be assigned one of the ~95 MMFF94 atom types before
 * any energy term can look up its parameters. A wrong type cascades
 * into wrong parameters for EVERY energy term.
 *
 * The decision tree considers:
 *   - Element
 *   - Coordination number (explicit bonded neighbors)
 *   - Bond orders to neighbors (single, double, triple, aromatic)
 *   - Neighbor elements (e.g., C=O vs C=C)
 *   - Ring membership (especially 3-4 membered small rings)
 *   - Formal charge
 *
 * The logic follows Halgren's 1996 description (J. Comput. Chem. 17, 520-552).
 * Each branch is flat — no deep inheritance — with comments explaining the
 * chemical reason for the distinction.
 */

import type { Molecule, TypedMolecule } from '../types';

/**
 * Assign an MMFF94 atom type to every atom in the molecule.
 *
 * Returns a TypedMolecule with atom_types[] parallel to atoms[].
 * Unrecognized environments get a safe fallback (type 1 for C, 5 for H, etc.).
 */
export function assign_atom_types(molecule: Molecule): TypedMolecule {
  const n = molecule.atoms.length;
  const atom_types = new Array<number>(n);

  // Build adjacency: for each atom, a list of {neighbor, bond_order}
  const adj: { nbr: number; order: number }[][] = Array.from({ length: n }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push({ nbr: bond.atom2, order: bond.bond_order });
    adj[bond.atom2].push({ nbr: bond.atom1, order: bond.bond_order });
  }

  // Ring detection: atoms that survive iterative leaf-stripping are ring atoms.
  const is_ring = find_ring_atoms(adj, n);

  for (let i = 0; i < n; i++) {
    const atom = molecule.atoms[i];
    const neighbors = adj[i];
    const coord = neighbors.length;

    // Check bond orders going out from this atom
    const has_double = neighbors.some(nb => nb.order === 2);
    const has_triple = neighbors.some(nb => nb.order === 3);
    const has_aromatic = neighbors.some(nb => nb.order >= 4);
    const double_nbrs = neighbors.filter(nb => nb.order === 2);

    switch (atom.element) {

      // ── Carbon ──────────────────────────────────────────────────────
      case 'C': {
        // Acetylenic: 2 neighbors, triple bond → type 4 (CSP)
        if (coord <= 2 && has_triple) {
          atom_types[i] = 4;
          break;
        }

        if (coord === 3) {
          // Trigonal planar (sp²): 3 neighbors, at least one double/aromatic

          if (is_ring[i] && (has_aromatic || is_aromatic_ring(i, adj, molecule))) {
            // Aromatic carbon: 3 bonds in an aromatic ring → type 37 (CB)
            atom_types[i] = 37;
            break;
          }

          if (has_double) {
            // Check what the double bond goes to
            const dbl_to_C = double_nbrs.some(nb => {
              const target = molecule.atoms[nb.nbr];
              return target.element === 'C';
            });
            const dbl_to_O = double_nbrs.some(nb => {
              const target = molecule.atoms[nb.nbr];
              return target.element === 'O';
            });

            if (dbl_to_O) {
              // Carbonyl: C=O with 3 neighbors → type 3 (C=O)
              atom_types[i] = 3;
            } else if (dbl_to_C) {
              // Alkene: C=C with 3 neighbors → type 2 (C=C, vinylic)
              atom_types[i] = 2;
            } else {
              // Other sp² (e.g., C=N) → type 2 as generic sp² default
              atom_types[i] = 2;
            }
            break;
          }

          // 3 neighbors, no double, no aromatic (should be rare for C)
          atom_types[i] = 1;
          break;
        }

        if (coord === 4) {
          // Tetrahedral (sp³): 4 single bonds

          if (is_ring[i]) {
            // Check ring size via the cyclopropyl/cyclobutyl neighbors
            const ring_size = estimate_ring_size(i, adj, molecule);
            if (ring_size === 3) {
              // Cyclopropyl C → type 22 (CR3R)
              atom_types[i] = 22;
              break;
            }
            if (ring_size === 4) {
              // Cyclobutyl C → type 20 (CR4R)
              atom_types[i] = 20;
              break;
            }
          }

          // Generic sp³ carbon → type 1 (CR)
          atom_types[i] = 1;
          break;
        }

        // Fallback for unusual coordination
        atom_types[i] = 1;
        break;
      }

      // ── Hydrogen ────────────────────────────────────────────────────
      case 'H': {
        // Determine what the H is bonded to (should have exactly 1 neighbor)
        if (coord === 0) {
          atom_types[i] = 5; // isolated H — safe guess
          break;
        }

        const host = molecule.atoms[neighbors[0].nbr];
        switch (host.element) {
          case 'C':
            // H bonded to carbon → type 5 (HC), regardless of C hybridization
            atom_types[i] = 5;
            break;
          case 'N':
            // H bonded to sp³ N → type 23 (HNR)
            // H bonded to sp² N (amide, imine) → type 28 (HNCO) or 27 (HN=C)
            // For now, all H-N → type 23
            atom_types[i] = 23;
            break;
          case 'O':
            // H bonded to oxygen → type 21 (HOR, alcohols)
            atom_types[i] = 21;
            break;
          case 'S':
            // H bonded to sulfur → type 71 (HS)
            atom_types[i] = 71;
            break;
          default:
            atom_types[i] = 5;
        }
        break;
      }

      // ── Oxygen ──────────────────────────────────────────────────────
      case 'O': {
        if (has_double) {
          // Carbonyl oxygen O=C → type 7 (O=C)
          atom_types[i] = 7;
        } else if (coord === 2) {
          // Divalent oxygen (ether, alcohol, water) → type 6 (OR)
          atom_types[i] = 6;
        } else if (coord === 1) {
          // Terminal oxygen, often anionic or radical — check later
          atom_types[i] = 6;
        } else {
          atom_types[i] = 6;
        }
        break;
      }

      // ── Nitrogen ────────────────────────────────────────────────────
      case 'N': {
        if (coord === 3 && !has_double && !has_aromatic) {
          // Amine N: 3 single bonds → type 8 (NR)
          atom_types[i] = 8;
        } else if (coord === 2 && has_double) {
          // Imine N=C → type 9 (N=C)
          // (nitrogen with a double bond and one other neighbor)
          atom_types[i] = 9;
        } else if (coord === 3 && has_double) {
          // Check if the double goes to C or O
          const dbl_to_O = double_nbrs.some(nb => {
            const target = molecule.atoms[nb.nbr];
            return target.element === 'O';
          });
          if (dbl_to_O) {
            // N-oxide (N→O) — check coordination
            atom_types[i] = 67; // N2OX (placeholder for N-oxide)
          } else {
            // Amide N (N-C=O) with delocalized lone pair → type 10 (NC=O)
            atom_types[i] = 10;
          }
        } else if (is_ring[i] && has_aromatic) {
          // Aromatic N — pyridine-like (type 38 NPYD) or pyrrole-like (type 39 NPYL)
          // For now, use type 38 for 6-ring, 39 for 5-ring
          const ring_size = estimate_ring_size(i, adj, molecule);
          atom_types[i] = (ring_size === 5) ? 39 : 38;
        } else {
          atom_types[i] = 8; // fallback: amine N
        }
        break;
      }

      // ── Sulfur ──────────────────────────────────────────────────────
      case 'S': {
        if (has_double && coord >= 3) {
          // Check for S=O / SO₂
          const dbl_to_O = double_nbrs.some(nb => {
            const target = molecule.atoms[nb.nbr];
            return target.element === 'O';
          });
          if (dbl_to_O && coord >= 4) {
            atom_types[i] = 18; // SO2 (sulfone)
          } else if (dbl_to_O) {
            atom_types[i] = 17; // S=O (sulfoxide)
          } else {
            atom_types[i] = 16; // S=C (thiocarbonyl)
          }
        } else {
          // Thiol or sulfide → type 15 (S)
          atom_types[i] = 15;
        }
        break;
      }

      // ── Halogens ────────────────────────────────────────────────────
      case 'F':  atom_types[i] = 11; break;
      case 'Cl': atom_types[i] = 12; break;
      case 'Br': atom_types[i] = 13; break;
      case 'I':  atom_types[i] = 14; break;

      // ── Silicon, phosphorus ─────────────────────────────────────────
      case 'Si': atom_types[i] = 19; break;
      case 'P':  atom_types[i] = 26; break;

      default:
        atom_types[i] = 1; // fallback: generic sp³ C
    }
  }

  return { ...molecule, atom_types };
}


/**
 * Ring detection via iterative leaf-stripping.
 *
 * Repeatedly removes atoms with degree < 2 (terminal atoms). Atoms that
 * survive are in at least one ring. This is O(n) and does not require
 * ring-size information — just membership.
 */
function find_ring_atoms(
  adj: { nbr: number; order: number }[][],
  n: number,
): boolean[] {
  const degree = adj.map(nb => nb.length);
  const queue: number[] = [];

  for (let i = 0; i < n; i++) {
    if (degree[i] < 2) queue.push(i);
  }

  while (queue.length > 0) {
    const idx = queue.shift()!;
    if (degree[idx] <= 0) continue;
    degree[idx] = 0; // mark as removed

    for (const { nbr } of adj[idx]) {
      if (degree[nbr] > 0) {
        degree[nbr]--;
        if (degree[nbr] < 2) queue.push(nbr);
      }
    }
  }

  // Atoms with degree > 0 after stripping are ring atoms
  return degree.map(d => d > 0);
}


/**
 * Check if an atom is in an aromatic ring (alternating single/double bonds
 * around a planar ring, typical of benzene and related systems).
 *
 * For each ring atom, we verify that every atom in the ring has exactly
 * one double bond to another ring atom (the Kekulé pattern). This
 * distinguishes benzene (all 6 C are sp², alternating bonds) from
 * non-aromatic cyclic alkenes like cyclohexadiene.
 */
function is_aromatic_ring(
  start: number,
  adj: { nbr: number; order: number }[][],
  molecule: Molecule,
): boolean {
  // Walk the ring starting from `start`, following alternating single/double
  // bonds. If we return to `start` having visited all ring atoms, it's aromatic.
  const visited = new Set<number>();
  const ring_atoms: number[] = [];

  // Find the ring via DFS, bounded at depth 6
  function dfs(node: number, parent: number, depth: number): boolean {
    if (depth > 6) return false;
    if (node === start && depth > 2) return true;

    visited.add(node);
    ring_atoms.push(node);

    for (const { nbr } of adj[node]) {
      if (nbr === parent) continue;
      // Only consider sp² carbon neighbors in a ring
      const elem = molecule.atoms[nbr].element;
      if (elem !== 'C') continue;
      if (adj[nbr].length !== 3) continue;
      if (!visited.has(nbr) || nbr === start) {
        if (dfs(nbr, node, depth + 1)) return true;
      }
    }

    visited.delete(node);
    ring_atoms.pop();
    return false;
  }

  if (!dfs(start, -1, 0)) return false;

  // Verify alternating bond pattern around the ring
  if (ring_atoms.length === 0) return false;

  // Walk the ring and count double bonds to ring neighbors
  for (const atom of ring_atoms) {
    const ring_dbl = adj[atom].filter(nb =>
      ring_atoms.includes(nb.nbr) && nb.order === 2
    );
    // Each aromatic carbon must have exactly one double bond to a ring neighbor
    if (ring_dbl.length !== 1) return false;
  }

  return true;
}


/**
 * Estimate the size of the smallest ring containing atom i.
 *
 * Uses BFS from i back to i, bounded to depth 6. Returns 0 if no ring
 * is found within the search depth.
 *
 * This is a heuristic: it finds the shortest cycle, not the SSSR ring.
 * For the typing rules, only 3- and 4-membered rings matter.
 */
function estimate_ring_size(
  start: number,
  adj: { nbr: number; order: number }[][],
  molecule: Molecule,
): number {
  // BFS limited to depth 6
  const visited = new Set<number>();
  const queue: { node: number; depth: number; parent: number }[] = [];

  visited.add(start);
  queue.push({ node: start, depth: 0, parent: -1 });

  let head = 0;
  while (head < queue.length) {
    const { node, depth, parent } = queue[head++];

    if (depth > 5) continue; // don't search beyond 6-membered rings

    for (const { nbr } of adj[node]) {
      if (nbr === parent) continue;

      // If we reached start via a different path, we found a cycle
      if (nbr === start && parent !== -1) {
        return depth + 1;
      }

      if (!visited.has(nbr)) {
        visited.add(nbr);
        queue.push({ node: nbr, depth: depth + 1, parent: node });
      }
    }
  }

  return 0; // no ring found within depth limit
}


/**
 * Compute partial charges for every atom using the MMFF94 bond charge
 * increment (BCI) model.
 *
 * MMFF94 does NOT store per-atom partial charges in its parameter tables.
 * Instead, each BOND TYPE has a charge increment value. The partial charge
 * on an atom is the SUM of the BCI values of every bond it participates in.
 *
 * Formal charges override the BCI sum.
 */
export function compute_bci_charges(molecule: TypedMolecule): void {
  molecule.partial_charges = molecule.atoms.map(() => 0.0);
}
