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

  // Aromatic ring perception: 5/6-membered rings with 6 π electrons
  // (Kekulé doubles + lone-pair heteroatoms), or V2000 aromatic bonds
  // (order 4). The per-atom rings feed the 5-ring α/β typing below.
  const aromatic = find_aromatic_rings(adj, molecule, is_ring);

  // Water detection: MMFF94 gives H2O dedicated types — O = 70
  // ("OXYGEN IN WATER"), H = 31 ("H-OH") — with r₀ 0.969 (bond
  // 0-31-70), distinct from alcohols (O 6, H 21, r₀ 0.972). An O
  // with exactly two H neighbors, each H bonded only to that O, is
  // water. Pre-scanned so the H case works regardless of atom order.
  const water_oxygens = new Set<number>();
  for (let i = 0; i < n; i++) {
    const atom = molecule.atoms[i];
    if (atom.element !== 'O' || adj[i].length !== 2) continue;
    const hs = adj[i].map(nb => nb.nbr);
    if (!hs.every(h => molecule.atoms[h].element === 'H' && adj[h].length === 1)) continue;
    water_oxygens.add(i);
  }

  // Carboxylate detection: C(=O)-O⁻ with a TERMINAL single-bonded oxygen
  // (no H or other bonds — the free acid's -OH has two neighbors and
  // stays 3/7/6). Carboxylate C → type 41 (CO₂M), both oxygens → 32
  // (O2CM). Pre-scanned so the O case works regardless of atom order.
  const carboxylate_carbons = new Set<number>();
  for (let i = 0; i < n; i++) {
    const atom = molecule.atoms[i];
    if (atom.element !== 'C') continue;
    const has_dbl_O = adj[i].some(nb => nb.order === 2 && molecule.atoms[nb.nbr].element === 'O');
    const has_terminal_O = adj[i].some(
      nb => nb.order === 1 && molecule.atoms[nb.nbr].element === 'O' && adj[nb.nbr].length === 1,
    );
    if (has_dbl_O && has_terminal_O) carboxylate_carbons.add(i);
  }

  // Amide N detection: N bonded to a carbonyl (or thiocarbonyl) carbon
  // — C=O/C=S, non-aromatic. The N itself holds only single bonds
  // (formamide's N: N–C(=O)H); its trigonal character is resonance, so
  // the test must inspect the NEIGHBOR's double bond, not this atom's.
  // Pre-scanned so both the N (type 10, NC=O) and its hydrogens (type
  // 28, HNCO) type correctly regardless of atom order.
  const amide_nitrogens = new Set<number>();
  for (let i = 0; i < n; i++) {
    const atom = molecule.atoms[i];
    if (atom.element !== 'N') continue;
    if (adj[i].some(nb => nb.order >= 2)) continue; // an N=double is imine/N-oxide, not amide
    // Sulfonamide N (NSO₂, type 43 — a roadmap item) takes precedence
    // over amide: an N on a sulfonyl S is never NC=O. Without this,
    // DIKGAF's N (bonded to both SO₂ and an ester carbonyl C) was
    // mis-typed 10, dragging a −1.6 kcal/mol oop term (the negative
    // amide k_oop) onto a center that has none.
    const bonded_to_sulfonyl = adj[i].some(nb => {
      const nbr = molecule.atoms[nb.nbr];
      if (nbr.element !== 'S' && nbr.element !== 'P') return false;
      let terminal_oxygens = 0;
      for (const b of adj[nb.nbr]) {
        if (molecule.atoms[b.nbr].element === 'O' && adj[b.nbr].length === 1) terminal_oxygens++;
      }
      return terminal_oxygens >= 2;
    });
    if (bonded_to_sulfonyl) continue;
    const bonded_to_carbonyl = adj[i].some(nb => {
      if (molecule.atoms[nb.nbr].element !== 'C') return false;
      return adj[nb.nbr].some(
        b => !(b.order >= 4) && b.order === 2 && ['O', 'S'].includes(molecule.atoms[b.nbr].element),
      );
    });
    if (bonded_to_carbonyl) amide_nitrogens.add(i);
  }

  for (let i = 0; i < n; i++) {
    const atom = molecule.atoms[i];
    const neighbors = adj[i];
    const n_neighbors = neighbors.length;

    // Check bond orders going out from this atom
    const has_double = neighbors.some(nb => nb.order === 2);
    const has_triple = neighbors.some(nb => nb.order === 3);
    // V2000 stores aromatic bonds as order 4; the SDF parser passes them through
    const has_aromatic = neighbors.some(nb => nb.order >= 4);
    const double_nbrs = neighbors.filter(nb => nb.order === 2);

    switch (atom.element) {

      // ── Carbon ──────────────────────────────────────────────────────
      case 'C': {
        // Aromatic ring carbon: 37 in a 6-ring (benzene/pyridine C),
        // 63/64/78 in a 5-ring by the α/β position relative to the
        // lone-pair heteroatoms (part II).
        if (aromatic.atoms.has(i)) {
          const in_5ring = (aromatic.rings_of.get(i) ?? []).some(r => r.path.length === 5);
          atom_types[i] = in_5ring
            ? type_aromatic_5ring_carbon(i, adj, molecule, aromatic)
            : 37;
          break;
        }

        // Acetylenic: 2 neighbors, triple bond → type 4 (CSP)
        if (n_neighbors <= 2 && has_triple) {
          atom_types[i] = 4;
          break;
        }

        if (n_neighbors === 3) {
          // Trigonal planar (sp²): 3 neighbors, at least one double/aromatic

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

            if (carboxylate_carbons.has(i)) {
              // Carboxylate carbon C(=O)O⁻ → type 41 (CO₂M)
              atom_types[i] = 41;
            } else if (dbl_to_O) {
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

        if (n_neighbors === 4) {
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

        // Fallback for unusual coordination numbers
        atom_types[i] = 1;
        break;
      }

      // ── Hydrogen ────────────────────────────────────────────────────
      case 'H': {
        // Determine what the H is bonded to (should have exactly 1 neighbor)
        if (n_neighbors === 0) {
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
            // H on an amide N (type 10) → 28 (HNCO). Other N–H stays
            // type 23 (HNR) for now; the imine-H (27, HN=C) and
            // sulfonamide-H (72) are roadmap items.
            atom_types[i] = amide_nitrogens.has(neighbors[0].nbr) ? 28 : 23;
            break;
          case 'O':
            // Water hydrogen → type 31 (H-OH); other O-H → type 21 (HOR)
            atom_types[i] = water_oxygens.has(neighbors[0].nbr) ? 31 : 21;
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
        if (water_oxygens.has(i)) {
          // Water oxygen → type 70 (OXYGEN IN WATER)
          atom_types[i] = 70;
        } else if ((aromatic.rings_of.get(i) ?? []).some(r => r.path.length === 5)) {
          // Aromatic 5-ring oxygen (furan, oxazole) → type 59 (OFUR)
          atom_types[i] = 59;
        } else if (n_neighbors === 1 && carboxylate_carbons.has(neighbors[0].nbr)) {
          // Carboxylate oxygen — the =O or the terminal -O⁻ → type 32 (O2CM)
          atom_types[i] = 32;
        } else if (has_double) {
          // Carbonyl oxygen O=C → type 7 (O=C)
          atom_types[i] = 7;
        } else {
          // Ether, alcohol, and terminal O are all type 6 for now;
          // the carboxylate (32) and oxide (35) distinctions are not yet typed.
          atom_types[i] = 6;
        }
        break;
      }

      // ── Nitrogen ────────────────────────────────────────────────────
      case 'N': {
        // Aromatic ring nitrogen: 38 in a 6-ring (pyridine), 39/65/66
        // in a 5-ring by position (pyrrole's N is 39; imidazole- and
        // thiazole-type N's are 65/66). The charged variants (58/76/81)
        // follow the same positions (part II).
        if (aromatic.atoms.has(i)) {
          const in_5ring = (aromatic.rings_of.get(i) ?? []).some(r => r.path.length === 5);
          atom_types[i] = in_5ring
            ? type_aromatic_5ring_nitrogen(i, adj, molecule, aromatic)
            : adj[i].length === 3 ? 58 : 38;
          break;
        }

        // Amide N (type 10, NC=O): single-bonded to a carbonyl carbon —
        // trigonal by resonance, not by an N=double bond (formamide).
        // See the amide_nitrogens pre-scan above.
        if (amide_nitrogens.has(i)) {
          atom_types[i] = 10;
          break;
        }

        if (n_neighbors === 3 && !has_double && !has_aromatic) {
          // Amine N: 3 single bonds → type 8 (NR)
          atom_types[i] = 8;
        } else if (n_neighbors === 2 && has_double) {
          // Imine N=C → type 9 (N=C)
          // (nitrogen with a double bond and one other neighbor)
          atom_types[i] = 9;
        } else if (n_neighbors === 3 && has_double) {
          // N with three neighbors AND its own double bond: N=O is the
          // N-oxide (type 67, placeholder). N=C/N=N (enamine, N-N=C
          // with delocalized lone pair) also take type 10 in the
          // reference — the true amide (N–C(=O), no N double bond) is
          // handled above by the amide_nitrogens pre-scan.
          const dbl_to_O = double_nbrs.some(nb => {
            const target = molecule.atoms[nb.nbr];
            return target.element === 'O';
          });
          if (dbl_to_O) {
            // N-oxide (N→O): placeholder — N2OX (type 67)
            atom_types[i] = 67;
          } else {
            atom_types[i] = 10;
          }
        } else {
          atom_types[i] = 8; // fallback: amine N
        }
        break;
      }

      // ── Sulfur ──────────────────────────────────────────────────────
      case 'S': {
        // Aromatic 5-ring sulfur (thiophene, thiazole) → type 44 (STHI);
        // 6-ring aromatic S falls through to the thioether typing (MMFF
        // has no aromatic type for it).
        const in_aromatic_5ring = (aromatic.rings_of.get(i) ?? []).some(r => r.path.length === 5);
        if (aromatic.atoms.has(i) && in_aromatic_5ring) {
          atom_types[i] = 44;
          break;
        }
        if (has_double && n_neighbors >= 3) {
          // Check for S=O / SO₂
          const dbl_to_O = double_nbrs.some(nb => {
            const target = molecule.atoms[nb.nbr];
            return target.element === 'O';
          });
          if (dbl_to_O && n_neighbors >= 4) {
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
      // All P → 26 for now (placeholder; phosphate P, type 25, not yet handled)
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
 * Aromatic ring perception (MMFF part II): planar 5- or 6-membered
 * rings with 6 π electrons — each ring double bond contributes 2, and
 * each ring N/O/S carrying NO ring double bond contributes its lone
 * pair (2): pyrrole's N, furan's O, thiophene's S. V2000 aromatic
 * bonds (order 4) settle the ring directly — the input already
 * resolved it. Replaces the old carbon-only alternating-bond walker,
 * which could not perceive heteroaromatic rings (pyridine typed
 * N=9, pyrrole's ring C's typed 2).
 *
 * Returns the aromatic atoms and, per atom, the aromatic rings
 * containing it (fused systems give more than one).
 */
interface AromaticRing {
  path: number[];
}

interface AromaticInfo {
  atoms: Set<number>;
  rings_of: Map<number, AromaticRing[]>;
}

function find_aromatic_rings(
  adj: { nbr: number; order: number }[][],
  molecule: Molecule,
  is_ring: boolean[],
): AromaticInfo {
  const atoms = new Set<number>();
  const rings_of = new Map<number, AromaticRing[]>();
  const seen = new Set<string>();
  const bond_order = (a: number, b: number): number => {
    const bond = molecule.bonds.find(
      bd => (bd.atom1 === a && bd.atom2 === b) || (bd.atom1 === b && bd.atom2 === a),
    );
    return bond ? bond.bond_order : 0;
  };

  // All simple cycles through ring atoms, bounded at 6-membered rings.
  const cycles: number[][] = [];
  for (let start = 0; start < molecule.atoms.length; start++) {
    if (!is_ring[start]) continue;
    const walk = (node: number, parent: number, path: number[]): void => {
      for (const { nbr } of adj[node]) {
        if (nbr === parent || !is_ring[nbr]) continue;
        if (nbr === start) {
          if (path.length >= 4 && path.length <= 6) {
            const key = [...path].sort((a, b) => a - b).join(',');
            if (!seen.has(key)) {
              seen.add(key);
              cycles.push([...path]);
            }
          }
        } else if (path.length < 6 && !path.includes(nbr)) {
          walk(nbr, node, [...path, nbr]);
        }
      }
    };
    walk(start, -1, [start]);
  }

  const mark = (path: number[]): void => {
    const ring = { path };
    for (const a of path) {
      atoms.add(a);
      const list = rings_of.get(a);
      if (list) list.push(ring);
      else rings_of.set(a, [ring]);
    }
  };

  for (const path of cycles) {
    const set = new Set(path);
    const size = path.length;

    // Chord check: a genuine ring has no bond between non-consecutive
    // members. Cage molecules (e.g. FUVDOP's triazine) are full of
    // chords; without this, the cycle walker could close a 6-cycle
    // through a cage bridge and call it a ring.
    let chorded = false;
    for (let p = 0; p < size && !chorded; p++) {
      for (let q = p + 2; q < size && !chorded; q++) {
        const adjacent = p === 0 && q === size - 1;
        if (adjacent) continue;
        if (adj[path[p]].some(nb => nb.nbr === path[q])) chorded = true;
      }
    }
    if (chorded) continue;

    let input_aromatic = false;
    for (let p = 0; p < size; p++) {
      const o = bond_order(path[p], path[(p + 1) % size]);
      if (o >= 4) input_aromatic = true;
    }
    if (input_aromatic) {
      mark(path);
      continue;
    }
    // Kekulé pattern: every ring atom carries exactly one ring double
    // bond, except that a 5-ring may have exactly one N/O/S with none
    // — the lone-pair donor (pyrrole's N, furan's O, thiophene's S).
    // The 6π count is implied (2 doubles + 1 lone pair in a 5-ring,
    // 3 doubles in a 6-ring); adjacent double bonds and saturated
    // rings (FUVDOP's all-single triazine cage, KAGBOJ's pyrone ring
    // with a bare carbon) fail the pattern.
    let zero_double_atoms = 0;
    let kekule = true;
    for (const a of path) {
      const ring_doubles = adj[a].filter(nb => set.has(nb.nbr) && nb.order === 2).length;
      if (ring_doubles > 1) { kekule = false; break; }
      if (ring_doubles === 0) {
        zero_double_atoms++;
        const el = molecule.atoms[a].element;
        if (el !== 'N' && el !== 'O' && el !== 'S') { kekule = false; break; }
      }
    }
    if (!kekule) continue;
    if (size === 6 && zero_double_atoms !== 0) continue;
    if (size === 5 && zero_double_atoms !== 1) continue;
    mark(path);
  }

  return { atoms, rings_of };
}

/**
 * The lone-pair heteroatoms at the α (ring-neighbor) and β (two bonds
 * away) positions of an aromatic 5-ring atom: S, O, or a pyrrole-type
 * N (3 explicit neighbors, not an N-oxide). Pyridine-type N's (2
 * neighbors) are NOT counted — they carry no lone pair for the ring's
 * π system — so in a thiazole only the S counts, which is why the
 * thiazole C's type differently from pyrrole's.
 */
function five_ring_alpha_beta(
  i: number,
  adj: { nbr: number; order: number }[][],
  molecule: Molecule,
  aromatic: AromaticInfo,
): { alpha: number[]; beta: number[] } {
  const is_hetero = (a: number): boolean => {
    const el = molecule.atoms[a].element;
    if (el === 'S' || el === 'O') return true;
    if (el !== 'N' || adj[a].length !== 3) return false;
    // N-oxide N (bonded to a terminal O) is not a lone-pair donor
    return !adj[a].some(nb => molecule.atoms[nb.nbr].element === 'O' && adj[nb.nbr].length === 1);
  };
  const alpha = new Set<number>();
  const beta = new Set<number>();
  for (const ring of aromatic.rings_of.get(i) ?? []) {
    if (ring.path.length !== 5) continue; // α/β positions live in 5-rings
    const set = new Set(ring.path);
    const idx = ring.path.indexOf(i);
    for (const a of [ring.path[(idx + 1) % 5], ring.path[(idx + 4) % 5]]) {
      if (is_hetero(a)) alpha.add(a);
      // β: a's other ring neighbors in this ring
      for (const b of adj[a]) {
        if (b.nbr !== i && set.has(b.nbr) && is_hetero(b.nbr)) beta.add(b.nbr);
      }
    }
  }
  return { alpha: [...alpha], beta: [...beta] };
}

function shares_ring(a: number, b: number, aromatic: AromaticInfo): boolean {
  const rings_a = aromatic.rings_of.get(a) ?? [];
  const rings_b = aromatic.rings_of.get(b) ?? [];
  return rings_a.some(ra => rings_b.includes(ra));
}

/**
 * Aromatic 5-ring carbon (C5A/C5B/C5): 63 when α to a lone-pair
 * heteroatom (pyrrole's α-C), 64 when β to one (pyrrole's β-C), 78
 * when no heteroatom is in reach — including fused rings whose α and
 * β heteroatoms live in different rings, and rings with only
 * pyridine-type N's (e.g. imidazole's C's beyond the α position).
 */
function type_aromatic_5ring_carbon(
  i: number,
  adj: { nbr: number; order: number }[][],
  molecule: Molecule,
  aromatic: AromaticInfo,
): number {
  const { alpha, beta } = five_ring_alpha_beta(i, adj, molecule, aromatic);
  if (alpha.length === 0 && beta.length === 0) return 78; // C5
  if (alpha.length > 0 && beta.length === 0) return 63; // C5A
  if (alpha.length === 0 && beta.length > 0) return 64; // C5B
  // Both: if an α and a β lie in different rings (fused systems) the
  // α/β concept breaks down — the general 5-ring type applies.
  for (const a of alpha) {
    for (const b of beta) {
      if (!shares_ring(a, b, aromatic)) return 78;
    }
  }
  const s_o = (list: number[]) => list.some(a => {
    const el = molecule.atoms[a].element;
    return el === 'S' || el === 'O';
  });
  if (s_o(alpha)) return 63;
  if (s_o(beta)) return 64;
  return 78;
}

/**
 * Aromatic 5-ring nitrogen: the lone-pair N (pyrrole-type, 3
 * neighbors) with no heteroatom in reach is 39 (pyrrole's N); the
 * pyridine-type N's (2 neighbors) are 65/66 by position; the charged
 * variants (81) and the anionic 76 follow the same positions.
 */
function type_aromatic_5ring_nitrogen(
  i: number,
  adj: { nbr: number; order: number }[][],
  molecule: Molecule,
  aromatic: AromaticInfo,
): number {
  const deg3 = adj[i].length === 3;
  const { alpha, beta } = five_ring_alpha_beta(i, adj, molecule, aromatic);
  if (alpha.length === 0 && beta.length === 0) return deg3 ? 39 : 76; // NPYL / N5M
  if (alpha.length > 0 && beta.length === 0) return deg3 ? 81 : 65; // N5A+ / N5A
  if (alpha.length === 0 && beta.length > 0) return deg3 ? 81 : 66; // N5B+ / N5B
  for (const a of alpha) {
    for (const b of beta) {
      if (!shares_ring(a, b, aromatic)) return 79; // N5
    }
  }
  const s_o = (list: number[]) => list.some(a => {
    const el = molecule.atoms[a].element;
    return el === 'S' || el === 'O';
  });
  if (s_o(alpha)) return 65;
  if (s_o(beta)) return 66;
  return 79;
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
  void molecule;  // accepted for the planned API; unused until implemented
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
