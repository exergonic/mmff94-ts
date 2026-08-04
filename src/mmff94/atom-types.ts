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
  // (O2CM). The thiocarboxylate C(=S)-S⁻ (CS2M) follows the same rule
  // with S in place of O (its =S is the S2CM type 72, not 16).
  // Pre-scanned so the O case works regardless of atom order.
  const carboxylate_carbons = new Set<number>();
  for (let i = 0; i < n; i++) {
    const atom = molecule.atoms[i];
    if (atom.element !== 'C') continue;
    const has_dbl_O = adj[i].some(nb => nb.order === 2 && molecule.atoms[nb.nbr].element === 'O');
    const has_terminal_O = adj[i].some(
      nb => nb.order === 1 && molecule.atoms[nb.nbr].element === 'O' && adj[nb.nbr].length === 1,
    );
    const has_dbl_S = adj[i].some(nb => nb.order === 2 && molecule.atoms[nb.nbr].element === 'S');
    const has_terminal_S = adj[i].some(
      nb => nb.order === 1 && molecule.atoms[nb.nbr].element === 'S' && adj[nb.nbr].length === 1,
    );
    if ((has_dbl_O && has_terminal_O) || (has_dbl_S && has_terminal_S)) carboxylate_carbons.add(i);
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

  // Nitrogen typing pass: the N branch is self-contained (connectivity,
  // elements, aromatic perception, and the order-independent prescans
  // above), so it runs as its own pass — the H-on-N rule needs the N's
  // type regardless of atom order.
  const n_types = new Array<number>(n).fill(-1);
  for (let i = 0; i < n; i++) {
    if (molecule.atoms[i].element === 'N') n_types[i] = type_nitrogen(i, adj, molecule, aromatic, amide_nitrogens);
  }

  for (let i = 0; i < n; i++) {
    const atom = molecule.atoms[i];
    const neighbors = adj[i];
    const n_neighbors = neighbors.length;

    // Check bond orders going out from this atom
    const has_double = neighbors.some(nb => nb.order === 2);
    const has_triple = neighbors.some(nb => nb.order === 3);
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
            ? type_aromatic_5ring_carbon(i, adj, molecule, aromatic, n_types)
            : 37;
          break;
        }

        // Acetylenic: 2 neighbors, triple bond → type 4 (CSP). The
        // cumulene/allene central C (two double bonds — JOFDUD's
        // N=C=C ketenimine) is also sp → 4.
        if (n_neighbors <= 2 && has_triple) {
          atom_types[i] = 4;
          break;
        }
        if (n_neighbors === 2 && double_nbrs.length === 2) {
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
            const dbl_to_N = double_nbrs.some(nb => {
              const target = molecule.atoms[nb.nbr];
              return target.element === 'N';
            });
            const dbl_to_S = double_nbrs.some(nb => {
              const target = molecule.atoms[nb.nbr];
              return target.element === 'S';
            });
            const dbl_to_P = double_nbrs.some(nb => {
              const target = molecule.atoms[nb.nbr];
              return target.element === 'P';
            });
            const n_N_neighbors = neighbors.filter(nb => molecule.atoms[nb.nbr].element === 'N').length;

            if (carboxylate_carbons.has(i)) {
              // Carboxylate carbon C(=O)O⁻ → type 41 (CO₂M)
              atom_types[i] = 41;
            } else if (dbl_to_N && n_N_neighbors >= 2
              && adj[double_nbrs.find(nb => molecule.atoms[nb.nbr].element === 'N')!.nbr].length === 3) {
              // The charged amidinium/guanidinium core: C(=N+)(N)(R)
              // with the =N carrying a third bond AND a second N
              // neighbor (JIFYUS's N+=C(NH2)...). The neutral
              // amidine's =N is 2-coordinate and its C is the plain
              // carbonyl type 3; an iminium C(=N+)(C)(C) with a
              // single N (CIJXOI10) is also 3.
              atom_types[i] = 57;
            } else if (dbl_to_O) {
              // Carbonyl: C=O with 3 neighbors → type 3 (C=O)
              atom_types[i] = 3;
            } else if (dbl_to_S) {
              // Thiocarbonyl C=S (thioamides, thioketones, sulfines)
              // shares the generic carbonyl type 3 (C=O covers C=S in
              // the MMFF94 type list).
              atom_types[i] = 3;
            } else if (dbl_to_N) {
              // The C=N carbon (amidines, imines) is also the general
              // carbonyl type 3 — the type list's C=O entry covers
              // C=N and C=S as well.
              atom_types[i] = 3;
            } else if (dbl_to_P) {
              // The C=P carbon (phosphoalkenes, the P=C(NH2)2 ylides)
              // shares the general carbonyl type 3.
              atom_types[i] = 3;
            } else if (dbl_to_C) {
              // Alkene: C=C with 3 neighbors → type 2 (C=C, vinylic);
              // an olefinic C in a 4-membered ring is CE4R (30).
              atom_types[i] = estimate_ring_size(i, adj, molecule) === 4 ? 30 : 2;
            } else {
              // Other sp² (e.g., C=N) → type 2 as generic sp² default
              atom_types[i] = 2;
            }
            break;
          }

          // 3 neighbors, no double, no aromatic (should be rare for C)
          // — the protonated guanidinium core C(NH₂)₃ is type 57 (CGD+)
          // even without the C=N double bond.
          atom_types[i] = n3_neighbor_count(i, adj, molecule) === 3 ? 57 : 1;
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
          case 'N': {
            // The H type follows the N's own type (computed in the
            // n_types pass): H on a cationic N (34 NR+, 54/55/56
            // iminium/guanidinium, 58 pyridinium, 76/81 imidazolium)
            // → 36 (HNR+); H on an amide (10), delocalized (40) or
            // sulfonamide/cyanamide (43) N → 28 (HNCO/HSP2/HNSO);
            // H on an imine =N (9) → 27 (HN=C); other N–H → 23 (HNR).
            const n_type = n_types[neighbors[0].nbr];
            if (n_type === 34 || n_type === 54 || n_type === 55
              || n_type === 56 || n_type === 58 || n_type === 76 || n_type === 81) {
              atom_types[i] = 36;
            } else if (n_type === 10 || n_type === 40 || n_type === 43 || n_type === 48) {
              atom_types[i] = 28;
            } else if (n_type === 9) {
              atom_types[i] = 27;
            } else {
              atom_types[i] = 23;
            }
            break;
          }
          case 'O': {
            // The O–H hydrogen type follows what ELSE the oxygen is
            // bonded to: water H is 31 (HOH); an H on an oxenium O=+
            // is 52 (HO=+); an acid H on P–OH is 24 (HOP, part I
            // Table III), on S–OH it is 33 (HOS), on a carboxylic
            // acid (the C also carries C=O) it is 24 (HOCO); an
            // enol/phenol H on an aromatic or vinylic C is 29 (HOCC);
            // the generic alcohol H is 21 (HOR).
            const o = neighbors[0].nbr;
            if (water_oxygens.has(o)) {
              atom_types[i] = 31; // H-OH
              break;
            }
            const other = adj[o].find(nb => nb.nbr !== i);
            if (other === undefined) {
              atom_types[i] = 21;
              break;
            }
            if (adj[o].some(nb => nb.order === 2)) {
              atom_types[i] = 52; // H on oxenium O=+
              break;
            }
            const other_atom = molecule.atoms[other.nbr];
            if (other_atom.element === 'P') {
              atom_types[i] = 24; // H-O-P (HOP)
            } else if (other_atom.element === 'S') {
              atom_types[i] = 33; // H-O-S (HOS)
            } else if (other_atom.element === 'C') {
              if (adj[other.nbr].some(b => b.order === 2 && molecule.atoms[b.nbr].element === 'O')) {
                atom_types[i] = 24; // carboxylic acid H (HOCO)
              } else if (
                aromatic.atoms.has(other.nbr)
                || adj[other.nbr].some(b => b.order === 2
                  && ['C', 'N'].includes(molecule.atoms[b.nbr].element))
              ) {
                atom_types[i] = 29; // enol/phenol H (HOCC)
              } else {
                atom_types[i] = 21; // alcohol H (HOR)
              }
            } else {
              atom_types[i] = 21;
            }
            break;
          }
          case 'S':
            // H bonded to sulfur → type 71 (HS)
            atom_types[i] = 71;
            break;
          case 'P':
            // H bonded to phosphorus (PH3, H-P=O) → type 71 (HS —
            // MMFF94 reuses the H-S type for P–H hydrogens)
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
        } else if (n_neighbors === 1 && !has_double) {
          // Terminal single-bonded oxygen. With explicit hydrogens this
          // is never a neutral O–H (that has two neighbors): on carbon it
          // is the oxide O⁻ (type 35, OM — alkoxide or phenoxide); on
          // nitrogen it is the nitro/nitrate O (32), an oxide on an
          // imine/pyridine N (35), or an N-oxide O (32) — the N's
          // environment decides (the O-N bond of -N(=O)-O⁻ etc.).
          const target = molecule.atoms[neighbors[0].nbr];
          if (target.element === 'H') {
            // A bare O–H (no other bonds) is the hydroxide anion
            // → type 35 (OM — oxide). A neutral O–H always has a
            // second bond (water has two H's → 70 above).
            atom_types[i] = 35;
          } else if (target.element === 'C') {
            atom_types[i] = 35; // OM — oxide O on sp3/sp2 C
          } else if (target.element === 'N') {
            let nbrTerminalO = 0;
            let nbrValence = 0;
            for (const b of adj[neighbors[0].nbr]) {
              const bn = molecule.atoms[b.nbr];
              if (bn.element === 'O' && adj[b.nbr].length === 1) nbrTerminalO++;
              nbrValence += b.order;
            }
            if (nbrTerminalO >= 2) atom_types[i] = 32; // nitro/nitrate O
            else if (adj[neighbors[0].nbr].length === 2 || nbrValence === 3) {
              atom_types[i] = 35; // oxide O on imine/pyridine N
            } else {
              atom_types[i] = 32; // N-oxide O (ONX)
            }
          } else if (target.element === 'P' || target.element === 'Cl') {
            // Terminal O on P or Cl is the O2CM type 32 (part I
            // Table III: OP/O2P/O3P/O4P for phosphates, phosphonates
            // and phosphine oxides; O4Cl for perchlorate). A
            // hypochlorite O⁻ is 32 too (OpenBabel probe).
            atom_types[i] = 32;
          } else if (target.element === 'S') {
            // Terminal O on S: 32 when the sulfur carries at least
            // two terminal oxygens (the O2S/O3S/O4S spec entries —
            // sulfonate, sulfate, sulfinate anion, sulfite dianion);
            // a lone terminal O⁻ on a sulfide-type S has no suite
            // case and stays the neutral 6.
            atom_types[i] = count_terminal_oxygens(neighbors[0].nbr, adj, molecule) >= 2 ? 32 : 6;
          } else {
            atom_types[i] = 6;
          }
        } else if (n_neighbors === 2 && has_double) {
          // Oxenium O⁺: one double bond and one single (pyrylium)
          // → type 51 (O=+)
          atom_types[i] = 51;
        } else if (has_double) {
          // Double-bonded oxygen. O=C and O=N stay the generic 7, but
          // terminal oxygens on P and on the sulfone family are the
          // O2CM type 32 (spec: OP and O2S/O3S/O4S — Halgren types
          // oxyanion O's as O2CM, not as carbonyl O). Sulfoxide-type
          // S=O (one terminal O, no S=N) keeps 7.
          const target = molecule.atoms[neighbors[0].nbr];
          if (target.element === 'P') {
            atom_types[i] = 32; // P=O — phosphine oxide, phosphate
          } else if (target.element === 'N') {
            // =O on N: 32 on the nitro/nitrate group (≥ 2 terminal
            // O's on the N — FUCTIG01's nitrate ion) and on a
            // 3-coordinate N-oxide N (pyridine N-oxide drawn with a
            // formal double); 7 on the 2-coordinate nitroso N=O
            // (nitroso compounds, nitrosamines).
            const n_atom = neighbors[0].nbr;
            const n_has_two_O = count_terminal_oxygens(n_atom, adj, molecule) >= 2;
            atom_types[i] = n_has_two_O || adj[n_atom].length === 3 ? 32 : 7;
          } else if (target.element === 'S') {
            // S=O: 32 on the sulfone family (≥ 2 terminal O's on the
            // S — sulfone, sulfonate, sulfate, sulfite/sulfinate
            // anions) and on sulfonyl imines (S=N double); 7 on the
            // sulfoxide family (one terminal O — sulfoxide, sulfinic
            // acid/ester, sulfine C=S=O, thiosulfinate).
            const s = neighbors[0].nbr;
            const s_has_dbl_N = adj[s].some(b => b.order === 2 && molecule.atoms[b.nbr].element === 'N');
            atom_types[i] = count_terminal_oxygens(s, adj, molecule) >= 2 || s_has_dbl_N ? 32 : 7;
          } else {
            atom_types[i] = 7; // O=C or O=N
          }
        } else {
          // Ether or alcohol oxygen → type 6
          atom_types[i] = 6;
        }
        break;
      }

      // ── Nitrogen ────────────────────────────────────────────────────
      // Assigned in the n_types pass above (its H rule needs the N
      // type regardless of atom order).
      case 'N':
        atom_types[i] = n_types[i];
        break;

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
        if (has_double && n_neighbors === 2) {
          // Two-coordinate S with a double bond: the sulfine C=S=O
          // (74, =S=O — "sulfinyl sulfur").
          atom_types[i] = 74;
        } else if (has_double && n_neighbors >= 3) {
          // Check for S=O / SO₂
          const dbl_to_O = double_nbrs.some(nb => {
            const target = molecule.atoms[nb.nbr];
            return target.element === 'O';
          });
          if (dbl_to_O && n_neighbors >= 4) {
            atom_types[i] = 18; // SO2 (sulfone)
          } else if (dbl_to_O && n_neighbors === 3) {
            // 3-coordinate S with an S=O double bond: the sulfoxide
            // (17) unless it also carries two terminal oxygens — the
            // anionic sulfinate R-S(=O)-O⁻ (73, SO2M); a terminal O
            // plus a terminal S is the thiosulfinate (also 73). The
            // R-S(=O)2⁻ form with a doubly bonded carbon (SURDOX02's
            // sulfonyl carbene) is the sulfone 18, not 73.
            let terminalO = 0;
            let terminalS = 0;
            for (const nb of neighbors) {
              const t = molecule.atoms[nb.nbr];
              if (t.element === 'O' && adj[nb.nbr].length === 1) terminalO++;
              if (t.element === 'S' && adj[nb.nbr].length === 1) terminalS++;
            }
            if (terminalO === 2) {
              const nonO = neighbors.find(nb => molecule.atoms[nb.nbr].element !== 'O');
              atom_types[i] = nonO !== undefined && nonO.order === 1 ? 73 : 18;
            } else if (terminalO >= 1 && terminalS >= 1) {
              atom_types[i] = 73; // SO2M — sulfinate S
            } else {
              atom_types[i] = 17; // S=O (sulfoxide)
            }
          } else if (dbl_to_O) {
            atom_types[i] = 17; // S=O (sulfoxide)
          } else if (n_neighbors === 3 && double_nbrs.some(nb => molecule.atoms[nb.nbr].element === 'N')) {
            // 3-coordinate S with an S=N double and no oxygen: the
            // sulfinimide S (FIZGEA's ring S) shares the sulfoxide
            // type 17.
            atom_types[i] = 17;
          } else {
            atom_types[i] = 16; // S=C (thiocarbonyl)
          }
        } else if (n_neighbors === 1) {
          // Terminal S: with a single bond it is the anionic 72 (SM).
          // A C=S double bond is the thiocarbonyl 16 (S=C — thioamides,
          // thioketones), while P=S and S=S stay 72 (thiophosphate and
          // disulfide anions use the same anionic type); the =S of a
          // thiocarboxylate C(=S)-S⁻ is also 72 (S2CM).
          const bond = neighbors[0];
          const is_thiocarbonyl = bond.order === 2 && molecule.atoms[bond.nbr].element === 'C'
            && !carboxylate_carbons.has(bond.nbr);
          atom_types[i] = is_thiocarbonyl ? 16 : 72;
        } else {
          // Thiol or sulfide → type 15 (S)
          atom_types[i] = 15;
        }
        break;
      }

      // ── Halogens ────────────────────────────────────────────────────
      // A bare halogen atom is the halide anion (89 F⁻, 90 Cl⁻, 91 Br⁻);
      // 4-coordinate Cl with four oxygens is perchlorate (77).
      case 'F':
        atom_types[i] = n_neighbors === 0 ? 89 : 11;
        break;
      case 'Cl':
        atom_types[i] = n_neighbors === 0 ? 90
          : n_neighbors === 4 && neighbors.every(nb => molecule.atoms[nb.nbr].element === 'O') ? 77
          : 12;
        break;
      case 'Br':
        atom_types[i] = n_neighbors === 0 ? 91 : 13;
        break;
      case 'I':  atom_types[i] = 14; break;

      // ── Silicon, phosphorus ─────────────────────────────────────────
      case 'Si': atom_types[i] = 19; break;
      case 'P': {
        // 4-coordinate P, or P with a double bond to O or S, is the
        // phosphate family (25: the spec's PO4/PO3/PO2/PO/PTET all
        // share it — phosphates, phosphonates, phosphine oxides);
        // tricoordinate P(III) is the phosphine 26; a C=P double
        // bond is the ylide 75 (-P=C).
        if (has_double && double_nbrs.some(nb => molecule.atoms[nb.nbr].element === 'C')) {
          atom_types[i] = 75; // -P=C
        } else if (n_neighbors >= 4 || has_double) {
          atom_types[i] = 25; // PO4 family
        } else {
          atom_types[i] = 26; // P (phosphine)
        }
        break;
      }

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

  const is_donor = (a: number): boolean => {
    const el = molecule.atoms[a].element;
    const has_terminal_O = adj[a].some(nb => {
      return molecule.atoms[nb.nbr].element === 'O' && adj[nb.nbr].length === 1;
    });
    if (el === 'N') return !has_terminal_O;
    if (el === 'O') return adj[a].length === 2;
    if (el === 'S') return adj[a].length === 2 && !has_terminal_O;
    return false;
  };

  // π evaluation: every ring atom must be conjugated — exactly one
  // ring double bond, OR a lone-pair donor, OR an exocyclic double to
  // a ring C/N (the fusion atom of a fused aromatic: the mmd writes
  // its π bond into the adjacent ring — benzimidazole's C3=C4/C7=C8,
  // benzofuran's C13=C14). The donor is an N (not an N-oxide), a
  // 2-coordinate O (furan), or a 2-coordinate S without terminal
  // oxygens — sulfone/sulfoxide sulfurs cannot donate, which keeps
  // the sulfolenes (FIZGOK/FIZGEA) non-aromatic. An exocyclic double
  // to a chain atom is a separate π system (COCXUN's allenyl-nitrile
  // substituents) and does not conjugate the ring. The π sum must be
  // 6 (donors count 2, conjugated atoms 1); anything else — an sp3
  // atom — breaks the ring.
  const eval_pi = (path: number[]): Map<number, number> | null => {
    const set = new Set(path);
    let pi = 0;
    const exo = new Map<number, number>();
    for (const a of path) {
      const ring_doubles = adj[a].filter(nb => set.has(nb.nbr) && nb.order === 2).length;
      if (ring_doubles > 1) return null;
      if (ring_doubles === 1) { pi += 1; continue; }
      if (is_donor(a)) { pi += 2; continue; }
      const e = adj[a].find(nb => !set.has(nb.nbr) && nb.order === 2);
      if (e && (molecule.atoms[e.nbr].element === 'C' || molecule.atoms[e.nbr].element === 'N')
        && is_ring[e.nbr]) {
        pi += 1;
        exo.set(a, e.nbr);
        continue;
      }
      return null;
    }
    return pi === 6 ? exo : null;
  };

  // The fused-aromatic exo allowance is mutually recursive (a
  // benzimidazole's two rings each carry the other's π bond), so the
  // candidate rings converge as a fixpoint: a ring survives only when
  // every exo partner lies in a ring that also survives. COGDEH's
  // first naphthalene ring is excluded this way — its exo partners
  // (the =N's) belong to the as-triazine ring, which is 7π and never
  // a candidate.
  const ring_candidates: { path: number[]; exo: Map<number, number> }[] = [];
  const ring_kept: boolean[] = [];
  for (const path of cycles) {
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
      ring_candidates.push({ path, exo: new Map() });
      ring_kept.push(true);
      continue;
    }
    const exo = eval_pi(path);
    if (exo) {
      ring_candidates.push({ path, exo });
      ring_kept.push(true);
    }
  }

  let changed = true;
  while (changed) {
    changed = false;
    for (let k = 0; k < ring_candidates.length; k++) {
      if (!ring_kept[k]) continue;
      for (const partner of ring_candidates[k].exo.values()) {
        const partner_kept = ring_candidates.some(
          (r2, k2) => ring_kept[k2] && r2.path.includes(partner),
        );
        if (!partner_kept) {
          ring_kept[k] = false;
          changed = true;
          break;
        }
      }
    }
  }

  for (let k = 0; k < ring_candidates.length; k++) {
    if (ring_kept[k]) mark(ring_candidates[k].path);
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
  n_types: number[],
): number {
  const { alpha, beta } = five_ring_alpha_beta(i, adj, molecule, aromatic);
  // The carbon between two 3-coordinate N's where at least one is an
  // imidazolium N (81) is the imidazolium carbon (CIM+). The n_types
  // gate replaces the old same-ring α check: CUDREY's 2-aminoimidazolium
  // C (ring N + exocyclic NH2) is 80, while FARMAM's bridgehead
  // (two pyrrole-type 39 N's in different rings) stays 63/64, and a
  // 2-aminopyrrole α-C (39 + exocyclic NH2) stays 63.
  const n3_nbrs = adj[i].filter(nb => {
    return molecule.atoms[nb.nbr].element === 'N' && adj[nb.nbr].length === 3;
  });
  if (n3_nbrs.length >= 2 && n3_nbrs.some(nb => n_types[nb.nbr] === 81)) {
    return 80; // CIM+
  }
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
 * Count the 3-coordinate N neighbors of atom i — the "N3 count" that
 * anchors the guanidinium rules (a carbon with three N3 neighbors is
 * the guanidinium core; an N3 attached to it is NGD+).
 */
function n3_neighbor_count(
  i: number,
  adj: { nbr: number; order: number }[][],
  molecule: Molecule,
): number {
  let count = 0;
  for (const nb of adj[i]) {
    if (molecule.atoms[nb.nbr].element === 'N' && adj[nb.nbr].length === 3) count++;
  }
  return count;
}

/**
 * The nitrogen typing decision tree, run as its own pass (the H-on-N
 * rule needs the N's type regardless of atom order). The branch is
 * self-contained: it reads only connectivity, elements, the aromatic
 * perception, and the order-independent amide prescan.
 */
function type_nitrogen(
  i: number,
  adj: { nbr: number; order: number }[][],
  molecule: Molecule,
  aromatic: AromaticInfo,
  amide_nitrogens: Set<number>,
): number {
  const neighbors = adj[i];
  const n_neighbors = neighbors.length;
  const has_double = neighbors.some(nb => nb.order === 2);
  const has_aromatic = neighbors.some(nb => nb.order >= 4);
  const double_nbrs = neighbors.filter(nb => nb.order === 2);

  // Aromatic ring nitrogen: 38 in a 6-ring (pyridine), 39/65/66
  // in a 5-ring by position (pyrrole's N is 39; imidazole- and
  // thiazole-type N's are 65/66). The charged variants (58/76/81)
  // follow the same positions (part II).
  if (aromatic.atoms.has(i)) {
    const in_5ring = (aromatic.rings_of.get(i) ?? []).some(r => r.path.length === 5);
    return in_5ring
      ? type_aromatic_5ring_nitrogen(i, adj, molecule, aromatic)
      : adj[i].length === 3 ? 58 : 38;
  }

  // Amide N (type 10, NC=O): single-bonded to a carbonyl carbon —
  // trigonal by resonance, not by an N=double bond (formamide).
  // See the amide_nitrogens pre-scan above.
  if (amide_nitrogens.has(i)) return 10;

  // Terminal nitrile N (C≡N) → 42 (NSP). The sp carbon itself is 4.
  if (n_neighbors === 1 && neighbors[0].order === 3) return 42;

  if (n_neighbors === 4) {
    // Quaternary N: 4 single bonds → type 34 (NR+). A terminal
    // oxygen neighbor makes it the sp3 N-oxide → type 68 (N3OX).
    const has_terminal_O = neighbors.some(nb => {
      const t = molecule.atoms[nb.nbr];
      return t.element === 'O' && adj[nb.nbr].length === 1;
    });
    return has_terminal_O ? 68 : 34;
  }

  if (n_neighbors === 3 && !has_double && !has_aromatic) {
    // Amine N with 3 single bonds. The delocalized NC=C family (40):
    // an N whose lone pair conjugates with an sp2 center — the
    // amidine/guanidine N-C=N, enamine/aniline N-C=C, and the
    // ylide-adjacent P=C(NH2)2. The guanidinium (56) and amidinium
    // (55) cations are the charged members; the neutral amidine's N
    // is plain 40. Sulfonamide (43, NSO2/NSO3 — also the
    // phosphonamide N-PO2 and the cyanamide H2N-C≡N) takes
    // precedence over the neutral delocalized type, but the charged
    // amidinium/guanidinium wins over the sulfonyl (JIFYUS's
    // amidinium N also carries the sulfonate counterion → 55, not 43).
    let charged = 0;
    let sulfonyl = false;
    let cyanamide = false;
    let azo_amide = false;
    let delocalized = false;
    for (const nb of neighbors) {
      const nbr = molecule.atoms[nb.nbr];
      if (nbr.element === 'N') {
        // The triazene N–N=N: the lone pair conjugates with the azo
        // group — the reference types it like an amide (10).
        if (adj[nb.nbr].some(b => b.order === 2 && (molecule.atoms[b.nbr].element === 'N' || molecule.atoms[b.nbr].element === 'O'))) {
          azo_amide = true;
        }
      }
      if ((nbr.element === 'S' || nbr.element === 'P')
        && count_terminal_oxygens(nb.nbr, adj, molecule) >= 2) {
        sulfonyl = true; // NSO2 — sulfonamide/sulfamate/phosphonamide N
      }
      if (nbr.element === 'C' && adj[nb.nbr].some(b => b.order === 3 && molecule.atoms[b.nbr].element === 'N')) {
        cyanamide = true; // H2N-C≡N
      }
      if (nbr.element !== 'C') continue;
      // The amidine core: the double-bonded N's coordination decides
      // the subtype — the charged amidinium carries a third bond on
      // its =N; the neutral amidine's =N is 2-coordinate.
      let n3 = 0;
      let dbl_N = -1;
      let dbl_to_C = false;
      let dbl_to_P = false;
      for (const b of adj[nb.nbr]) {
        const bn = molecule.atoms[b.nbr];
        if (bn.element === 'N' && adj[b.nbr].length === 3) n3++;
        if (b.order === 2 && bn.element === 'N') dbl_N = b.nbr;
        if (b.order === 2 && bn.element === 'C') dbl_to_C = true;
        if (b.order === 2 && bn.element === 'P') dbl_to_P = true;
      }
      if (n3 === 3) { charged = 56; break; } // NGD+ — guanidinium N
      if (dbl_N >= 0) {
        // NCN+ — amidinium N: the =N carries a third bond (JIFYUS's
        // acyclic amidinium, CUDREY's 5-ring imidazolium). A 6-ring
        // aromatic =N (the pyridinium-type ring N of a protonated
        // adenine, COJFIQ) is not an amidinium — its amino N stays
        // the delocalized 40.
        const dblN_in_6ring = (aromatic.rings_of.get(dbl_N) ?? []).some(r => r.path.length === 6);
        if (adj[dbl_N].length === 3 && !dblN_in_6ring) {
          charged = 55;
        } else {
          delocalized = true; // neutral amidine NC=N
        }
      }
      if (dbl_to_C) delocalized = true; // enamine/aniline NC=C
      if (dbl_to_P) delocalized = true; // the P=C(NH2)2 ylide N
    }
    if (charged) return charged;
    if (sulfonyl || cyanamide) return 43; // NSO2 / cyanamide
    // The triazene N–N=N types like an amide (10) — but only when it
    // does not also conjugate with an sp2 carbon: SEMXOX's
    // tetrazolyl-alkene N (azo + C=C) is the delocalized 40.
    if (azo_amide && !delocalized) return 10;
    if (delocalized) return 40; // NC=C — delocalized lone pair
    return 8; // NR — plain amine
  }

  if (n_neighbors === 2 && has_double) {
    // The nitroso N=O (46, NNO) — nitrosomethane, the dinitroso
    // ring N's.
    if (double_nbrs[0].nbr !== undefined && molecule.atoms[double_nbrs[0].nbr].element === 'O') {
      return 46;
    }
    // The N=S family: an N double-bonded to a sulfonyl/sulfoximine S
    // is the sulfonyl imine 48 (FIZGOK's =N–S(=O)2); an N
    // double-bonded to a sulfoxide-type S while single-bonded to a
    // true sulfonyl is the sulfonamide 43 (FIZGEA's =N–S–SO2). An
    // N=C imine stays 9 even with an SO2 substituent (CAGREH10's
    // ring N — the imine wins over the sulfonyl).
    const dbl_nbr = double_nbrs[0].nbr;
    const dbl_el = molecule.atoms[dbl_nbr].element;
    if (dbl_el === 'C' || dbl_el === 'N') return 9; // imine / azo
    // The double goes to S.
    const to = count_terminal_oxygens(dbl_nbr, adj, molecule);
    const s_has_dbl_N = adj[dbl_nbr].some(b => b.order === 2 && molecule.atoms[b.nbr].element === 'N');
    if (to >= 2 || (to >= 1 && s_has_dbl_N)) return 48; // NS2O — N=S(=O)2 / sulfoximine
    for (const nb of neighbors) {
      if (nb.nbr === dbl_nbr) continue;
      const t = molecule.atoms[nb.nbr];
      if (t.element === 'S' && count_terminal_oxygens(nb.nbr, adj, molecule) >= 2) {
        return 43; // NSO2 — N bonded to a sulfonyl
      }
    }
    return 9;
  }

  if (n_neighbors === 3 && has_double) {
    // The sulfonyl imine with an extra H/substituent (the
    // sulfoximine's =N–S(=O)) is also 48.
    for (const nb of neighbors) {
      const t = molecule.atoms[nb.nbr];
      if (nb.order === 2 && t.element === 'S' && count_terminal_oxygens(nb.nbr, adj, molecule) >= 1) {
        const s_has_dbl_N = adj[nb.nbr].some(b => b.order === 2 && molecule.atoms[b.nbr].element === 'N');
        if (s_has_dbl_N) return 48;
      }
    }
    // N with three neighbors AND its own double bond. Count the
    // terminal oxygens first: one → the sp2 N-oxide (67, N2OX),
    // two or more → nitro/nitrate N (45, NO2/NO3). Otherwise the
    // double goes to C or N and the iminium family applies: the
    // N3-neighbors of the double-bonded carbon fix the subtype —
    // 1 (or 0) → iminium (54, N+=C), 2 → the NCN+ pair (55),
    // 3 → guanidinium N (56, NGD+); N+=N is also 54.
    let terminalO = 0;
    for (const nb of neighbors) {
      const t = molecule.atoms[nb.nbr];
      if (t.element === 'O' && adj[nb.nbr].length === 1) terminalO++;
    }
    if (terminalO === 1) {
      return 67; // N2OX — sp2 N-oxide
    }
    if (terminalO >= 2) {
      return 45; // NO2/NO3 — nitro/nitrate N
    }
    const dbl_nbr = double_nbrs[0].nbr;
    if (molecule.atoms[dbl_nbr].element === 'N') {
      return 54; // N+=N — diazenium
    }
    let n3 = 0;
    for (const b of adj[dbl_nbr]) {
      const bn = molecule.atoms[b.nbr];
      if (bn.element === 'N' && adj[b.nbr].length === 3) n3++;
    }
    return n3 === 3 ? 56 : n3 === 2 ? 55 : 54;
  }

  return 8; // fallback: amine N
}

/**
 * Count the terminal oxygens on atom i — O's whose only bond is to i.
 * This is the count that decides the O2CM assignment for S=O and S-O⁻:
 * a sulfone/sulfonate/sulfate/sulfinate-anion sulfur carries ≥ 2
 * terminal oxygens (32), a sulfoxide/sulfinic-acid/sulfite-ester
 * sulfur carries exactly 1 (its S=O stays the carbonyl-like 7).
 */
function count_terminal_oxygens(
  i: number,
  adj: { nbr: number; order: number }[][],
  molecule: Molecule,
): number {
  let count = 0;
  for (const nb of adj[i]) {
    if (molecule.atoms[nb.nbr].element === 'O' && adj[nb.nbr].length === 1) count++;
  }
  return count;
}

/**
 * Estimate the size of the smallest ring containing atom i.
 *
 * BFS from i with a depth map. A cycle closes when an edge joins two
 * BFS-tree nodes (length = d(x) + d(y) + 1) — including the case where
 * two children of the start node are bonded to each other, which the
 * naive visited-set BFS misses (a cyclopropane's triangle closes
 * exactly that way: both ring neighbors are depth-1 children of the
 * start). The tree-child edge (the reverse of a parent edge) is
 * skipped, so no false cycle is counted. Bounded to depth 6; returns
 * 0 if no ring is found within the bound.
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
  const depth = new Map<number, number>();
  const parent = new Map<number, number>();
  const queue: number[] = [];
  depth.set(start, 0);
  parent.set(start, -1);
  queue.push(start);

  let best = 7; // bound: only 3- and 4-membered rings matter
  let head = 0;
  while (head < queue.length) {
    const node = queue[head++];
    const d = depth.get(node)!;
    if (d > 5) continue; // don't search beyond 6-membered rings

    for (const { nbr } of adj[node]) {
      if (nbr === parent.get(node)) continue; // the tree edge back up
      const dn = depth.get(nbr);
      if (dn === undefined) {
        depth.set(nbr, d + 1);
        parent.set(nbr, node);
        queue.push(nbr);
      } else if (parent.get(nbr) !== node) {
        // A second path to an already-visited node closes a cycle
        // through start: start→...→node + node–nbr + nbr→...→start.
        best = Math.min(best, d + dn + 1);
      }
    }
  }

  return best <= 6 ? best : 0;
}
