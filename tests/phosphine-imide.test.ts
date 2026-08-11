import { describe, it, expect } from 'vitest';
import type { Molecule } from '../src/types.js';
import { assign_atom_types } from '../src/mmff94/assign-atom-types.js';
import { torsion_terms } from '../src/mmff94/energy/torsion.js';
import { make_class_context } from '../src/mmff94/parameters/parameter-classes.js';

// Methoxyiminophosphine (H–P=N–O–CH₃, PubChem CID 129800975) — the
// 2026-08-10 P=N regression. The planar π system (P=N–O–C) is the
// chemistry; the bug was that P=N typed as 25 (PO4, crd 4 — a
// saturated type) routed the H–P=N–O torsion into the par file's
// zero row 0-0-9-25-0, leaving the P=N rotation free (measured
// 0.11 kcal/mol across 180°) — the refined geometry kept whatever
// twist the embedding start had, and the P's p orbital rendered off
// the N/O π plane. Coordinates are the PubChem 3D conformer (which
// is itself a flat z=0 2D artifact — planar by construction).
const METHOXYIMINOPHOSPHINE: Molecule = {
  atoms: [
    { index: 0, element: 'P', x: 1.8632, y: 0.1259, z: 0.0 },
    { index: 1, element: 'O', x: -0.5303, y: 0.5016, z: 0.0 },
    { index: 2, element: 'N', x: 0.4636, y: -0.4849, z: 0.0 },
    { index: 3, element: 'C', x: -1.7965, y: -0.1426, z: 0.0 },
    { index: 4, element: 'H', x: -2.5751, y: 0.6247, z: 0.0001 },
    { index: 5, element: 'H', x: -1.9099, y: -0.7588, z: 0.8974 },
    { index: 6, element: 'H', x: -1.9100, y: -0.7586, z: -0.8975 },
    { index: 7, element: 'H', x: 2.6057, y: -1.0786, z: -0.0001 },
  ],
  bonds: [
    { atom1: 0, atom2: 2, bond_order: 2 }, // P=N
    { atom1: 0, atom2: 7, bond_order: 1 }, // P–H
    { atom1: 2, atom2: 1, bond_order: 1 }, // N–O
    { atom1: 1, atom2: 3, bond_order: 1 }, // O–C
    { atom1: 3, atom2: 4, bond_order: 1 },
    { atom1: 3, atom2: 5, bond_order: 1 },
    { atom1: 3, atom2: 6, bond_order: 1 },
  ],
};

describe('phosphine imide typing (P=N → 75/62, 2026-08-10)', () => {
  it('types P=N as the ylide pair 75/62, matching OpenBabel', () => {
    const typed = assign_atom_types(METHOXYIMINOPHOSPHINE);
    expect(typed.atom_types[0]).toBe(75); // P — the ylide P (was 25)
    expect(typed.atom_types[2]).toBe(62); // N — the ylidic N⁻ (was 9)
    expect(typed.atom_types[1]).toBe(6);  // O (OR)
    expect(typed.atom_types[3]).toBe(1);  // C (CR)
    expect(typed.atom_types[7]).toBe(71); // H(P) (HS fallback)
  });

  it('gives H–P=N–O the empirical rule-(c) V2 = 3.795 — OpenBabel\'s measured barrier', () => {
    // Rule (c): formal order 2, π = 0.4 (75 mltb 2, 62 mltb 0) →
    // V2 = 6·0.4·√(U_P·U_N) = 6·0.4·√(1.25·2.0) = 3.7947. OpenBabel's
    // HIGH-verbosity log lists exactly 3.795 for this dihedral; with
    // the old 25/9 typing the par zero row 0-0-9-25-0 short-circuited
    // the empirical rules and the rotation was free.
    const typed = assign_atom_types(METHOXYIMINOPHOSPHINE);
    const adj: number[][] = Array.from({ length: typed.atoms.length }, () => []);
    for (const b of typed.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
    const ctx = make_class_context(typed, adj);
    const t = torsion_terms(ctx, typed, 7, 0, 2, 1); // H–P=N–O
    expect(t?.v1).toBe(0);
    expect(t?.v2).toBeCloseTo(3.7947, 3);
    expect(t?.v3).toBe(0);
  });

  it('zeroes the C–O–N=P torsion — the Tinker-arbitrated case (1)', () => {
    // The C–O–N=P torsion (central N(62)–O(6) single bond): both
    // central types carry pilp and no mltb, so rule (g) case (1)
    // applies — no torsion at all. Arbitrated 2026-08-10 with
    // Tinker's ktors (lenovo build): its branch chain zeroes any
    // both-pilp pair (tors1/2/3 = 0) with no mltb requirement, and
    // its prm has no mmfftorsion row for original 169 (class 62).
    // OpenBabel's rule-(g) reading gives this pair V2 = 4.8 (its log
    // lists 4.800) — wrong per Tinker; our old mltb-gated fallback
    // gave rule (h)'s V3 = 0.548 — also wrong. The suite cannot
    // arbitrate: every both-pilp-no-mltb bond in the 761 resolves via
    // a table row or never forms a dihedral.
    const typed = assign_atom_types(METHOXYIMINOPHOSPHINE);
    const adj: number[][] = Array.from({ length: typed.atoms.length }, () => []);
    for (const b of typed.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
    const ctx = make_class_context(typed, adj);
    const t = torsion_terms(ctx, typed, 3, 1, 2, 0); // C–O–N=P
    expect(t).toBeUndefined();
  });
});
