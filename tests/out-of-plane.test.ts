import { describe, it, expect } from 'vitest';
import { calc_oop_energy } from '../src/mmff94/energy/out-of-plane';
import { wilson_oop_angle } from '../src/utils/vector';
import type { TypedMolecule } from '../src/types';

// Symmetric trigonal pyramid: central atom j at the origin, three
// substituents at azimuths 0°/120°/240° and polar angle α from the +z
// axis. By symmetry all three Wilson angles are equal, with
//
//   sin χ = (3/2) · sin α · cos α / √(cos²α + sin²α/4)
//
// (unit-normal projection of one bond onto the plane of the other two).
// All three angles share one force constant, so
//
//   E = 3 · 0.043844 · (k_oop/2) · χ²
function pyramid_at_alpha(alpha_deg: number): { chi_deg: number } {
  const alpha = (alpha_deg * Math.PI) / 180;
  const s = Math.sin(alpha);
  const c = Math.cos(alpha);
  const sin_chi = (1.5 * s * c) / Math.sqrt(c * c + (s * s) / 4);
  return { chi_deg: (Math.asin(sin_chi) * 180) / Math.PI };
}

// Trigonal pyramid with polar angle α, typed via an explicit type list.
function pyramid_molecule(alpha_deg: number, atom_types: number[]): TypedMolecule {
  const alpha = (alpha_deg * Math.PI) / 180;
  const s = Math.sin(alpha);
  const c = Math.cos(alpha);
  return {
    atoms: [
      { index: 0, element: 'C', x: s, y: 0, z: c },                 // i
      { index: 1, element: 'C', x: 0, y: 0, z: 0 },                 // j (central)
      { index: 2, element: 'C', x: -s / 2, y: (Math.sqrt(3) * s) / 2, z: c }, // k
      { index: 3, element: 'C', x: -s / 2, y: (-Math.sqrt(3) * s) / 2, z: c }, // l
    ],
    bonds: [
      { atom1: 0, atom2: 1, bond_order: 1 },
      { atom1: 1, atom2: 2, bond_order: 1 },
      { atom1: 1, atom2: 3, bond_order: 1 },
    ],
    atom_types,
  };
}

const OOP_UNIT = 0.043844;

describe('wilson_oop_angle', () => {
  it('returns 0 for a planar center', () => {
    const j: [number, number, number] = [0, 0, 0];
    const i: [number, number, number] = [1, 0, 0];
    const k: [number, number, number] = [0, 1, 0];
    const l: [number, number, number] = [-0.6, -0.8, 0]; // in the xy plane
    expect(wilson_oop_angle(i, j, k, l)).toBeCloseTo(0.0, 10);
  });

  it('measures the angle between the bond and the plane (60° case)', () => {
    // Plane (i, j, k) is the xy plane; l sits 0.866 above it with |jl| = 1,
    // so the bond j→l makes exactly 60° with the plane.
    const j: [number, number, number] = [0, 0, 0];
    const i: [number, number, number] = [1, 0, 0];
    const k: [number, number, number] = [0, 1, 0];
    const l: [number, number, number] = [0.3, 0.4, Math.sqrt(0.75)];
    expect(wilson_oop_angle(i, j, k, l)).toBeCloseTo(60.0, 9);
  });

  it('sign records which side of the plane the atom sits on', () => {
    const j: [number, number, number] = [0, 0, 0];
    const i: [number, number, number] = [1, 0, 0];
    const k: [number, number, number] = [0, 1, 0];
    const l: [number, number, number] = [0.3, 0.4, -Math.sqrt(0.75)];
    expect(wilson_oop_angle(i, j, k, l)).toBeCloseTo(-60.0, 9);
  });

  it('returns 0 when the two in-plane bonds are collinear', () => {
    const j: [number, number, number] = [0, 0, 0];
    const i: [number, number, number] = [1, 0, 0];
    const k: [number, number, number] = [-1, 0, 0];
    const l: [number, number, number] = [0, 0, 1];
    expect(wilson_oop_angle(i, j, k, l)).toBeCloseTo(0.0, 10);
  });
});

describe('Out-of-Plane Energy', () => {
  it('returns zero when no atom has three neighbors', () => {
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'C', x: 0, y: 0, z: 0 },
        { index: 1, element: 'C', x: 1.5, y: 0, z: 0 },
      ],
      bonds: [{ atom1: 0, atom2: 1, bond_order: 1 }],
      atom_types: [1, 1],
    };
    expect(calc_oop_energy(mol)).toBeCloseTo(0.0, 10);
  });

  it('returns zero for a perfectly planar center (formaldehyde-like)', () => {
    // C=O with two H, all in the z=0 plane. The planar geometry makes
    // every Wilson angle zero regardless of k_oop.
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'O', x: 1.2, y: 0, z: 0 },     // type 7 (O=C)
        { index: 1, element: 'C', x: 0, y: 0, z: 0 },       // type 3 (carbonyl C)
        { index: 2, element: 'H', x: -0.6, y: 1.0, z: 0 },  // type 5
        { index: 3, element: 'H', x: -0.6, y: -1.0, z: 0 }, // type 5
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 2 },
        { atom1: 1, atom2: 2, bond_order: 1 },
        { atom1: 1, atom2: 3, bond_order: 1 },
      ],
      atom_types: [7, 3, 5, 5],
    };
    expect(calc_oop_energy(mol)).toBeCloseTo(0.0, 10);
  });

  it('applies the specific force constant for a pyramidal carbonyl carbon', () => {
    // Carbonyl C (type 3) with substituents {5, 5, 7} matches the
    // specific parameter '5-3-5-7' (k_oop = 0.103).
    const { chi_deg } = pyramid_at_alpha(80);
    const mol = pyramid_molecule(80, [5, 3, 5, 7]);

    const expected = 3 * OOP_UNIT * (0.103 / 2) * chi_deg * chi_deg;
    expect(calc_oop_energy(mol)).toBeCloseTo(expected, 6);
  });

  it('falls back to the wildcard force constant when no specific entry exists', () => {
    // Vinylic C (type 2) with substituents {1, 1, 5} has no specific
    // parameter, so the wildcard default '0-2-0-0' (k_oop = 0.02) applies.
    const { chi_deg } = pyramid_at_alpha(80);
    const mol = pyramid_molecule(80, [1, 2, 1, 5]);

    const expected = 3 * OOP_UNIT * (0.02 / 2) * chi_deg * chi_deg;
    expect(calc_oop_energy(mol)).toBeCloseTo(expected, 6);
  });

  it('gives zero for a pyramidal amine nitrogen (k_oop = 0)', () => {
    // Amine N (type 8) has an explicit zero force constant: the angle-bend
    // reference angles alone decide how pyramidal the center is. This is
    // the tri-coordinate-but-not-planar case: the term IS evaluated, it
    // just contributes nothing.
    const mol = pyramid_molecule(70, [1, 8, 1, 1]);
    expect(calc_oop_energy(mol)).toBeCloseTo(0.0, 10);
  });

  it('gives negative energy for a pyramidal amide nitrogen', () => {
    // Amide N (type 10) with substituents {1, 3, 28} matches the specific
    // parameter '1-10-3-28' with a NEGATIVE k_oop (-0.02). MMFF94 (unlike
    // MMFF94s) deliberately favors pyramidalized amide nitrogen to match
    // the nonplanar MP2 geometries used in parameterization.
    const { chi_deg } = pyramid_at_alpha(70);
    const mol = pyramid_molecule(70, [1, 10, 3, 28]);

    const expected = 3 * OOP_UNIT * (-0.02 / 2) * chi_deg * chi_deg;
    const energy = calc_oop_energy(mol);
    expect(energy).toBeCloseTo(expected, 6);
    expect(energy).toBeLessThan(0);
  });

  it('skips tri-coordinate centers that have no out-of-plane parameters', () => {
    // Alkane C (type 1) has no oop parameters at all, so a hypothetical
    // tri-coordinate sp3 carbon contributes nothing.
    const mol = pyramid_molecule(80, [1, 1, 1, 5]);
    expect(calc_oop_energy(mol)).toBeCloseTo(0.0, 10);
  });

  it('is invariant to substituent ordering in the adjacency list', () => {
    // Same geometry and types, but the bonds are stored in a different
    // order, so the center's neighbor list arrives permuted. The lookup
    // matches the sorted multiset, so the energy must not change.
    const alpha = (80 * Math.PI) / 180;
    const s = Math.sin(alpha);
    const c = Math.cos(alpha);
    const base = pyramid_molecule(80, [5, 3, 5, 7]);
    const permuted: TypedMolecule = {
      ...base,
      bonds: [
        { atom1: 1, atom2: 3, bond_order: 1 }, // l first
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 1, atom2: 2, bond_order: 1 },
      ],
    };
    expect(calc_oop_energy(permuted)).toBeCloseTo(calc_oop_energy(base), 10);
    // Sanity: the permuted geometry is still a pyramid of the same size
    expect(s).toBeGreaterThan(0);
    expect(c).toBeGreaterThan(0);
  });
});
