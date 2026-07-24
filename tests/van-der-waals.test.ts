import { describe, it, expect } from 'vitest';
import { calc_vdw_energy } from '../src/mmff94/energy/van-der-waals';
import type { TypedMolecule } from '../src/types';

describe('Van der Waals Energy', () => {
  it('returns zero for a single atom', () => {
    const mol: TypedMolecule = {
      atoms: [{ index: 0, element: 'H', x: 0, y: 0, z: 0 }],
      bonds: [],
      atom_types: [5],
    };
    expect(calc_vdw_energy(mol)).toBeCloseTo(0.0, 10);
  });

  it('equals −ε at r = R_ij for two identical type-5 H atoms', () => {
    // For type 5 (H): A_i = 4.2, alpha = 0.25, N = 0.8, G = 1.209
    //
    // Per-atom reduced radius: R_i = A_i * alpha^0.25 = 4.2 * 0.7071 = 2.970
    // For identical atoms: R_ij = R_i (Waldman-Hagler with g=0 gives no correction)
    //
    // At r = R_ij, the buffered 14-7 simplifies to E = -epsilon:
    //
    //   epsilon = 181.16 * G^2 * alpha^2 / (2 * sqrt(alpha/N)) / R_ij^6
    //           = 181.16 * 1.209^2 * 0.25^2
    //             / (2 * sqrt(0.25/0.8)) / 2.970^6
    //           = 16.550 / 1.118 / 686.1
    //           = 0.0216 kcal/mol

    const R_ij = 4.2 * Math.pow(0.25, 0.25); // = 2.970
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'H', x: 0, y: 0, z: 0 },
        { index: 1, element: 'H', x: R_ij, y: 0, z: 0 },
      ],
      bonds: [],
      atom_types: [5, 5],
    };

    const energy = calc_vdw_energy(mol);
    expect(energy).toBeCloseTo(-0.0216, 3);
  });

  it('gives positive energy inside the repulsive wall (r < R_ij)', () => {
    // Two type-5 H atoms at 1.5 A — well inside R_ij = 2.97 A.
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'H', x: 0, y: 0, z: 0 },
        { index: 1, element: 'H', x: 1.5, y: 0, z: 0 },
      ],
      bonds: [],
      atom_types: [5, 5],
    };

    const energy = calc_vdw_energy(mol);
    expect(energy).toBeGreaterThan(0);
  });

  it('approaches zero at large separation (r >> R_ij)', () => {
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'H', x: 0, y: 0, z: 0 },
        { index: 1, element: 'H', x: 50, y: 0, z: 0 },
      ],
      bonds: [],
      atom_types: [5, 5],
    };

    const energy = calc_vdw_energy(mol);
    expect(Math.abs(energy)).toBeLessThan(0.001);
  });

  it('matches OpenBabel reference for ethane (total vdW = 0.20278)', () => {
    // Avogadro-optimized ethane, using all 7 H-H and C-H and C-C pairs.
    // Reference from: obabel ethane.sdf -otxt --ff mmff94 --energy --append "Energy" --log
    const ethane: TypedMolecule = {
      atoms: [
        { index: 0, element: 'C', x: 0.7560, y: 0.0000, z: -0.0000 },
        { index: 1, element: 'C', x: -0.7560, y: -0.0000, z: -0.0000 },
        { index: 2, element: 'H', x: 1.1404, y: -0.5122, z: 0.8871 },
        { index: 3, element: 'H', x: 1.1404, y: -0.5122, z: -0.8871 },
        { index: 4, element: 'H', x: 1.1404, y: 1.0244, z: 0.0000 },
        { index: 5, element: 'H', x: -1.1404, y: 0.5122, z: 0.8871 },
        { index: 6, element: 'H', x: -1.1404, y: 0.5122, z: -0.8871 },
        { index: 7, element: 'H', x: -1.1404, y: -1.0244, z: 0.0000 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 0, atom2: 2, bond_order: 1 },
        { atom1: 0, atom2: 3, bond_order: 1 },
        { atom1: 0, atom2: 4, bond_order: 1 },
        { atom1: 1, atom2: 5, bond_order: 1 },
        { atom1: 1, atom2: 6, bond_order: 1 },
        { atom1: 1, atom2: 7, bond_order: 1 },
      ],
      atom_types: [1, 1, 5, 5, 5, 5, 5, 5],
    };

    const energy = calc_vdw_energy(ethane);
    // Our factor convention (Halgren) vs OpenBabel (MM2) means our
    // bond-stretch and angle-bend values are ~2x OpenBabel's. The vdW
    // formula should match directly since the buffered 14-7 uses the
    // same conversion factors in both implementations.
    expect(energy).toBeCloseTo(0.20278, 2);
  });
});
