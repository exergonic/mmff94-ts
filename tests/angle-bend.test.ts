import { describe, it, expect } from 'vitest';
import { calc_angle_bend_energy } from '../src/mmff94/energy/angle-bend';
import type { TypedMolecule } from '../src/types';

describe('Angle Bend Energy', () => {
  it('calculates energy for a simple 3-atom system', () => {
    // Using parameter '1-1-1': { k_a: 0.851, theta0: 109.608 }
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'C', x: 1.0, y: 0, z: 0 },
        { index: 1, element: 'C', x: 0, y: 0, z: 0 },
        { index: 2, element: 'C', x: 0, y: 1.0, z: 0 }
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 1, atom2: 2, bond_order: 1 }
      ],
      atom_types: [1, 1, 1]
    };

    const energy = calc_angle_bend_energy(mol);
    // Angle is 90 degrees
    // theta0 = 109.608
    // d_theta = 90 - 109.608 = -19.608
    // E = 0.043844 * 0.851 * (-19.608)^2 = 14.345
    expect(energy).toBeCloseTo(14.345, 3);
  });
});
