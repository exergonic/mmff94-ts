import { describe, it, expect } from 'vitest';
import { calc_bond_stretch_energy } from '../src/mmff94/energy/bond-stretch';
import type { TypedMolecule } from '../src/types';

describe('Bond Stretch Energy', () => {
  it('calculates energy for a simple 2-atom system', () => {
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'C', x: 0, y: 0, z: 0 },
        { index: 1, element: 'C', x: 1.6, y: 0, z: 0 }
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 }
      ],
      atom_types: [1, 1]
    };

    const energy = calc_bond_stretch_energy(mol);
    // r = 1.6, r0 = 1.508, dr = 0.092, k_b = 4.258, cs = -2
    // ½k_b        = 2.129
    // Harmonic:   143.9325 * 2.129 * 0.092^2 = 2.5938
    // Anharmonic: 1 + (-2)*0.092 + 7/12*4*0.092^2 = 0.83575
    // E = 2.5938 * 0.83575 = 2.1677
    expect(energy).toBeCloseTo(2.1677, 3);
  });

  it('handles reversed type keys (min-max ordering)', () => {
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'C', x: 0, y: 0, z: 0 },
        { index: 1, element: 'C', x: 1.582, y: 0, z: 0 }
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 }
      ],
      atom_types: [2, 1]
    };

    const energy = calc_bond_stretch_energy(mol);
    // r = 1.582, r0 = 1.482, dr = 0.1, k_b = 4.539, cs = -2
    // ½k_b        = 2.2695
    // Harmonic:   143.9325 * 2.2695 * 0.1^2 = 3.2665
    // Anharmonic: 1 + (-2)*0.1 + 7/12*4*0.01 = 0.82333
    // E = 3.2665 * 0.82333 = 2.6895
    expect(energy).toBeCloseTo(2.6895, 3);
  });
});
