import { describe, it, expect } from 'vitest';
import { calc_bond_stretch_energy } from '../src/mmff94/energy/bond-stretch';
import type { TypedMolecule } from '../src/types';

describe('Bond Stretch Energy', () => {
  it('calculates energy for a simple 2-atom system', () => {
    // Using parameter: '0-1-1': { k_b: 4.258, r0: 1.508 }
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
    // r = 1.6, r0 = 1.508, dr = 0.092
    // Harmonic:    143.9325 * 4.258 * 0.092^2 = 5.1875
    // Anharmonic:  1 + cs*dr + 7/12*cs^2*dr^2  with cs=-2
    //            = 1 - 0.184 + 0.01975 = 0.83575
    // E = 5.1875 * 0.83575 = 4.3353
    expect(energy).toBeCloseTo(4.3353, 3);
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
      atom_types: [2, 1] // Deliberately reversed
    };

    const energy = calc_bond_stretch_energy(mol);
    // r = 1.582, r0 = 1.482, dr = 0.1
    // Harmonic:    143.9325 * 4.539 * 0.1^2 = 6.5330
    // Anharmonic:  1 + cs*dr + 7/12*cs^2*dr^2  with cs=-2
    //            = 1 - 0.2 + 0.02333 = 0.82333
    // E = 6.5330 * 0.82333 = 5.3789
    expect(energy).toBeCloseTo(5.3789, 3);
  });
});
