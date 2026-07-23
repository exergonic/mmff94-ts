import { describe, it, expect } from 'vitest';
import { calc_bond_stretch_energy } from '../src/mmff94/energy/bond-stretch';
import type { TypedMolecule } from '../src/types';

describe('Bond Stretch Energy', () => {
  it('calculates energy for a simple 2-atom system', () => {
    // Using parameter: '1-1': { k_b: 4.258, r0: 1.508 }
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
    // E = 143.88 * 4.258 * (0.092)^2 = 5.1856
    expect(energy).toBeCloseTo(5.1856, 3);
  });

  it('handles reversed type keys (min-max ordering)', () => {
    // Assuming '1-2': { k_b: 4.539, r0: 1.482 } from mmffbond.par
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
    // E = 143.88 * 4.539 * (0.1)^2 = 6.5307
    expect(energy).toBeCloseTo(6.5307, 3);
  });
});
