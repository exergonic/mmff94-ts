import { describe, it, expect } from 'vitest';
import { calc_angle_bend_energy } from '../src/mmff94/energy/angle-bend';
import type { TypedMolecule } from '../src/types';

describe('Angle Bend Energy', () => {
  it('calculates energy for a simple 3-atom system (eq. 3, harmonic + cubic)', () => {
    // Parameter '0-1-1-1': { k_a: 0.851, theta0: 109.608 }
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
    // θ = 90°, θ₀ = 109.608°, Δθ = −19.608°
    // k_a = 0.851, half_k_a = 0.4255, cb = −0.007
    // Harmonic:  0.043844 * 0.4255 * (−19.608)² = 7.1726
    // Anharmonic: 1 + (−0.007)(−19.608) = 1.1373
    // E = 7.1726 * 1.1373 = 8.1571
    expect(energy).toBeCloseTo(8.1571, 3);
  });

  it('uses cosine form for near-linear angles (eq. 4)', () => {
    // Acetylenic C angles use eq. (4): E = 143.9325 * k_a * (1 + cos θ)
    // Parameter '0-1-4-4': { k_a: 0.423, theta0: 180.0 }
    // At θ = 180° (collinear), cos 180° = −1, so E = 0
    // At θ = 150°, cos 150° = −√3/2 ≈ −0.8660
    //   E = 143.9325 * 0.423 * (1 − 0.8660) = 8.156
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'C', x: 0, y: 0, z: 1.2 },
        { index: 1, element: 'C', x: 0, y: 0, z: 0 },
        { index: 2, element: 'C', x: 0, y: 0, z: -1.2 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 1, atom2: 2, bond_order: 3 },
      ],
      // Parameter lookup: t_min=1, tj=4, t_max=4 → key "0-1-4-4"
      atom_types: [1, 4, 4],
    };

    // Collinear: θ = 180°, energy is zero
    expect(calc_angle_bend_energy(mol)).toBeCloseTo(0.0, 4);

    // Bent at 150°
    const sin150 = Math.sin(150 * Math.PI / 180);
    const cos150 = Math.cos(150 * Math.PI / 180);
    const bent: TypedMolecule = {
      ...mol,
      atoms: [
        { index: 0, element: 'C', x: 0, y: sin150, z: cos150 },
        { index: 1, element: 'C', x: 0, y: 0, z: 0 },
        { index: 2, element: 'C', x: 0, y: 0, z: -1 },
      ],
    };
    const energy = calc_angle_bend_energy(bent);
    expect(energy).toBeGreaterThan(0);
  });
});
