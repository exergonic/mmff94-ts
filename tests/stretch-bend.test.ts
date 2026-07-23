import { describe, it, expect } from 'vitest';
import { calc_stretch_bend_energy } from '../src/mmff94/energy/stretch-bend';
import type { TypedMolecule } from '../src/types';

describe('Stretch-Bend Energy', () => {
  it('returns zero for equilibrium geometry', () => {
    // C-C-C angle at equilibrium (109.608°, parameter '0-1-1-1')
    // Bonds at equilibrium (1.508 Å, parameter '0-1-1')
    // Place central atom at origin, atom0 on +x axis, atom2 at angle theta from +x
    const r0 = 1.508;
    const theta_deg = 109.608;
    const theta_rad = theta_deg * Math.PI / 180;

    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'C', x: r0, y: 0, z: 0 },
        { index: 1, element: 'C', x: 0, y: 0, z: 0 },
        { index: 2, element: 'C', x: r0 * Math.cos(theta_rad), y: r0 * Math.sin(theta_rad), z: 0 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 1, atom2: 2, bond_order: 1 },
      ],
      atom_types: [1, 1, 1],
    };

    const energy = calc_stretch_bend_energy(mol);
    expect(energy).toBeCloseTo(0.0, 10);
  });

  it('calculates non-zero energy for distorted geometry', () => {
    // Both bonds stretched by 0.1 Å, angle closed by 10°
    const r0 = 1.508;
    const dr = 0.1;
    const theta0_deg = 109.608;
    const d_theta = -10.0;
    const actual_theta_deg = theta0_deg + d_theta;
    const actual_theta_rad = actual_theta_deg * Math.PI / 180;

    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'C', x: r0 + dr, y: 0, z: 0 },
        { index: 1, element: 'C', x: 0, y: 0, z: 0 },
        { index: 2, element: 'C', x: (r0 + dr) * Math.cos(actual_theta_rad), y: (r0 + dr) * Math.sin(actual_theta_rad), z: 0 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 1, atom2: 2, bond_order: 1 },
      ],
      atom_types: [1, 1, 1],
    };

    const energy = calc_stretch_bend_energy(mol);
    // Both bonds stretched by 0.1 Å, d_theta = -10°
    // Parameter '0-1-1-1': k_sb_IJK: 0.206, k_sb_KJI: 0.206
    // E per angle = 2.51210 * (0.206 * 0.1 + 0.206 * 0.1) * (-10)
    //            = 2.51210 * 0.0412 * (-10) = -1.0350
    expect(energy).toBeCloseTo(-1.035, 2);
  });
});
