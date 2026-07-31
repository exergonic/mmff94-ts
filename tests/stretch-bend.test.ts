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

  it('calculates energy for an asymmetric H-C-C angle', () => {
    // Regression test: the bond r₀ lookups must use each bond's OWN type
    // pair (sorted). The old code looked up [t_min, tj] and [t_max, tj],
    // which missed for H-C (types [5,1] stored as '0-1-5') and silently
    // skipped the angle — ethane's stretch-bend came out 0.0000.
    //
    // Parameter '0-1-1-5': k_sb_IJK = 0.227 (min-type side = C, i.e. the
    // C-C bond), k_sb_KJI = 0.070 (max-type side = H, the C-H bond).
    // θ₀ = 110.549° for H-C-C (not the 109.608° of C-C-C).
    // H at 1.193 Å (dr = +0.1), C-C at 1.608 Å (dr = +0.1), θ = 120.549°
    // (dθ = +10°).
    // E = 2.51210 * (0.070*0.1 + 0.227*0.1) * 10 = 2.51210 * 0.0297 * 10
    //   = 0.7461
    const rCH = 1.093 + 0.1;
    const rCC = 1.508 + 0.1;
    const theta_deg = 110.549 + 10.0;
    const theta_rad = theta_deg * Math.PI / 180;

    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'H', x: rCH * Math.cos(theta_rad), y: rCH * Math.sin(theta_rad), z: 0 },
        { index: 1, element: 'C', x: 0, y: 0, z: 0 },
        { index: 2, element: 'C', x: rCC, y: 0, z: 0 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 1, atom2: 2, bond_order: 1 },
      ],
      atom_types: [5, 1, 1],
    };

    const energy = calc_stretch_bend_energy(mol);
    expect(energy).toBeCloseTo(0.7461, 3);
  });
});
