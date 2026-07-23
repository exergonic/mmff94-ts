import { describe, it, expect } from 'vitest';
import { calc_torsion_energy } from '../src/mmff94/energy/torsion';
import type { TypedMolecule } from '../src/types';

describe('Torsion Energy', () => {
  it('returns zero for a molecule with no dihedrals', () => {
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'C', x: 0, y: 0, z: 0 },
        { index: 1, element: 'C', x: 1.5, y: 0, z: 0 },
      ],
      bonds: [{ atom1: 0, atom2: 1, bond_order: 1 }],
      atom_types: [1, 1],
    };
    expect(calc_torsion_energy(mol)).toBeCloseTo(0.0, 10);
  });

  it('calculates energy for H-C-C-H dihedral at tau = 0° (eclipsed)', () => {
    // Parameter '0-5-1-1-5': V1=0.284 (n=1, gamma=0), V2=-1.386 (n=2, gamma=180), V3=0.314 (n=3, gamma=0)
    // tau = 0°:  E = (0.284/2)*(1+cos(0)) + (-1.386/2)*(1+cos(-180)) + (0.314/2)*(1+cos(0))
    //          = 0.142*2 + (-0.693)*(1+(-1)) + 0.157*2 = 0.284 + 0 + 0.314 = 0.598
    //
    // Geometry: both H atoms on the same side of the C-C axis gives τ = 0°.
    const rCC = 1.5;
    const rCH = 1.1;

    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'H', x: 0, y: rCH, z: 0 },
        { index: 1, element: 'C', x: 0, y: 0, z: 0 },
        { index: 2, element: 'C', x: rCC, y: 0, z: 0 },
        { index: 3, element: 'H', x: rCC, y: -rCH, z: 0 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 1, atom2: 2, bond_order: 1 },
        { atom1: 2, atom2: 3, bond_order: 1 },
      ],
      atom_types: [5, 1, 1, 5],
    };

    expect(calc_torsion_energy(mol)).toBeCloseTo(0.598, 3);
  });

  it('calculates energy for H-C-C-H dihedral at tau = 180° (staggered)', () => {
    // tau = 180°: E = (0.284/2)*(1+cos(180)) + (-1.386/2)*(1+cos(180)) + (0.314/2)*(1+cos(540))
    //           = 0.142*0 + (-0.693)*0 + 0.157*0 = 0.0
    //
    // Geometry: H atoms on opposite sides of the C-C axis gives τ = 180°.
    const rCC = 1.5;
    const rCH = 1.1;

    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'H', x: 0, y: rCH, z: 0 },
        { index: 1, element: 'C', x: 0, y: 0, z: 0 },
        { index: 2, element: 'C', x: rCC, y: 0, z: 0 },
        { index: 3, element: 'H', x: rCC, y: rCH, z: 0 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 1, atom2: 2, bond_order: 1 },
        { atom1: 2, atom2: 3, bond_order: 1 },
      ],
      atom_types: [5, 1, 1, 5],
    };

    expect(calc_torsion_energy(mol)).toBeCloseTo(0.0, 3);
  });
});
