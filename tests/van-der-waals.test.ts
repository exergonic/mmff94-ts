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

  it('equals −ε at r = R* for two identical type-5 H atoms', () => {
    // For type 5 (H):  R* = 4.2 Å, α = 0.25 Å³, N = 0.8, G = 1.209
    //
    // At r = R*, the buffered 14-7 simplifies to exactly E = −ε:
    //
    //   (1.07·R* / (R* + 0.07·R*))⁷ = (1.07 / 1.07)⁷ = 1
    //   (1.12·R*⁷ / (R*⁷ + 0.12·R*⁷) − 2) = (1.12 / 1.12) − 2 = −1
    //   ∴ E = ε × 1 × (−1) = −ε
    //
    // Slater-Kirkwood well depth for two identical type-5 H:
    //
    //   ε = 181.16 · G² · α² / (2 · α/√N)
    //     = 181.16 × 1.209² × 0.25² / (2 × 0.25/√0.8)
    //     = 181.16 × 1.4617 × 0.0625 / (2 × 0.27951)
    //     = 181.16 × 0.091356 / 0.55902
    //     = 29.61 kcal/mol
    //
    // So at r = R* = 4.2 Å, E must equal −29.61.

    const R_star = 4.2;
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'H', x: 0, y: 0, z: 0 },
        { index: 1, element: 'H', x: R_star, y: 0, z: 0 },
      ],
      bonds: [],
      atom_types: [5, 5],
    };

    const energy = calc_vdw_energy(mol);
    expect(energy).toBeCloseTo(-29.61, 1);
  });

  it('gives positive energy inside the repulsive wall (r < R*)', () => {
    // Two type-5 H atoms at 2.0 Å — well inside R* = 4.2 Å.
    // The buffered 14-7 should give a large positive repulsion.
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'H', x: 0, y: 0, z: 0 },
        { index: 1, element: 'H', x: 2.0, y: 0, z: 0 },
      ],
      bonds: [],
      atom_types: [5, 5],
    };

    const energy = calc_vdw_energy(mol);
    expect(energy).toBeGreaterThan(0);
  });

  it('approaches zero at large separation (r >> R*)', () => {
    // Two type-5 H atoms at 50 Å — the buffered 14-7 should be near zero.
    const mol: TypedMolecule = {
      atoms: [
        { index: 0, element: 'H', x: 0, y: 0, z: 0 },
        { index: 1, element: 'H', x: 50, y: 0, z: 0 },
      ],
      bonds: [],
      atom_types: [5, 5],
    };

    const energy = calc_vdw_energy(mol);
    // Should be a very small negative number
    expect(Math.abs(energy)).toBeLessThan(0.01);
  });
});
