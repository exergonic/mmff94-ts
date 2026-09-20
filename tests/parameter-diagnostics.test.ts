import { describe, it, expect, afterEach } from 'vitest';
import { readFileSync } from 'fs';
import type { Molecule } from '../src/types.js';
import { parse_sdf } from '../src/sdf.js';
import {
  diagnose_molecule,
  set_parameter_warning_handler,
  type ParameterDiagnostics,
} from '../src/mmff94/diagnostics.js';

// diagnose_molecule: the three things that can happen to an interaction —
// a stored row, the part V empirical rules, or a drop (no row and no
// rule). The first two are fine; the third is the one the library must
// never leave silent. These tests measure the counts on real molecules
// and check that anything not clean reaches the warning handler.

function fixture(name: string): Molecule {
  return parse_sdf(readFileSync(`tests/fixtures/sdf/${name}`, 'utf-8'));
}

function total(map: Record<string, number>): number {
  return Object.values(map).reduce((a, b) => a + b, 0);
}

afterEach(() => set_parameter_warning_handler((s) => console.warn(s)));

describe('diagnose_molecule', () => {
  it('reports a well-parameterized molecule as clean, with no empirical use', () => {
    const d = diagnose_molecule(fixture('ethane.sdf'), { quiet: true });
    expect(d.clean).toBe(true);
    expect(d.atoms).toEqual([]);
    expect(total(d.empirical)).toBe(0);
    expect(total(d.dropped)).toBe(0);
    expect(d.summary).toContain('clean');
  });

  it('counts empirical-rule interactions without calling them gaps', () => {
    // Vinyl phosphine exercises the part V rules: the phosphine angle
    // wildcard has k_a = 0 and the P-C-H force constant comes from
    // eq. 20 (see docs/vinyl-phosphine-divergence.md).
    const d = diagnose_molecule(fixture('vinylphosphine.sdf'), { quiet: true });
    console.log(`vinylphosphine: ${d.summary}`);
    expect(d.empirical.angle).toBeGreaterThan(0);
    expect(d.clean).toBe(true);
    expect(total(d.dropped)).toBe(0);
  });

  it('drops nothing in a 304-atom protein fragment (trp-cage)', () => {
    const d = diagnose_molecule(fixture('trpcage.sdf'), { quiet: true });
    console.log(`trpcage: ${d.summary}, empirical ${JSON.stringify(d.empirical)}`);
    expect(total(d.dropped)).toBe(0);
    expect(d.clean).toBe(true);
  });

  it('flags a coordination gap and warns through the handler', () => {
    const pcl5: Molecule = {
      atoms: [
        { index: 0, element: 'P', x: 0, y: 0, z: 0 },
        { index: 1, element: 'Cl', x: 0, y: 0, z: 2 },
        { index: 2, element: 'Cl', x: 0, y: 0, z: -2 },
        { index: 3, element: 'Cl', x: 2, y: 0, z: 0 },
        { index: 4, element: 'Cl', x: -1, y: 1.7, z: 0 },
        { index: 5, element: 'Cl', x: -1, y: -1.7, z: 0 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 0, atom2: 2, bond_order: 1 },
        { atom1: 0, atom2: 3, bond_order: 1 },
        { atom1: 0, atom2: 4, bond_order: 1 },
        { atom1: 0, atom2: 5, bond_order: 1 },
      ],
    };
    const seen: string[] = [];
    set_parameter_warning_handler((s) => seen.push(s));

    const d: ParameterDiagnostics = diagnose_molecule(pcl5);
    expect(d.clean).toBe(false);
    expect(d.atoms.length).toBe(1);
    expect(seen.length).toBe(1);
    expect(seen[0]).toContain('coordination');

    // quiet suppresses the handler but not the report
    const quiet = diagnose_molecule(pcl5, { quiet: true });
    expect(seen.length).toBe(1);
    expect(quiet.atoms.length).toBe(1);
  });

  it('distinguishes a rule-declined torsion from a parameter gap', () => {
    // Cyanamide H2N-C#N: the H-N-C#N paths run through a linear centre, so
    // MMFF94 rules (a)/(e)/(f) declare those torsions absent. The library
    // reports that as a note, not as a gap — and never silently.
    const cyanamide: Molecule = {
      atoms: [
        { index: 0, element: 'N', x: 0, y: 0, z: 0 },
        { index: 1, element: 'H', x: 0.8, y: 0.6, z: 0 },
        { index: 2, element: 'H', x: -0.8, y: 0.6, z: 0 },
        { index: 3, element: 'C', x: 1.35, y: 0, z: 0 },
        { index: 4, element: 'N', x: 2.5, y: 0, z: 0 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 0, atom2: 2, bond_order: 1 },
        { atom1: 0, atom2: 3, bond_order: 1 },
        { atom1: 3, atom2: 4, bond_order: 3 },
      ],
    };
    const seen: string[] = [];
    set_parameter_warning_handler((s) => seen.push(s));

    const d = diagnose_molecule(cyanamide);
    expect(d.dropped.torsion).toBeGreaterThan(0);
    expect(d.atoms).toEqual([]);
    expect(d.clean).toBe(false);
    expect(seen.length).toBe(1);
    expect(seen[0]).toContain('notes');
    expect(seen[0]).not.toContain('GAPS');
  });

  it('does not warn about a clean molecule', () => {
    const seen: string[] = [];
    set_parameter_warning_handler((s) => seen.push(s));
    diagnose_molecule(fixture('benzene.sdf'));
    expect(seen).toEqual([]);
  });
});
