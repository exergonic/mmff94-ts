import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import type { Molecule, TypedMolecule } from '../src/types.js';
import { parse_sdf } from '../src/sdf.js';
import { parameter_gap_report } from '../src/mmff94/parameter-gaps.js';

// parameter_gap_report: which atoms run on generic MMFF94 parameters.
// The signal is the type's expected coordination vs the atom's actual
// σ-neighbor count — Tinker's own "Unusual Number of Attached Atoms"
// check — plus the type-1 fallback for elements outside the MMFF94
// type space (the case other toolkits refuse outright).

function fixture(name: string): Molecule {
  return parse_sdf(readFileSync(`tests/fixtures/sdf/${name}`, 'utf-8'));
}

describe('parameter_gap_report', () => {
  it('is empty for a well-typed molecule (ethane)', () => {
    const report = parameter_gap_report(fixture('ethane.sdf'));
    expect(report.atoms).toEqual([]);
    expect(report.untyped).toEqual([]);
  });

  it('flags pentacoordinate P in PCl5 (type 25 expects crd 4, has 5)', () => {
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
    const report = parameter_gap_report(pcl5);
    expect(report.atoms).toEqual([
      { index: 0, element: 'P', type: 25, coordination: 5, expected_crd: 4 },
    ]);
    expect(report.untyped).toEqual([]);
  });

  it('flags hexacoordinate S in SF6 (type 15 expects crd 2, has 6)', () => {
    const sf6: Molecule = {
      atoms: [
        { index: 0, element: 'S', x: 0, y: 0, z: 0 },
        { index: 1, element: 'F', x: 1.56, y: 0, z: 0 },
        { index: 2, element: 'F', x: -1.56, y: 0, z: 0 },
        { index: 3, element: 'F', x: 0, y: 1.56, z: 0 },
        { index: 4, element: 'F', x: 0, y: -1.56, z: 0 },
        { index: 5, element: 'F', x: 0, y: 0, z: 1.56 },
        { index: 6, element: 'F', x: 0, y: 0, z: -1.56 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 0, atom2: 2, bond_order: 1 },
        { atom1: 0, atom2: 3, bond_order: 1 },
        { atom1: 0, atom2: 4, bond_order: 1 },
        { atom1: 0, atom2: 5, bond_order: 1 },
        { atom1: 0, atom2: 6, bond_order: 1 },
      ],
    };
    const report = parameter_gap_report(sf6);
    expect(report.atoms).toEqual([
      { index: 0, element: 'S', type: 15, coordination: 6, expected_crd: 2 },
    ]);
    expect(report.untyped).toEqual([]);
  });

  it('reports elements outside the MMFF94 type space as untyped', () => {
    // Al has no MMFF94 type: the typer's last resort is the generic
    // sp³ C type (1), which is what other toolkits refuse outright.
    const al: Molecule = {
      atoms: [{ index: 0, element: 'Al', x: 0, y: 0, z: 0 }],
      bonds: [],
    };
    const report = parameter_gap_report(al);
    expect(report.atoms).toEqual([]);
    expect(report.untyped).toEqual([0]);
  });

  it('is empty for the ylide (P type 25 crd 4, 4 σ — the crd-aware typing)', () => {
    const report = parameter_gap_report(fixture('methylenetriphenylphosphorane.sdf'));
    expect(report.atoms).toEqual([]);
    expect(report.untyped).toEqual([]);
  });

  it('does NOT flag fewer-neighbors-than-crd atoms (the validated unsaturated variants)', () => {
    // The FIZGEA pattern: an imide N typed 43 (crd 3) with only 2 σ
    // neighbors. The suite validates this chemistry exactly (the
    // type's angular model still fits), so the report must stay
    // silent — only the EXCEEDS direction is a gap.
    const imideN: TypedMolecule = {
      atoms: [
        { index: 0, element: 'N', x: 0, y: 0, z: 0 },
        { index: 1, element: 'S', x: 1.5, y: 0, z: 0 },
        { index: 2, element: 'C', x: -1.2, y: 0.9, z: 0 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 2 },
        { atom1: 0, atom2: 2, bond_order: 1 },
      ],
      atom_types: [43, 74, 3],
    };
    const report = parameter_gap_report(imideN);
    expect(report.atoms).toEqual([]);
    expect(report.untyped).toEqual([]);
  });
});
