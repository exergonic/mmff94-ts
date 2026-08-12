import type { Molecule, TypedMolecule } from '../types.js';
import { assign_atom_types } from './assign-atom-types.js';
import { ATOM_TYPE_PROPERTIES } from './parameters/index.js';

/**
 * One atom whose assigned MMFF94 type cannot represent its actual
 * coordination.
 *
 * MMFF94's type space has no entries for hypervalent main-group
 * centers: hexacoordinate S (SF₆) falls back to type 15 (crd 2 —
 * the thiol/sulfide type), pentacoordinate P (PCl₅) to type 25
 * (crd 4 — the PO4 family), and the empirical angle protocol then
 * emits reference angles from a coordination branch that cannot
 * produce the real geometry. The resulting structure is generic,
 * not MMFF94-quality. Tinker flags exactly this in its own output
 * ("Atoms with an Unusual Number of Attached Atoms: Expected 2,
 * Found 6"); this report is the same check.
 *
 * Only the EXCEEDS direction is flagged. Fewer neighbors than the
 * type's crd is the validated "unsaturated variant" pattern — the
 * imide N of FIZGEA (type 43, crd 3, 2 neighbors) and the sulfonyl
 * carbene S of SURDOX02 (type 18, crd 4, 3) both trip a plain
 * mismatch yet reproduce BatchMin's suite energies exactly, because
 * the type's angular model still fits. A coordination ABOVE the
 * type's crd means the type space has no representation of the
 * atom's environment at all (measured: 0 of the 761 suite
 * molecules exceed their type's crd — the signal never false-fires
 * on validated chemistry).
 */
export interface AtomParameterGap {
  index: number;
  element: string;
  type: number;
  coordination: number;   // actual σ neighbors in the graph
  expected_crd: number;   // the type's crd from the properties table
}

/**
 * Which atoms of a molecule run on generic MMFF94 parameters.
 *
 * - `atoms`: atoms whose coordination EXCEEDS their type's crd
 *   (hypervalent centers, and the metal-cation types whose crd is 0
 *   — the parameter-inert salts of the validation suite).
 * - `untyped`: atoms whose element is outside the MMFF94 type space
 *   and fell to the generic type-1 (sp³ C) fallback — the case other
 *   toolkits refuse outright (openchemlib cannot type hexavalent S
 *   and rejects the molecule; OpenBabel silently types it 0).
 *
 * The report is diagnostic only — it never changes typing or
 * energies. Callers decide what to do (e.g. warn the user that a
 * locally refined geometry is approximate).
 */
export interface ParameterGapReport {
  atoms: AtomParameterGap[];
  untyped: number[];
}

export function parameter_gap_report(molecule: Molecule | TypedMolecule): ParameterGapReport {
  const typed = 'atom_types' in molecule ? molecule : assign_atom_types(molecule);

  const neighbors = new Array<number>(typed.atoms.length).fill(0);
  for (const bond of typed.bonds) {
    neighbors[bond.atom1]++;
    neighbors[bond.atom2]++;
  }

  const atoms: AtomParameterGap[] = [];
  const untyped: number[] = [];

  for (let i = 0; i < typed.atoms.length; i++) {
    const type = typed.atom_types[i];
    const element = typed.atoms[i].element;

    // The typer's last resort for elements outside its switch and
    // metal map is the generic sp³ C type (1). A non-C atom carrying
    // type 1 got here through that fallback.
    if (type === 1 && element !== 'C') {
      untyped.push(i);
      continue;
    }

    const expected_crd = ATOM_TYPE_PROPERTIES[type]?.crd ?? -1;
    if (neighbors[i] > expected_crd) {
      atoms.push({ index: i, element, type, coordination: neighbors[i], expected_crd });
    }
  }

  return { atoms, untyped };
}
