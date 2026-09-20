/**
 * Parameter diagnostics: what the force field had to work around for a
 * molecule, in one report.
 *
 * Three things can happen to an interaction, and only the first is
 * silent in the energy alone:
 *
 *   1. a stored parameter row is used (the normal path);
 *   2. no row exists and the part V empirical rules build the parameters
 *      (legitimate — Halgren designed it that way, and this library
 *      matches the reference behaviour);
 *   3. neither exists and the interaction is DROPPED from the energy.
 *
 * diagnose_molecule() measures (2) and (3) by arming the resolution
 * counters and running one readable energy evaluation — that evaluation
 * enumerates exactly the interactions the field uses, so the counts are
 * what the energy really did, not a re-derivation. It costs one extra
 * energy evaluation, and only when called; the hot path is untouched.
 *
 * Coordination gaps (an atom whose type cannot represent its
 * coordination — hypervalent centres such as SF6 or PCl5) come from
 * parameter_gap_report(), which this wraps.
 *
 * Anything not clean is reported through the warning handler
 * (set_parameter_warning_handler, console.warn by default). Loud by
 * design: the library would rather say "this interaction was dropped"
 * than return a quietly incomplete energy.
 */
import { calc_energy } from './energy/total.js';
import { parameter_gap_report, type AtomParameterGap } from './parameter-gaps.js';
import {
  arm_resolution_counters,
  disarm_resolution_counters,
  type InteractionKind,
  type ResolutionCounters,
} from './resolution-counters.js';
import type { Molecule, TypedMolecule } from '../types.js';

const KINDS_ORDER: InteractionKind[] = [
  'bond',
  'angle',
  'stretch_bend',
  'torsion',
  'out_of_plane',
];

export interface ParameterDiagnostics extends ResolutionCounters {
  /** Atoms whose type cannot represent their coordination. */
  atoms: AtomParameterGap[];
  /** True when no interaction was omitted and no atom exceeds its type's
   *  coordination. Empirical-rule use alone keeps a molecule clean;
   *  rule-declined torsion paths (see dropped.torsion) do not. */
  clean: boolean;
  /** One-line human summary, also what the warning handler receives. */
  summary: string;
}

export type ParameterWarningHandler = (summary: string) => void;

let handler: ParameterWarningHandler | null =
  typeof console !== 'undefined' && typeof console.warn === 'function'
    ? (s) => console.warn(s)
    : null;

/** Install (or clear, with null) the handler that receives the summary of
 *  any molecule that is not clean. Defaults to console.warn. */
export function set_parameter_warning_handler(next: ParameterWarningHandler | null): void {
  handler = next;
}

function count(map: Record<InteractionKind, number>): number {
  return KINDS_ORDER.reduce((a, k) => a + map[k], 0);
}

export function diagnose_molecule(
  molecule: Molecule | TypedMolecule,
  options: { quiet?: boolean } = {},
): ParameterDiagnostics {
  // The live object; the note_* calls in the terms mutate it during the
  // evaluation below.
  const counters = arm_resolution_counters();
  try {
    calc_energy(molecule as Molecule);
  } finally {
    disarm_resolution_counters();
  }
  const gap = parameter_gap_report(molecule);

  const empiricalTotal = count(counters.empirical);
  const droppedTotal = count(counters.dropped);
  // Torsion paths are the one kind the empirical rules can decline outright
  // (rules (a)/(e)/(f): linear centres and unsaturated-sp2 paths carry no
  // torsion). That is MMFF94 reporting the interaction does not exist, not a
  // missing parameter — so the two classes read differently below.
  const torsionOmitted = counters.dropped.torsion;
  const strictDropped = droppedTotal - torsionOmitted;
  const strictGap = gap.atoms.length > 0 || strictDropped > 0;
  const clean = gap.atoms.length === 0 && droppedTotal === 0;

  const parts: string[] = [];
  if (gap.atoms.length > 0) {
    const list = gap.atoms
      .slice(0, 5)
      .map((a) => `${a.element}${a.index + 1} (type ${a.type}, ${a.coordination} neighbours)`)
      .join(', ');
    parts.push(
      `${gap.atoms.length} atom(s) exceed their type's coordination: ${list}` +
        (gap.atoms.length > 5 ? ', …' : ''),
    );
  }
  if (torsionOmitted > 0) {
    parts.push(
      `${torsionOmitted} torsion path(s) omitted — MMFF94 rules (a)/(e)/(f): linear centres and unsaturated-sp2 paths carry no torsion`,
    );
  }
  if (strictDropped > 0) {
    const per = KINDS_ORDER.filter((k) => k !== 'torsion' && counters.dropped[k] > 0)
      .map((k) => `${counters.dropped[k]} ${k}`)
      .join(', ');
    parts.push(
      `${strictDropped} interaction(s) DROPPED — no parameter row and no empirical rule: ${per}`,
    );
  }
  const summary =
    parts.length === 0
      ? `parameters: clean${
          empiricalTotal > 0
            ? ` (${empiricalTotal} interaction(s) built by the part V empirical rules)`
            : ''
        }`
      : `${strictGap ? 'MMFF94 parameter GAPS' : 'MMFF94 parameter notes (no gaps)'} — ${parts.join('; ')}`;

  if (parts.length > 0 && !options.quiet && handler) handler(summary);

  return { ...counters, atoms: gap.atoms, clean, summary };
}
