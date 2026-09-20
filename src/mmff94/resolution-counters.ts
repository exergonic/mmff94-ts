/**
 * Resolution counters for diagnose_molecule().
 *
 * The energy terms resolve every interaction through the parameter tables
 * first and the part V empirical rules second; when both fail the
 * interaction is dropped. Those events are invisible in the energy, so
 * diagnose_molecule() arms these counters and runs one readable energy
 * evaluation, which enumerates exactly the interactions the field uses.
 *
 * Every note_* call checks the flag first, so the hot path pays nothing;
 * the counters are armed only during a diagnostic pass.
 */

export type InteractionKind =
  | 'bond'
  | 'angle'
  | 'stretch_bend'
  | 'torsion'
  | 'out_of_plane';

export interface ResolutionCounters {
  /** Interactions resolved by the part V empirical rules (legitimate). */
  empirical: Record<InteractionKind, number>;
  /** Interactions dropped: no table row and no empirical rule. */
  dropped: Record<InteractionKind, number>;
}

const KINDS: InteractionKind[] = ['bond', 'angle', 'stretch_bend', 'torsion', 'out_of_plane'];

let armed: ResolutionCounters | null = null;

function zeroed(): Record<InteractionKind, number> {
  const out = {} as Record<InteractionKind, number>;
  for (const k of KINDS) out[k] = 0;
  return out;
}

/** Arm the counters; returns the live object. Diagnostics only. */
export function arm_resolution_counters(): ResolutionCounters {
  armed = { empirical: zeroed(), dropped: zeroed() };
  return armed;
}

/** Disarm. Safe to call even when not armed. */
export function disarm_resolution_counters(): void {
  armed = null;
}

/** Called from a term's cold path when the empirical rules supplied the
 *  parameters. No-op unless a diagnostic pass is running. */
export function note_empirical(kind: InteractionKind): void {
  if (armed) armed.empirical[kind]++;
}

/** Called from a term's cold path when an interaction was dropped. No-op
 *  unless a diagnostic pass is running. */
export function note_dropped(kind: InteractionKind): void {
  if (armed) armed.dropped[kind]++;
}
