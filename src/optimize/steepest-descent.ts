/**
 * Steepest descent optimizer — simple fallback minimizer.
 *
 * L-BFGS is the primary optimizer; this is the robust-but-slow
 * fallback for when L-BFGS misbehaves (e.g. a pathological starting
 * geometry). Each iteration moves the atoms in the direction of the
 * negative gradient, with the step size set by an Armijo backtracking
 * line search:
 *
 *   x_{k+1} = x_k − α·∇E(x_k)
 *
 *   E(x_k − α·∇E) < E(x_k) − c₁·α·‖∇E‖²   (c₁ = 1e-4, Armijo)
 *
 * Start with α = initial_step_size and halve it until the condition
 * holds. The first trial is capped in physical space (no atom moves
 * more than MAX_TRIAL_STEP = 2 Å) — the same wall-guard that protects
 * L-BFGS's first trial on MMFF94's stiff vdW surfaces.
 *
 * Armijo backtracking guarantees monotone descent, but the method is
 * only linearly convergent and zig-zags down valleys (MMFF94 surfaces
 * are stiff — condition ~10⁴), so it needs far more iterations than
 * L-BFGS to reach the same tolerance. Use it as a fallback, not as
 * the default.
 */

import type { Molecule, TypedMolecule, EnergyComponents, OptimizationResult } from '../types.js';
import type { EnergyGradientFn } from './l-bfgs.js';
import { prepare_molecule } from '../mmff94/prepare.js';
import { calc_energy } from '../mmff94/energy/total.js';
import { calc_gradient } from '../mmff94/gradient/total.js';
import { create_fast_system, FastSystem } from './fast-system.js';

export interface SteepestDescentOptions {
  max_iterations?: number;
  /** Gradient threshold (kcal/mol/Å) — meaning depends on `criterion`
   *  (see LbfgsOptions: 'max' bounds max |g_i|, 'rms' bounds the RMS
   *  gradient). Default 0.05. */
  gradient_tolerance?: number;
  /** Convergence criterion — 'max' | 'rms' | 'either' (default 'either',
   *  see LbfgsOptions). */
  criterion?: 'max' | 'rms' | 'either';
  /** RMS threshold for criterion 'either' (default 0.02 — see
   *  LbfgsOptions). */
  rms_gradient_tolerance?: number;
  initial_step_size?: number;  // first line-search trial (default 1.0)
}

const C1 = 1e-4;            // Armijo sufficient-decrease constant
const MAX_LINE_SEARCH = 40; // backtracking budget for one line search
const MAX_TRIAL_STEP = 2.0; // Å — the first trial never moves an atom farther
const BACKTRACK = 0.5;      // halving factor on a rejected trial

/** Largest |force component| over all atoms (kcal/mol/Å). */
function max_gradient_norm(gradient: ArrayLike<number>): number {
  let m = 0;
  for (let i = 0; i < gradient.length; i++) m = Math.max(m, Math.abs(gradient[i]));
  return m;
}

/** RMS of the force components (kcal/mol/Å) — the TINKER-style signal. */
function rms_gradient_norm(gradient: ArrayLike<number>): number {
  let s = 0;
  for (let i = 0; i < gradient.length; i++) s += gradient[i] * gradient[i];
  return gradient.length > 0 ? Math.sqrt(s / gradient.length) : NaN;
}

/**
 * Optimize geometry using steepest descent with an Armijo line search.
 *
 * Primarily a fallback for when L-BFGS fails to converge. The input
 * molecule is not mutated; the result carries the optimized geometry.
 * A bare Molecule is typed and charged on demand, and the
 * energy-and-gradient oracle defaults to the built-in one — the
 * simple path is a single call: optimize_steepest_descent(molecule).
 */
export function optimize_steepest_descent(
  molecule: Molecule,
  calc_energy_gradient_or_options?: EnergyGradientFn | SteepestDescentOptions,
  options?: SteepestDescentOptions,
): OptimizationResult {
  // The callback is optional, so the second slot accepts either the
  // oracle or the options (same convention as L-BFGS).
  const has_oracle = typeof calc_energy_gradient_or_options === 'function';
  const calc_energy_gradient = has_oracle ? calc_energy_gradient_or_options : undefined;
  const opts = has_oracle ? options : (calc_energy_gradient_or_options ?? options);
  const max_iterations = opts?.max_iterations ?? 1000;
  const gradient_tolerance = opts?.gradient_tolerance ?? 0.05;
  const criterion = opts?.criterion ?? 'either';
  const rms_effective =
    criterion === 'rms' ? (opts?.rms_gradient_tolerance ?? gradient_tolerance)
    : criterion === 'either' ? (opts?.rms_gradient_tolerance ?? 0.02)
    : undefined;
  const initial_step_size = opts?.initial_step_size ?? 1.0;

  /** The active convergence gates — mirrors L-BFGS's converged_now. */
  function converged_at(g: Float64Array): boolean {
    const mx = max_gradient_norm(g);
    const rms = rms_gradient_norm(g);
    if (criterion === 'max') return mx < gradient_tolerance;
    if (criterion === 'rms') return rms < rms_effective!;
    return mx < gradient_tolerance || rms < rms_effective!;
  }

  // Simple path: type + charge on demand, built-in oracle by default
  // (same as L-BFGS). The built-in path runs on the compiled fast
  // system — never materializing atom objects for trials (the generic
  // path below clones the atoms per trial; the fast path evaluates the
  // flat coordinate buffer directly).
  const prepared = prepare_molecule(molecule);
  const oracle =
    calc_energy_gradient ??
    ((m: TypedMolecule) => ({ energy: calc_energy(m), gradient: calc_gradient(m) }));
  const fast: FastSystem | null = calc_energy_gradient ? null : create_fast_system(prepared);

  // Working copy — the optimizer owns its coordinates and never
  // mutates the caller's molecule (same pattern as L-BFGS: the spread
  // keeps atom_types and partial_charges; only the atoms are cloned).
  const work: TypedMolecule = {
    ...prepared,
    atoms: prepared.atoms.map(a => ({ ...a })),
  };
  const n = 3 * work.atoms.length;
  const x: Float64Array = new Float64Array(n);
  for (let a = 0; a < work.atoms.length; a++) {
    x[3 * a] = work.atoms[a].x;
    x[3 * a + 1] = work.atoms[a].y;
    x[3 * a + 2] = work.atoms[a].z;
  }

  /** Flat energy+gradient: fast path directly, generic path via work. */
  function evaluate_at_coords(coords: Float64Array): {
    total: number; grad: Float64Array; components: EnergyComponents;
  } {
    if (fast) {
      const grad = new Float64Array(n);
      fast.evaluate(coords, grad);
      return { total: fast.total, grad, components: { ...fast.components } };
    }
    for (let a = 0; a < work.atoms.length; a++) {
      work.atoms[a].x = coords[3 * a];
      work.atoms[a].y = coords[3 * a + 1];
      work.atoms[a].z = coords[3 * a + 2];
    }
    const s = oracle(work);
    const grad = new Float64Array(n);
    for (let a = 0; a < work.atoms.length; a++) {
      grad[3 * a] = s.gradient[a][0];
      grad[3 * a + 1] = s.gradient[a][1];
      grad[3 * a + 2] = s.gradient[a][2];
    }
    return { total: s.energy.total, grad, components: s.energy };
  }

  let state = evaluate_at_coords(x);
  let iterations = 0;
  let converged = converged_at(state.grad);

  while (!converged && iterations < max_iterations) {
    iterations++;

    const g = state.grad;
    // The Armijo slope: ‖∇E‖² (the directional derivative along −g).
    let g2 = 0;
    for (let i = 0; i < n; i++) g2 += g[i] * g[i];
    // The first trial is capped in physical space: the 2-norm of the
    // step is at most MAX_TRIAL_STEP Å, so no atom moves farther than
    // that (smaller steps on steep slopes), and never larger than the
    // requested initial step size.
    let alpha = Math.min(initial_step_size, MAX_TRIAL_STEP / Math.sqrt(g2));

    // Trial buffer — reused across backtracking scans (no per-trial
    // atom cloning on the fast path).
    const trial: Float64Array = new Float64Array(n);
    let accepted = false;
    for (let ls = 0; ls < MAX_LINE_SEARCH; ls++) {
      for (let i = 0; i < n; i++) trial[i] = x[i] - alpha * g[i];
      const trial_state = evaluate_at_coords(trial);
      if (trial_state.total <= state.total - C1 * alpha * g2) {
        // Sufficient decrease: accept, commit the trial coords.
        for (let i = 0; i < n; i++) x[i] = trial[i];
        state = trial_state;
        accepted = true;
        break;
      }
      alpha *= BACKTRACK;
    }

    if (!accepted) {
      // No Armijo-acceptable step (e.g. a degenerate flat direction):
      // give up rather than oscillate. The geometry is left at the
      // best point found.
      break;
    }

    converged = converged_at(state.grad);
  }

  // Sync the flat coordinates back into the returned molecule.
  for (let a = 0; a < work.atoms.length; a++) {
    work.atoms[a].x = x[3 * a];
    work.atoms[a].y = x[3 * a + 1];
    work.atoms[a].z = x[3 * a + 2];
  }

  return {
    molecule: work,
    energy: state.components,
    iterations,
    converged,
    final_max_gradient: max_gradient_norm(state.grad),
    final_rms_gradient: rms_gradient_norm(state.grad),
  };
}
