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

import type { TypedMolecule, OptimizationResult } from '../types';
import type { EnergyGradientFn } from './l-bfgs';

export interface SteepestDescentOptions {
  max_iterations?: number;
  gradient_tolerance?: number; // kcal/mol/Å
  initial_step_size?: number;  // first line-search trial (default 1.0)
}

const C1 = 1e-4;            // Armijo sufficient-decrease constant
const MAX_LINE_SEARCH = 40; // backtracking budget for one line search
const MAX_TRIAL_STEP = 2.0; // Å — the first trial never moves an atom farther
const BACKTRACK = 0.5;      // halving factor on a rejected trial

/** Largest |force component| over all atoms (kcal/mol/Å). */
function max_gradient_norm(gradient: number[][]): number {
  let m = 0;
  for (const g of gradient) {
    for (const v of g) m = Math.max(m, Math.abs(v));
  }
  return m;
}

/**
 * Optimize geometry using steepest descent with an Armijo line search.
 *
 * Primarily a fallback for when L-BFGS fails to converge. The input
 * molecule is not mutated; the result carries the optimized geometry.
 */
export function optimize_steepest_descent(
  molecule: TypedMolecule,
  calc_energy_gradient: EnergyGradientFn,
  options?: SteepestDescentOptions,
): OptimizationResult {
  const max_iterations = options?.max_iterations ?? 1000;
  const gradient_tolerance = options?.gradient_tolerance ?? 0.05;
  const initial_step_size = options?.initial_step_size ?? 1.0;

  // Working copy — the optimizer owns its coordinates and never
  // mutates the caller's molecule (same pattern as L-BFGS: the spread
  // keeps atom_types and partial_charges; only the atoms are cloned).
  const work: TypedMolecule = {
    ...molecule,
    atoms: molecule.atoms.map(a => ({ ...a })),
  };

  let state = calc_energy_gradient(work);
  let iterations = 0;
  let converged = max_gradient_norm(state.gradient) < gradient_tolerance;

  while (!converged && iterations < max_iterations) {
    iterations++;

    const g = state.gradient;
    // The Armijo slope: ‖∇E‖² (the directional derivative along −g).
    let g2 = 0;
    for (const grad of g) {
      for (const v of grad) g2 += v * v;
    }
    // The first trial is capped in physical space: the 2-norm of the
    // step is at most MAX_TRIAL_STEP Å, so no atom moves farther than
    // that (smaller steps on steep slopes), and never larger than the
    // requested initial step size.
    let alpha = Math.min(initial_step_size, MAX_TRIAL_STEP / Math.sqrt(g2));

    let accepted = false;
    for (let ls = 0; ls < MAX_LINE_SEARCH; ls++) {
      const trial: TypedMolecule = {
        ...work,
        atoms: work.atoms.map((a, idx) => ({
          ...a,
          x: a.x - alpha * g[idx][0],
          y: a.y - alpha * g[idx][1],
          z: a.z - alpha * g[idx][2],
        })),
      };
      const trial_state = calc_energy_gradient(trial);
      if (trial_state.energy.total <= state.energy.total - C1 * alpha * g2) {
        // Sufficient decrease: accept, commit the trial as the new
        // working geometry (its atoms are already fresh copies).
        work.atoms = trial.atoms;
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

    converged = max_gradient_norm(state.gradient) < gradient_tolerance;
  }

  return {
    molecule: work,
    energy: state.energy,
    iterations,
    converged,
    final_max_gradient: max_gradient_norm(state.gradient),
  };
}
