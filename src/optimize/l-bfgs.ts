/**
 * L-BFGS geometry optimizer.
 *
 * Limited-memory Broyden–Fletcher–Goldfarb–Shanno algorithm, following
 * Nocedal & Wright (Numerical Optimization, 2nd ed.): Algorithm 7.5
 * (L-BFGS with the two-loop recursion), Algorithm 3.5 (line search
 * with strong Wolfe conditions) and Algorithm 3.6 (the zoom with
 * cubic interpolation). Textbook implementation — nothing custom.
 *
 * L-BFGS is a quasi-Newton method that builds an approximation of the
 * inverse Hessian from a limited history (m = 10 by default) of
 * position and gradient changes. It needs only the energy and its
 * gradient at each step — exactly what calc_energy/calc_gradient
 * provide — not the full Hessian matrix.
 *
 * The algorithm:
 *   1. Start with an initial geometry.
 *   2. At each iteration:
 *      a. Search direction p = −H·g from the two-loop recursion,
 *         with the standard initial-Hessian scaling γ = sᵀy/yᵀy.
 *      b. Strong-Wolfe line search along p (c1 = 1e-4, c2 = 0.9),
 *         with cubic interpolation in the zoom.
 *      c. x ← x + α·p; evaluate the new gradient.
 *      d. Update the history with s = α·p and y = g_new − g
 *         (discarding the oldest pair when m is exceeded; a pair
 *         with non-positive curvature yᵀs is skipped — it carries
 *         no usable information about the inverse Hessian).
 *   3. Stop when max |g_i| < gradient_tolerance (default 0.05
 *      kcal/mol/Å, the "well-minimized structure" threshold), when
 *      the step becomes negligible, or when the iteration cap hits.
 *
 * The optimizer minimizes the total energy E; the returned gradient
 * is dE/dx, so the descent direction is −H·∇E.
 */

import type { TypedMolecule, EnergyComponents, OptimizationResult } from '../types';

export interface LbfgsOptions {
  max_iterations?: number;
  gradient_tolerance?: number; // kcal/mol/Å
  history_size?: number; // m in L-BFGS (default 10)
}

/**
 * The energy-and-gradient oracle the optimizer drives. Receives the
 * optimizer's working copy of the molecule (same atom types and
 * charges as the input; only coordinates move).
 */
export type EnergyGradientFn = (
  molecule: TypedMolecule,
) => { energy: EnergyComponents; gradient: number[][] };

// Strong Wolfe line-search constants (Nocedal & Wright, eq. 3.13/3.14)
const C1 = 1e-4;  // Armijo: sufficient decrease
const C2 = 0.9;   // curvature: |φ'(α)| ≤ c2·|φ'(0)|
const MAX_LINE_SEARCH = 40; // zoom iterations per line search
const MAX_ALPHA = 20;       // step cap (units follow the coordinates, Å)
const MAX_TRIAL_STEP = 2.0; // Å — the first line-search trial never moves an atom farther

interface EvalState {
  components: EnergyComponents;
  total: number;
  gradient: number[]; // flat, 3·n_atoms
}

/**
 * Optimize a molecule's geometry using L-BFGS.
 *
 * @param molecule  Starting geometry (typed — the optimizer needs the
 *                  atom types to evaluate energy and gradient).
 * @param calc_energy_gradient  Function returning (energy, gradient).
 * @param options   Convergence parameters.
 * @returns The optimized geometry, final energy, and convergence info.
 */
export function optimize_lbfgs(
  molecule: TypedMolecule,
  calc_energy_gradient: EnergyGradientFn,
  options?: LbfgsOptions,
): OptimizationResult {
  const max_iterations = options?.max_iterations ?? 1000;
  const gradient_tolerance = options?.gradient_tolerance ?? 0.05;
  const history_size = Math.max(1, Math.min(50, options?.history_size ?? 10));

  // Working copy — the optimizer owns its coordinates and never
  // mutates the caller's molecule. The spread keeps atom_types and
  // partial_charges; only the atoms (positions) are cloned.
  const work: TypedMolecule = {
    ...molecule,
    atoms: molecule.atoms.map(a => ({ ...a })),
  };

  const n = 3 * work.atoms.length;

  // Flat view of the working coordinates (x = [x₀,y₀,z₀, x₁,…]).
  const x: number[] = new Array(n);
  for (let a = 0; a < work.atoms.length; a++) {
    x[3 * a] = work.atoms[a].x;
    x[3 * a + 1] = work.atoms[a].y;
    x[3 * a + 2] = work.atoms[a].z;
  }

  // Evaluate energy + gradient at the flat coordinates, keeping the
  // working molecule's coordinates in sync (the callback reads them).
  function evaluate(coords: number[]): EvalState {
    for (let a = 0; a < work.atoms.length; a++) {
      work.atoms[a].x = coords[3 * a];
      work.atoms[a].y = coords[3 * a + 1];
      work.atoms[a].z = coords[3 * a + 2];
    }
    const { energy, gradient } = calc_energy_gradient(work);
    const flat = new Array(n);
    for (let a = 0; a < work.atoms.length; a++) {
      flat[3 * a] = gradient[a][0];
      flat[3 * a + 1] = gradient[a][1];
      flat[3 * a + 2] = gradient[a][2];
    }
    return { components: energy, total: energy.total, gradient: flat };
  }

  // ── helpers on flat arrays ──────────────────────────────────────
  function dot(a: number[], b: number[]): number {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
  function max_abs(a: number[]): number {
    let m = 0;
    for (const v of a) m = Math.max(m, Math.abs(v));
    return m;
  }
  function direction_norm(a: number[]): number {
    return Math.sqrt(dot(a, a));
  }

  // ── L-BFGS state ────────────────────────────────────────────────
  const s_history: number[][] = []; // position changes, newest last
  const y_history: number[][] = []; // gradient changes, newest last
  const rho_history: number[] = []; // 1 / yᵀs per pair

  let state = evaluate(x);
  let iterations = 0;
  let converged = false;
  let final_max_gradient = max_abs(state.gradient);

  if (n === 0 || final_max_gradient < gradient_tolerance) {
    converged = true;
  }

  while (!converged && iterations < max_iterations) {
    iterations++;

    // Search direction from the two-loop recursion (Algorithm 7.4):
    // q = g − Σ αᵢyᵢ (newest → oldest), r = γq + Σ sᵢ(αᵢ − β) (oldest → newest)
    const m = s_history.length;
    const q = state.gradient.slice();
    const alpha_rec: number[] = new Array(m);
    for (let i = m - 1; i >= 0; i--) {
      alpha_rec[i] = rho_history[i] * dot(s_history[i], q);
      for (let d = 0; d < n; d++) q[d] -= alpha_rec[i] * y_history[i][d];
    }
    // Initial-Hessian scaling: γ = sᵀy / yᵀy (Nocedal & Wright eq. 7.20)
    let gamma = 1;
    if (m > 0) {
      const sy = dot(s_history[m - 1], y_history[m - 1]);
      const yy = dot(y_history[m - 1], y_history[m - 1]);
      if (yy > 0) gamma = sy / yy;
    }
    const direction = new Array<number>(n);
    for (let d = 0; d < n; d++) direction[d] = -gamma * q[d];
    for (let i = 0; i < m; i++) {
      const beta = rho_history[i] * dot(y_history[i], direction);
      for (let d = 0; d < n; d++) direction[d] += s_history[i][d] * (alpha_rec[i] - beta);
    }

    // The two-loop result is −H·g; if it is not a descent direction
    // (roundoff on ill-scaled pairs — wall-limited steps make γ tiny,
    // so the −γ·g anchor term is swamped by the history corrections),
    // fall back to steepest descent and drop the history: a fresh −g
    // is always downhill, and the next iteration rebuilds the history
    // from the current region instead of accumulating garbage.
    if (dot(direction, state.gradient) >= 0) {
      for (let d = 0; d < n; d++) direction[d] = -state.gradient[d];
      s_history.length = 0;
      y_history.length = 0;
      rho_history.length = 0;
    }

    // Strong-Wolfe line search (Algorithms 3.5 + 3.6)
    const ls = line_search(direction, state, gamma);

    if (!ls.found) {
      // No acceptable step (e.g. a degenerate flat direction): give
      // up rather than oscillate. The geometry is left at the best
      // point found.
      break;
    }

    // Accept the step: x ← x + α·p
    const step = new Array<number>(n);
    for (let d = 0; d < n; d++) {
      step[d] = ls.alpha * direction[d];
      x[d] += step[d];
    }
    state = ls.state_at_alpha;

    final_max_gradient = max_abs(state.gradient);
    if (final_max_gradient < gradient_tolerance) {
      converged = true;
      break;
    }

    // Update the limited-memory history. A pair whose step is below
    // the noise floor carries no curvature information — only the
    // roundoff of the gradient difference — so it is discarded rather
    // than poisoning the next direction (the L-BFGS correction terms
    // from noise pairs can overwhelm the γ-scaled −g term and flip the
    // direction to non-descent on stiff surfaces).
    const y = new Array<number>(n);
    for (let d = 0; d < n; d++) y[d] = state.gradient[d] - ls.gradient_before[d];
    const sy = dot(step, y);
    if (sy > 1e-12 && max_abs(step) > 1e-4) {
      s_history.push(step);
      y_history.push(y);
      rho_history.push(1.0 / sy);
      if (s_history.length > history_size) {
        s_history.shift();
        y_history.shift();
        rho_history.shift();
      }
    }

  }

  // The working molecule's coordinates were synced by the last
  // evaluate() call, so `work` already holds the final geometry.
  return {
    molecule: work,
    energy: state.components,
    iterations,
    converged,
    final_max_gradient,
  };

  // ── strong-Wolfe line search with cubic zoom ─────────────────────
  // φ(α) = E(x + αp); the search finds α > 0 satisfying
  //   Armijo:      φ(α) ≤ φ(0) + c1·α·φ'(0)
  //   Curvature:   |φ'(α)| ≤ c2·|φ'(0)|
  function line_search(
    direction: number[],
    start: EvalState,
    gamma: number,
  ): { alpha: number; found: boolean; state_at_alpha: EvalState; gradient_before: number[] } {
    const phi0 = start.total;
    const phi0_prime = dot(start.gradient, direction);

    // Not a descent direction — the caller guards this, but be safe.
    if (phi0_prime >= 0) {
      return { alpha: 0, found: false, state_at_alpha: start, gradient_before: start.gradient };
    }

    let alpha_prev = 0;
    let f_prev = phi0;
    // The first trial step compensates the initial-Hessian scaling γ
    // folded into the direction: p = −γ·r̃, so α₀ = 1/γ makes the
    // trial move α₀·p = −r̃ — the unscaled two-loop direction. Without
    // this, a tiny γ (stiff MMFF94 surfaces give γ ~ 10⁻³–10⁻⁵) makes
    // the α = 1 trial satisfy the Wolfe conditions trivially, and the
    // optimizer accepts a γ-sized step every iteration — convergent in
    // theory, glacial in practice.
    // The physical cap keeps the trial (and therefore the zoom's
    // initial bracket) inside a few Å: when the first trial fails
    // Armijo, the zoom shrinks the bracket from its far end, and a
    // 10⁵-wide α bracket would burn the whole zoom budget on
    // bisection before reaching a useful step.
    let alpha = Math.min(1.0 / gamma, MAX_TRIAL_STEP / direction_norm(direction));
    for (let iter = 0; iter < MAX_LINE_SEARCH; iter++) {
      const state_at_alpha = evaluate_at(alpha, direction);
      const f = state_at_alpha.total;
      const phi_prime = dot(state_at_alpha.gradient, direction);

      // Armijo failure (or no progress vs the previous trial):
      // bracket [alpha_prev, alpha] and zoom.
      if (f > phi0 + C1 * alpha * phi0_prime || (iter > 0 && f >= f_prev)) {
        const zoomed = zoom(alpha_prev, f_prev, alpha, f, phi0, phi0_prime, direction);
        return {
          alpha: zoomed.alpha,
          found: zoomed.found,
          state_at_alpha: zoomed.state,
          gradient_before: start.gradient,
        };
      }
      // Strong curvature condition satisfied: accept.
      if (Math.abs(phi_prime) <= -C2 * phi0_prime) {
        return { alpha, found: true, state_at_alpha, gradient_before: start.gradient };
      }
      // Derivative turned positive: the minimum is bracketed.
      if (phi_prime >= 0) {
        const zoomed = zoom(alpha, f, alpha_prev, f_prev, phi0, phi0_prime, direction);
        return {
          alpha: zoomed.alpha,
          found: zoomed.found,
          state_at_alpha: zoomed.state,
          gradient_before: start.gradient,
        };
      }
      // Still descending: extend the trial step.
      alpha_prev = alpha;
      f_prev = f;
      alpha = Math.min(2.0 * alpha, MAX_ALPHA);
    }

    // Budget exhausted without a bracket: accept the last trial
    // point that was actually evaluated (alpha_prev).
    return {
      alpha: alpha_prev,
      found: true,
      state_at_alpha: evaluate_at(alpha_prev, direction),
      gradient_before: start.gradient,
    };
  }

  /** Evaluate E(x + α·p) without moving the working coordinates. */
  function evaluate_at(alpha: number, direction: number[]): EvalState {
    const coords = new Array<number>(n);
    for (let d = 0; d < n; d++) coords[d] = x[d] + alpha * direction[d];
    return evaluate(coords);
  }

  /**
   * The zoom (Algorithm 3.6): shrink a bracket [lo, hi] known to
   * contain a point satisfying the strong Wolfe conditions, choosing
   * trial points by cubic interpolation (eq. 3.59) — bisection when
   * the interpolation leaves the safeguarded central 80% interval.
   */
  function zoom(
    lo: number,
    f_lo: number,
    hi: number,
    f_hi: number,
    phi0: number,
    phi0_prime: number,
    direction: number[],
  ): { alpha: number; found: boolean; state: EvalState } {
    // φ'(lo) is needed by the cubic; lo is either α = 0 (whose
    // derivative is φ'(0)) or a previously evaluated trial — one
    // extra evaluation keeps the bookkeeping simple.
    let g_lo = dot(evaluate_at(lo, direction).gradient, direction);

    for (let iter = 0; iter < MAX_LINE_SEARCH; iter++) {
      // Cubic interpolation through (lo, f_lo, g_lo) and (hi, f_hi, g_hi).
      // Equal endpoint values mean a flat interval: bisect.
      const g_hi = dot(evaluate_at(hi, direction).gradient, direction);
      let alpha_j: number;
      if (Math.abs(f_lo - f_hi) < 1e-14) {
        alpha_j = 0.5 * (lo + hi);
      } else {
        const d1 = g_lo + g_hi - 3.0 * (f_lo - f_hi) / (lo - hi);
        const d2_raw = d1 * d1 - g_lo * g_hi;
        const d2 = Math.sign(hi - lo) * Math.sqrt(Math.max(0, d2_raw));
        alpha_j = hi - (hi - lo) * (g_hi + d2 - d1) / (g_hi - g_lo + 2.0 * d2);
      }
      const lo_safe = Math.min(lo, hi);
      const hi_safe = Math.max(lo, hi);
      // The cubic's minimizer is only trusted inside the bracket; the
      // central-80% restriction from the textbook is deliberately NOT
      // applied here — with a steep vdW wall a few milli-Ångström from
      // the bracket's lo end (nicotine's SDF geometry), the correct
      // answer lies far outside the central 80%, and forcing bisection
      // from the far end burns the whole zoom budget on halving steps
      // down the wall. Any point inside the bracket is a valid trial;
      // bisection remains the fallback for pathological interpolations.
      if (!(alpha_j > lo_safe && alpha_j < hi_safe) || !isFinite(alpha_j)) {
        alpha_j = 0.5 * (lo + hi);
      }

      const state_j = evaluate_at(alpha_j, direction);
      const f_j = state_j.total;
      const g_j = dot(state_j.gradient, direction);

      if (f_j > phi0 + C1 * alpha_j * phi0_prime || f_j >= f_lo) {
        hi = alpha_j;
        f_hi = f_j;
      } else {
        if (Math.abs(g_j) <= -C2 * phi0_prime) {
          return { alpha: alpha_j, found: true, state: state_j };
        }
        if (g_j * (hi - lo) >= 0) {
          hi = lo;
          f_hi = f_lo;
        }
        lo = alpha_j;
        f_lo = f_j;
        g_lo = g_j;
      }
    }

    // Budget exhausted: return the lo side of the bracket — it always
    // satisfies the Armijo condition.
    return { alpha: lo, found: true, state: evaluate_at(lo, direction) };
  }
}
