/**
 * Steepest descent optimizer — simple fallback minimizer.
 *
 * This is a robust but slow algorithm used when L-BFGS fails
 * (e.g., for highly strained starting geometries). It moves each
 * atom in the direction of the negative gradient (i.e., downhill).
 *
 * STATUS: NOT YET IMPLEMENTED — returns the input geometry with
 * converged: false (placeholder). The description below is the
 * design for the real implementation (AGENTS.md Phase 6).
 *
 *   x_{k+1} = x_k − α · ∇E(x_k)
 *
 * The step size α is determined by an Armijo backtracking line search:
 * start with α = 1.0 and halve it until the Armijo condition is satisfied:
 *
 *   E(x_k − α·∇E) < E(x_k) − c₁ · α · |∇E|²
 *
 * where c₁ = 1e-4 (standard value).
 *
 * This algorithm is slow (it zig-zags down valleys), but it is
 * guaranteed to converge to a local minimum if the line search
 * satisfies the Wolfe conditions.
 */

import type { Molecule, EnergyComponents, OptimizationResult } from '../types';

export interface SteepestDescentOptions {
  max_iterations?: number;
  gradient_tolerance?: number;
  initial_step_size?: number;
}

/**
 * Optimize geometry using steepest descent with Armijo line search.
 *
 * Primarily a fallback for when L-BFGS fails to converge.
 */
export function optimize_steepest_descent(
  molecule: Molecule,
  calc_energy_gradient: (mol: Molecule) => { energy: EnergyComponents; gradient: number[][] },
  options?: SteepestDescentOptions
): OptimizationResult {
  void options;  // part of the Phase 6 API; unused until the algorithm lands
  // TODO: implement steepest descent with Armijo line search.
  //
  // For each iteration:
  //   1. Compute energy + gradient.
  //   2. Check convergence: max |gradient| < tolerance.
  //   3. Line search: backtrack from a starting step size.
  //   4. Update positions: x = x - step * gradient.
  //   5. Repeat.
  //
  // Return the result.
  const energy = calc_energy_gradient(molecule).energy;
  return {
    molecule,
    energy,
    iterations: 0,
    converged: false,
    final_max_gradient: Infinity,
  };
}
