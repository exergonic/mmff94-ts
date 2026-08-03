/**
 * L-BFGS geometry optimization tests.
 *
 * Per the AGENTS.md Phase 6 spec: optimize each fixture to
 * max |gradient| < 0.05 kcal/mol/Å — from the SDF geometry AND from
 * a perturbed geometry — and check that the final energy is lower
 * (or equal, for fixtures already at a minimum) and convergence is
 * reported honestly.
 *
 * The perturbation matters: a minimizer that only works from a good
 * guess is not a minimizer. The SDF fixtures are the "given" starting
 * points; the perturbed runs shake every atom by ~0.2 Å.
 */

import { describe, it, expect } from 'vitest';
import { readFileSync, readdirSync } from 'fs';
import { join, parse } from 'path';
import { parse_sdf } from '../src/sdf';
import { assign_atom_types } from '../src/mmff94/atom-types';
import { compute_bci_charges } from '../src/mmff94/charges';
import { calc_energy } from '../src/mmff94/energy/total';
import { calc_gradient } from '../src/mmff94/gradient/total';
import { optimize_lbfgs, type EnergyGradientFn } from '../src/optimize/l-bfgs';
import { optimize_steepest_descent } from '../src/optimize/steepest-descent';
import type { TypedMolecule } from '../src/types';

const SDF_DIR = join(__dirname, 'fixtures', 'sdf');

const GRADIENT_TOL = 0.05; // kcal/mol/Å — the Phase 6 spec

// All fixtures are asserted at the Phase 6 spec (max |gradient| < 0.05).
// formamide needed a caveat until the amide-N typing landed (types
// 10/28): its wrong-typed surface had an artificial minimum. With the
// typing fixed, the SDF geometry IS the MMFF94 minimum and the fixture
// converges like the others.
const OPTIMIZER_CAVEATS: Record<string, { tolerance?: number; skip?: string }> = {};

// Deterministic pseudo-random perturbation (LCG) — the test must be
// reproducible across runs and platforms.
function perturb(molecule: TypedMolecule, seed: number, amplitude: number): TypedMolecule {
  let s = seed >>> 0;
  const rand = () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return (s / 0x100000000) * 2 - 1; // [-1, 1)
  };
  const clone: TypedMolecule = {
    ...molecule,
    atoms: molecule.atoms.map(a => ({ ...a })),
  };
  for (const atom of clone.atoms) {
    atom.x += amplitude * rand();
    atom.y += amplitude * rand();
    atom.z += amplitude * rand();
  }
  return clone;
}

function energy_gradient_fn(): EnergyGradientFn {
  return mol => ({ energy: calc_energy(mol), gradient: calc_gradient(mol) });
}

function max_gradient(gradient: number[][]): number {
  let m = 0;
  for (const g of gradient) for (const v of g) m = Math.max(m, Math.abs(v));
  return m;
}

describe('L-BFGS optimization', () => {
  for (const file of readdirSync(SDF_DIR).filter(f => f.endsWith('.sdf'))) {
    const name = parse(file).name;
    const raw = parse_sdf(readFileSync(join(SDF_DIR, file), 'utf-8'));
    raw.name = name;
    const typed = assign_atom_types(raw);
    // The charged copy flows into the optimizer (and into perturb()
    // below, whose spread keeps the charges — valid because they are
    // geometry-independent).
    const charged = compute_bci_charges(typed);

    const starting_gradient = max_gradient(calc_gradient(charged));
    const starting_energy = calc_energy(charged).total;

    it(`${name}: converges from the SDF geometry and from a perturbed one`, () => {
      const caveat = OPTIMIZER_CAVEATS[name];
      const tol = caveat?.tolerance ?? GRADIENT_TOL;
      if (caveat?.skip) {
        // Documented caveat (typing gap / pathological surface): the
        // optimizer still runs and must descend, but no convergence
        // claim is made. See OPTIMIZER_CAVEATS above.
        const from_sdf = optimize_lbfgs(charged, energy_gradient_fn(), { gradient_tolerance: tol, max_iterations: 500 });
        expect(from_sdf.energy.total).toBeLessThan(starting_energy + 1e-6);
        return;
      }

      // Run 1: from the SDF geometry — the "given" starting point.
      const from_sdf = optimize_lbfgs(charged, energy_gradient_fn(), {
        gradient_tolerance: tol,
      });
      expect(from_sdf.converged).toBe(true);
      expect(from_sdf.final_max_gradient).toBeLessThan(tol);
      // The energy must not go uphill (fixtures already at a minimum
      // — e.g. ethane — may only match the starting value).
      expect(from_sdf.energy.total).toBeLessThanOrEqual(starting_energy + 1e-6);
      // A molecule starting far from a minimum must actually descend
      // (the gate scales with the effective tolerance — nicotine's
      // 0.57 kcal/mol/Å start at tol = 1.0 is not "far").
      if (starting_gradient > 10 * tol) {
        expect(from_sdf.energy.total).toBeLessThan(starting_energy - 1e-3);
      }

      // Run 2: from a perturbed geometry — every atom shaken by ~0.2 Å.
      const perturbed = perturb(charged, name.length * 7919 + 17, 0.2);
      const perturbed_start = calc_energy(perturbed).total;
      const from_perturbed = optimize_lbfgs(perturbed, energy_gradient_fn(), {
        gradient_tolerance: tol,
      });
      expect(from_perturbed.converged).toBe(true);
      expect(from_perturbed.final_max_gradient).toBeLessThan(tol);
      // Must descend meaningfully from the shaken start.
      expect(from_perturbed.energy.total).toBeLessThan(perturbed_start - 1e-3);
      // And land in the same basin as the SDF run. The comparison is
      // deliberately loose: the tolerance lets a trajectory stop
      // anywhere inside a small ball around the minimum, so the two
      // finals can differ by ~0.001 kcal/mol even for identical basins
      // (pyridine: 4e-5; trimethylamine: 5e-4).
      expect(Math.abs(from_perturbed.energy.total - from_sdf.energy.total)).toBeLessThan(0.5);
    });
  }

  it('respects a tighter gradient tolerance', () => {
    const raw = parse_sdf(readFileSync(join(SDF_DIR, 'butane.sdf'), 'utf-8'));
    const typed = assign_atom_types(raw);
    const charged = compute_bci_charges(typed);
    const result = optimize_lbfgs(charged, energy_gradient_fn(), {
      gradient_tolerance: 0.005,
    });
    expect(result.converged).toBe(true);
    expect(result.final_max_gradient).toBeLessThan(0.005);
  });

  it('does not mutate the input molecule', () => {
    const raw = parse_sdf(readFileSync(join(SDF_DIR, 'propane.sdf'), 'utf-8'));
    const typed = assign_atom_types(raw);
    const charged = compute_bci_charges(typed);
    const before = typed.atoms.map(a => [a.x, a.y, a.z]);
    optimize_lbfgs(charged, energy_gradient_fn());
    expect(typed.atoms.map(a => [a.x, a.y, a.z])).toEqual(before);
  });

  it('handles an ill-conditioned quadratic (condition 10⁴)', () => {
    // E = x² + 10⁴·y² on a single fake atom. L-BFGS must converge in a
    // handful of iterations despite the 10⁴ curvature ratio. (Condition
    // 10⁶ is the documented boundary where the 2 Å trial-step cap makes
    // the small-curvature mode crawl — real MMFF94 surfaces sit at
    // condition ~10⁴ or better, e.g. the bond/angle stiffness ratio.)
    const mol: TypedMolecule = {
      atoms: [{ index: 0, element: 'C', x: 3, y: 4, z: 0 }],
      bonds: [],
      atom_types: [1],
      partial_charges: [0],
    };
    const cond = 1e4;
    const result = optimize_lbfgs(
      mol,
      m => {
        const x = m.atoms[0].x;
        const y = m.atoms[0].y;
        const energy = x * x + cond * y * y;
        const gradient = [[2 * x, 2 * cond * y, 0]];
        return { energy: { total: energy, bond_stretch: 0, angle_bend: 0, stretch_bend: 0, torsion: 0, van_der_waals: 0, electrostatic: 0, out_of_plane: 0 }, gradient };
      },
      { gradient_tolerance: 1e-8, max_iterations: 100 },
    );
    expect(result.converged).toBe(true);
    expect(result.final_max_gradient).toBeLessThan(1e-8);
    expect(Math.abs(result.molecule.atoms[0].x)).toBeLessThan(1e-6);
    expect(Math.abs(result.molecule.atoms[0].y)).toBeLessThan(1e-6);
  });
});

describe('Steepest descent fallback', () => {
  // The fallback converges the same fixtures at the spec, but only
  // linearly: the valley zig-zag costs 20-400+ iterations where
  // L-BFGS needs a handful. nicotine is the honest boundary — its vdW
  // canyon defeats 1000 SD iterations (|g|∞ stalls near 1.1; L-BFGS
  // escapes in 441), so it is asserted to descend but not to converge.
  const SD_CAVEATS: Record<string, { skip?: string }> = {
    nicotine: {
      skip:
        'vdW canyon defeats the linearly-convergent fallback in 1000 iterations — descent only',
    },
  };

  for (const file of readdirSync(SDF_DIR).filter(f => f.endsWith('.sdf'))) {
    const name = parse(file).name;
    const raw = parse_sdf(readFileSync(join(SDF_DIR, file), 'utf-8'));
    raw.name = name;
    const charged = compute_bci_charges(assign_atom_types(raw));
    const starting_energy = calc_energy(charged).total;

    it(`${name}: descends and converges at the spec from both starts`, () => {
      const caveat = SD_CAVEATS[name];

      // Run 1: from the SDF geometry.
      const from_sdf = optimize_steepest_descent(charged, energy_gradient_fn());
      // Armijo guarantees descent; at the very least the energy never
      // goes uphill (fixtures already at a minimum may not move).
      expect(from_sdf.energy.total).toBeLessThan(starting_energy + 1e-6);
      if (caveat) return; // descent asserted; no convergence claim
      expect(from_sdf.converged).toBe(true);
      expect(from_sdf.final_max_gradient).toBeLessThan(GRADIENT_TOL);

      // Run 2: from a perturbed geometry — every atom shaken by ~0.2 Å.
      const perturbed = perturb(charged, name.length * 7919 + 17, 0.2);
      const perturbed_start = calc_energy(perturbed).total;
      const from_perturbed = optimize_steepest_descent(perturbed, energy_gradient_fn());
      expect(from_perturbed.energy.total).toBeLessThan(perturbed_start - 1e-3);
      expect(from_perturbed.converged).toBe(true);
      expect(from_perturbed.final_max_gradient).toBeLessThan(GRADIENT_TOL);
    });
  }

  it('does not mutate the input molecule', () => {
    const raw = parse_sdf(readFileSync(join(SDF_DIR, 'propane.sdf'), 'utf-8'));
    const typed = assign_atom_types(raw);
    const charged = compute_bci_charges(typed);
    const before = typed.atoms.map(a => [a.x, a.y, a.z]);
    optimize_steepest_descent(charged, energy_gradient_fn());
    expect(typed.atoms.map(a => [a.x, a.y, a.z])).toEqual(before);
  });
});
