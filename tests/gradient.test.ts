/**
 * Finite-difference cross-checks for the analytical gradients.
 *
 * For every atom of every fixture, every Cartesian coordinate is
 * perturbed by ±δ = 10⁻⁶ Å and the energy term is evaluated at both
 * points; the central difference (E(x+δ) − E(x−δ))/(2δ) must match
 * the analytical gradient to a relative error < 10⁻⁵.
 *
 * The tolerance floor of 10⁻⁵ (kcal/mol/Å) covers gradient components
 * that are genuinely zero by symmetry — a relative error is undefined
 * there, but an absolute error of 10⁻⁵ is still well below the finite-
 * difference noise floor (~10⁻⁷ for energies in the 10s of kcal/mol).
 *
 * This test is the safety net for the one rule that matters in the
 * gradient phase: the gradient of a term must visit EXACTLY the same
 * interactions as the energy term (same parameter lookups, same pair
 * lists, same 1-4 scaling) — any divergence shows up here as a
 * mismatch on the affected atoms.
 */

import { describe, it, expect } from 'vitest';
import { readFileSync, readdirSync } from 'fs';
import { join, parse } from 'path';
import { parse_sdf } from '../src/sdf';
import { parse_mmd } from '../src/utils/mmd-parser';
import { assign_atom_types } from '../src/mmff94/assign-atom-types';
import { assign_bci_charges } from '../src/mmff94/charges';
import type { TypedMolecule } from '../src/types';

import { calc_bond_stretch_energy } from '../src/mmff94/energy/bond-stretch';
import { calc_angle_bend_energy } from '../src/mmff94/energy/angle-bend';
import { calc_stretch_bend_energy } from '../src/mmff94/energy/stretch-bend';
import { calc_torsion_energy } from '../src/mmff94/energy/torsion';
import { calc_vdw_energy } from '../src/mmff94/energy/van-der-waals';
import { calc_electrostatic_energy } from '../src/mmff94/energy/electrostatic';
import { calc_oop_energy } from '../src/mmff94/energy/out-of-plane';
import { calc_bond_stretch_gradient } from '../src/mmff94/gradient/bond-stretch';
import { calc_angle_bend_gradient } from '../src/mmff94/gradient/angle-bend';
import { calc_stretch_bend_gradient } from '../src/mmff94/gradient/stretch-bend';
import { calc_torsion_gradient } from '../src/mmff94/gradient/torsion';
import { calc_vdw_gradient } from '../src/mmff94/gradient/van-der-waals';
import { calc_electrostatic_gradient } from '../src/mmff94/gradient/electrostatic';
import { calc_oop_gradient } from '../src/mmff94/gradient/out-of-plane';
import { calc_gradient } from '../src/mmff94/gradient/total';
import { calc_energy } from '../src/mmff94/energy/total';

const SDF_DIR = join(__dirname, 'fixtures', 'sdf');
const SUITE_DIR = join(__dirname, 'fixtures', 'validation-suite');

// All SDF fixtures, plus the three regression-pinned suite molecules
// (FUVDOP's 3-ring closure, FILNOD's fused 5-ring, JIYJAC's class-2
// strbnd) — the molecules where a parameter-resolution subtlety was
// once wrong are the ones most likely to hide a gradient subtlety.
const SUITE_MOLECULES = ['FUVDOP', 'FILNOD', 'JIYJAC'];

const DELTA = 1e-6; // Å — small enough for 10⁻⁵ relative accuracy, big enough to beat roundoff
const REL_TOL = 1e-5;
const ABS_FLOOR = 1e-5; // kcal/mol/Å

function load_fixtures(): TypedMolecule[] {
  const molecules: TypedMolecule[] = [];
  for (const file of readdirSync(SDF_DIR).filter(f => f.endsWith('.sdf'))) {
    const molecule = parse_sdf(readFileSync(join(SDF_DIR, file), 'utf-8'));
    molecule.name = parse(file).name;
    molecules.push(assign_atom_types(molecule));
  }
  const suite_text = readFileSync(join(SUITE_DIR, 'MMFF94.mmd'), 'utf-8');
  const suite_molecules = parse_mmd(suite_text);
  for (const code of SUITE_MOLECULES) {
    const mol = suite_molecules.find(m => m.name === code);
    if (mol) molecules.push(assign_atom_types(mol));
  }
  return molecules;
}

/**
 * Central finite-difference check of one term over all fixtures.
 * Returns the worst relative error found (or 0 if every component
 * passed), printing per-fixture worst errors as it goes.
 */
function check_term(
  name: string,
  energy_fn: (m: TypedMolecule) => number,
  gradient_fn: (m: TypedMolecule) => number[][],
): number {
  const molecules = load_fixtures();
  let worst_overall = 0;
  console.log(`\n=== ${name} gradient vs finite differences (δ = ${DELTA} Å) ===`);
  for (const molecule of molecules) {
    // Charges are geometry-independent; compute once (the returned
    // charged molecule shares the atoms, so the FD perturbations
    // below move the same geometry) so the electrostatic term (and
    // its gradient) stay consistent.
    const charged = assign_bci_charges(molecule);

    const analytic = gradient_fn(charged);
    let worst = 0;
    let worst_atom = -1;

    for (let a = 0; a < molecule.atoms.length; a++) {
      for (let axis = 0; axis < 3; axis++) {
        const original = molecule.atoms[a][axis === 0 ? 'x' : axis === 1 ? 'y' : 'z'];
        molecule.atoms[a][axis === 0 ? 'x' : axis === 1 ? 'y' : 'z'] = original + DELTA;
        const e_plus = energy_fn(charged);
        molecule.atoms[a][axis === 0 ? 'x' : axis === 1 ? 'y' : 'z'] = original - DELTA;
        const e_minus = energy_fn(charged);
        molecule.atoms[a][axis === 0 ? 'x' : axis === 1 ? 'y' : 'z'] = original;

        const fd = (e_plus - e_minus) / (2 * DELTA);
        const exact = analytic[a][axis];
        const scale = Math.max(Math.abs(fd), Math.abs(exact), 1.0);
        const rel = Math.abs(fd - exact) / scale;
        if (rel > worst) {
          worst = rel;
          worst_atom = a;
        }
        if (rel > REL_TOL || Math.abs(fd - exact) > ABS_FLOOR) {
          throw new Error(
            `${name}: ${molecule.name} atom ${a} axis ${axis}: ` +
            `finite difference ${fd.toExponential(4)} vs analytical ${exact.toExponential(4)}`,
          );
        }
      }
    }
    worst_overall = Math.max(worst_overall, worst);
    const flag = worst > 1e-7 ? '  ← largest' : '';
    console.log(`  ${(molecule.name ?? '?').padEnd(14)} worst rel err ${worst.toExponential(2)} (atom ${worst_atom})${flag}`);
  }
  console.log(`  worst over all fixtures: ${worst_overall.toExponential(2)}`);
  return worst_overall;
}

describe('gradient finite-difference checks', () => {
  it('bond stretch', () => {
    const worst = check_term('bond_stretch', calc_bond_stretch_energy, calc_bond_stretch_gradient);
    expect(worst).toBeLessThan(1e-5);
  });

  it('angle bend', () => {
    const worst = check_term('angle_bend', calc_angle_bend_energy, calc_angle_bend_gradient);
    expect(worst).toBeLessThan(1e-5);
  });

  it('stretch-bend', () => {
    const worst = check_term('stretch_bend', calc_stretch_bend_energy, calc_stretch_bend_gradient);
    expect(worst).toBeLessThan(1e-5);
  });

  it('torsion', () => {
    const worst = check_term('torsion', calc_torsion_energy, calc_torsion_gradient);
    expect(worst).toBeLessThan(1e-5);
  });

  it('van der Waals', () => {
    const worst = check_term('van_der_waals', calc_vdw_energy, calc_vdw_gradient);
    expect(worst).toBeLessThan(1e-5);
  });

  it('electrostatic (with 1-4 ×0.75)', () => {
    const worst = check_term('electrostatic', calc_electrostatic_energy, calc_electrostatic_gradient);
    expect(worst).toBeLessThan(1e-5);
  });

  it('out-of-plane', () => {
    const worst = check_term('out_of_plane', calc_oop_energy, calc_oop_gradient);
    expect(worst).toBeLessThan(1e-5);
  });

  it('total gradient = sum of the seven term gradients', () => {
    // The wiring in gradient/total.ts is trivial, but a missing term
    // would silently vanish from the public API — pin the identity.
    const molecules = load_fixtures();
    for (const molecule of molecules) {
      const charged = assign_bci_charges(molecule);
      const total = calc_gradient(charged);
      const sum = molecule.atoms.map(() => [0, 0, 0]);
      for (const term of [
        calc_bond_stretch_gradient(charged),
        calc_angle_bend_gradient(charged),
        calc_stretch_bend_gradient(charged),
        calc_torsion_gradient(charged),
        calc_vdw_gradient(charged),
        calc_electrostatic_gradient(charged),
        calc_oop_gradient(charged),
      ]) {
        for (let a = 0; a < molecule.atoms.length; a++) {
          for (let axis = 0; axis < 3; axis++) sum[a][axis] += term[a][axis];
        }
      }
      for (let a = 0; a < molecule.atoms.length; a++) {
        for (let axis = 0; axis < 3; axis++) {
          expect(Math.abs(total[a][axis] - sum[a][axis])).toBeLessThan(1e-9);
        }
      }
    }
  });

  it('total gradient matches finite differences of the total energy', () => {
    const worst = check_term('total', m => calc_energy(m).total, calc_gradient);
    expect(worst).toBeLessThan(1e-5);
  });

  it('angle-bend gradient at an exactly-linear angle is the FD limit, not zero', () => {
    // A sketcher/embedder can hand the force field an exactly-linear
    // water (H-O-H = 180°). The 1/sin θ analytical form is a 0/0 limit
    // there; the FD fallback must give the true force — which bends
    // the molecule, not stalls it.
    const water = assign_atom_types({
      name: 'linear-water',
      atoms: [
        { index: 0, element: 'O', x: 0, y: 0, z: 0 },
        { index: 1, element: 'H', x: 1, y: 0, z: 0 },
        { index: 2, element: 'H', x: -1, y: 0, z: 0 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 1 },
        { atom1: 0, atom2: 2, bond_order: 1 },
      ],
    });
    const charged = assign_bci_charges(water);
    const analytic = calc_angle_bend_gradient(charged);

    // The transverse force is nonzero — the angle wants to bend.
    expect(Math.abs(analytic[1][1])).toBeGreaterThan(1e-3);

    // The energy is CUSPED in the transverse coordinate (like |x|), so
    // a central energy FD averages the two sides to zero — compare
    // against the forward difference instead, which matches the
    // one-sided angle derivative the fallback computes.
    const e0 = calc_angle_bend_energy(charged);
    for (let a = 0; a < 3; a++) {
      for (let axis = 0; axis < 3; axis++) {
        const key = axis === 0 ? 'x' : axis === 1 ? 'y' : 'z';
        const original = water.atoms[a][key];
        water.atoms[a][key] = original + DELTA;
        const e_plus = calc_angle_bend_energy(charged);
        water.atoms[a][key] = original;
        const fd = (e_plus - e0) / DELTA;
        // The suite's relative convention (a forward difference carries
        // an O(δ) relative truncation — the same scale check_term uses).
        const scale = Math.max(Math.abs(fd), Math.abs(analytic[a][axis]), 1.0);
        expect(Math.abs(fd - analytic[a][axis]) / scale).toBeLessThan(1e-5);
      }
    }
  });

  it('oop gradient at χ = 90° is the FD limit, not zero', () => {
    // A tri-coordinate center (carbonyl C) with a substituent exactly
    // perpendicular to the reference plane: the 1/cos χ form is a 0/0
    // limit; the FD fallback must push the substituent back toward
    // the plane. (An amine N would be a vacuous test — its k_oop is
    // zero by design.)
    const carbonyl = assign_atom_types({
      name: 'perpendicular-carbonyl',
      atoms: [
        { index: 0, element: 'C', x: 0, y: 0, z: 0 },
        { index: 1, element: 'O', x: 1, y: 0, z: 0 },
        { index: 2, element: 'H', x: 0, y: 1, z: 0 },
        { index: 3, element: 'H', x: 0, y: 0, z: 1 },
      ],
      bonds: [
        { atom1: 0, atom2: 1, bond_order: 2 },
        { atom1: 0, atom2: 2, bond_order: 1 },
        { atom1: 0, atom2: 3, bond_order: 1 },
      ],
    });
    const charged = assign_bci_charges(carbonyl);
    const analytic = calc_oop_gradient(charged);

    // The force is nonzero — the perpendicular H is pushed sideways
    // back toward the plane (the component along the plane normal is
    // genuinely zero: moving along the normal keeps χ = 90°).
    expect(Math.abs(analytic[3][0])).toBeGreaterThan(1e-4);

    // Same cusp convention as the angle test: forward differences.
    const e0 = calc_oop_energy(charged);
    for (let a = 0; a < 4; a++) {
      for (let axis = 0; axis < 3; axis++) {
        const key = axis === 0 ? 'x' : axis === 1 ? 'y' : 'z';
        const original = carbonyl.atoms[a][key];
        carbonyl.atoms[a][key] = original + DELTA;
        const e_plus = calc_oop_energy(charged);
        carbonyl.atoms[a][key] = original;
        const fd = (e_plus - e0) / DELTA;
        const scale = Math.max(Math.abs(fd), Math.abs(analytic[a][axis]), 1.0);
        expect(Math.abs(fd - analytic[a][axis]) / scale).toBeLessThan(1e-5);
      }
    }
  });
});
