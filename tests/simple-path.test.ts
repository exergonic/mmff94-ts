/**
 * The simple path: a bare Molecule straight from parse_sdf() is all
 * the top-level functions need.
 *
 * The library works two ways, and neither penalizes the other:
 *   - simple: calc_energy(mol) / calc_gradient(mol) /
 *     optimize_lbfgs(mol) — typing and charges happen on demand, and
 *     the results are identical to the rich path's;
 *   - rich: assign_atom_types → assign_bci_charges → calc_energy /
 *     calc_gradient / optimize with a custom oracle.
 *
 * These tests pin the contract: the simple path returns the SAME full
 * per-term breakdown, the optimizer results carry the typed/charged
 * molecule for follow-up work at the minimum, and the input is never
 * mutated.
 */

import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { parse_sdf } from '../src/sdf';
import { assign_atom_types } from '../src/mmff94/atom-types';
import { assign_bci_charges } from '../src/mmff94/charges';
import { calc_energy, calc_gradient } from '../src/mmff94';
import { optimize_lbfgs } from '../src/optimize/l-bfgs';
import { optimize_steepest_descent } from '../src/optimize/steepest-descent';
import type { TypedMolecule } from '../src/types';

const SDF_DIR = join(__dirname, 'fixtures', 'sdf');

describe('simple path — a bare Molecule is enough', () => {
  it('calc_energy on a bare Molecule equals the rich path, per term', () => {
    // nicotine exercises the charged typing + formal-charge model.
    const raw = parse_sdf(readFileSync(join(SDF_DIR, 'nicotine.sdf'), 'utf-8'));

    // Rich path: explicit typing + charges.
    const rich = calc_energy(assign_bci_charges(assign_atom_types(raw)));

    // Simple path: one call on the bare molecule.
    const simple = calc_energy(raw);

    // Identical full breakdown — the simple path is never a subset.
    expect(simple).toEqual(rich);
    expect(simple.total).toBe(rich.total);
    // The input is untouched by the on-demand preparation.
    expect((raw as TypedMolecule).atom_types).toBeUndefined();
    expect((raw as TypedMolecule).partial_charges).toBeUndefined();
  });

  it('calc_gradient on a bare Molecule equals the rich path', () => {
    const raw = parse_sdf(readFileSync(join(SDF_DIR, 'butane.sdf'), 'utf-8'));
    const rich = calc_gradient(assign_bci_charges(assign_atom_types(raw)));
    expect(calc_gradient(raw)).toEqual(rich);
  });

  it('optimize_lbfgs works with a bare Molecule and no callback', () => {
    const raw = parse_sdf(readFileSync(join(SDF_DIR, 'butane.sdf'), 'utf-8'));
    const result = optimize_lbfgs(raw);

    expect(result.converged).toBe(true);
    expect(result.final_max_gradient).toBeLessThan(0.05);
    // The result exposes everything: the full per-term breakdown …
    expect(result.energy.bond_stretch).toBeTypeOf('number');
    // … and the typed/charged molecule at the minimum, so the rich
    // path (per-term energies, charges, gradients) is reachable from
    // the simple path without re-doing anything.
    expect(result.molecule.atom_types).toBeDefined();
    expect(result.molecule.partial_charges).toBeDefined();
    expect(calc_energy(result.molecule).total).toBe(result.energy.total);
    // The input is untouched (the optimizer works on a prepared copy).
    expect((raw as TypedMolecule).atom_types).toBeUndefined();
  });

  it('optimize_steepest_descent works with a bare Molecule and no callback', () => {
    const raw = parse_sdf(readFileSync(join(SDF_DIR, 'ethane.sdf'), 'utf-8'));
    const result = optimize_steepest_descent(raw);
    expect(result.converged).toBe(true);
    expect(result.molecule.atom_types).toBeDefined();
    expect(result.molecule.partial_charges).toBeDefined();
  });
});
