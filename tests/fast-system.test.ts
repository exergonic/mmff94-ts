/**
 * Differential test: the compiled fast path vs the readable terms.
 *
 * The optimizer runs on the fast system (see src/optimize/fast-system.ts)
 * — one-time parameter resolution into typed arrays, zero-allocation
 * evaluation. This test is the truth-discipline guard for that kernel:
 *
 *   - ENERGY, bit-for-bit: bond/angle/stretch-bend/OOP/electrostatic
 *     (the fast kernel hoists their unit factors with the readable
 *     multiplication order — identical doubles, computed once).
 *   - ENERGY, ≤1e-9: torsion and van der Waals. Those two carry
 *     deliberate arithmetic rewrites (double/triple-angle identities
 *     replace cos(2τ)/cos(3τ); r⁷ is a multiplication chain instead of
 *     Math.pow(r,7)) — mathematically equal, ULP-rounded differently.
 *   - GRADIENT: fast vs readable per term, |Δ| ≤ 1e-8 absolute. The
 *     angle/dihedral/OOP accumulators are scalar expansions of the
 *     readable axis-loop algebra; reassociation moves last-ULP rounding
 *     only (and the readable side itself mixes Math.hypot and sqrt
 *     paths the fused kernel resolves to one).
 *
 * Both are checked at the fixtures' own geometries AND at deterministic
 * perturbations (0.05/0.2 Å per coordinate), so the comparison covers
 * off-equilibrium regions of every term — the cubic bond term, the
 * angle cubic near large deviations, torsion at arbitrary phases, and
 * the FD-fallback cusps stay out of scope by construction (never
 * reached away from exactly 0°/180°).
 */

import { describe, it, expect } from 'vitest';
import { readFileSync, readdirSync } from 'fs';
import { join } from 'path';
import { parse_sdf } from '../src/sdf';
import { parse_mmd } from '../src/utils/mmd-parser';
import { assign_atom_types } from '../src/mmff94/assign-atom-types';
import { assign_bci_charges } from '../src/mmff94/charges';
import { calc_energy } from '../src/mmff94/energy/total';
import { calc_bond_stretch_gradient } from '../src/mmff94/gradient/bond-stretch';
import { calc_angle_bend_gradient } from '../src/mmff94/gradient/angle-bend';
import { calc_stretch_bend_gradient } from '../src/mmff94/gradient/stretch-bend';
import { calc_torsion_gradient } from '../src/mmff94/gradient/torsion';
import { calc_vdw_gradient } from '../src/mmff94/gradient/van-der-waals';
import { calc_electrostatic_gradient } from '../src/mmff94/gradient/electrostatic';
import { calc_oop_gradient } from '../src/mmff94/gradient/out-of-plane';
import { create_fast_system } from '../src/optimize/fast-system';
import type { TypedMolecule } from '../src/types';

const SDF_DIR = join(__dirname, 'fixtures', 'sdf');
const SUITE_DIR = join(__dirname, 'fixtures', 'validation-suite');

// The readable per-term gradients, in the same order as the fast
// kernel's term_mask indexing (bond, angle, strbnd, torsion, oop, vdw,
// elec).
const READABLE_GRADIENTS = [
  calc_bond_stretch_gradient,
  calc_angle_bend_gradient,
  calc_stretch_bend_gradient,
  calc_torsion_gradient,
  calc_oop_gradient,
  calc_vdw_gradient,
  calc_electrostatic_gradient,
];

function flat(m: TypedMolecule): Float64Array {
  const x = new Float64Array(3 * m.atoms.length);
  for (let a = 0; a < m.atoms.length; a++) {
    x[3*a] = m.atoms[a].x; x[3*a+1] = m.atoms[a].y; x[3*a+2] = m.atoms[a].z;
  }
  return x;
}

/** Deterministic pseudo-random perturbation of every coordinate. */
function perturb(m: TypedMolecule, amplitude: number): TypedMolecule {
  let seed = m.atoms.length * 7919 + Math.round(amplitude * 1e6);
  const atoms = m.atoms.map(a => ({ ...a }));
  for (const a of atoms) {
    for (const key of ['x', 'y', 'z'] as const) {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      a[key] = (a[key] as number) + ((seed / 0x7fffffff) - 0.5) * 2 * amplitude;
    }
  }
  return { ...m, atoms };
}

function check(name: string, prepared: TypedMolecule): void {
  const sys = create_fast_system(prepared);
  const coords = flat(prepared);
  const grad = new Float64Array(coords.length);

  sys.evaluate(coords, grad, new Uint8Array(7).fill(1));

  // 1. Total energy: bit-for-bit where the kernel kept the readable
//    multiplication order, ≤1e-9 for the two ULP-rewritten terms
//    (torsion trig identities, vdW pow chains).
  const ref = calc_energy(prepared);
  const kit: [string, number][] = [
    ['bond_stretch', ref.bond_stretch],
    ['angle_bend', ref.angle_bend],
    ['stretch_bend', ref.stretch_bend],
    ['out_of_plane', ref.out_of_plane],
    ['electrostatic', ref.electrostatic],
  ];
  for (const [key, rv] of kit) {
    expect(sys.components[key as keyof typeof sys.components], `${name}: ${key} bitwise`).toBe(rv);
  }
  for (const [key, rv] of [['torsion', ref.torsion], ['van_der_waals', ref.van_der_waals]] as [string, number][]) {
    const d = Math.abs(sys.components[key as keyof typeof sys.components] - rv);
    expect(d, `${name}: ${key} ≤1e-9`).toBeLessThan(1e-9);
  }
  expect(Math.abs(sys.total - ref.total), `${name}: total`).toBeLessThan(1e-8);

  // 2. Per-term gradients within the reassociation tolerance.
  const one = new Uint8Array(7);
  for (let t = 0; t < 7; t++) {
    one.fill(0); one[t] = 1;
    sys.evaluate(coords, grad, one);
    const refG = READABLE_GRADIENTS[t](prepared);
    for (let a = 0; a < prepared.atoms.length; a++) {
      for (let ax = 0; ax < 3; ax++) {
        const d = Math.abs(grad[3*a + ax] - refG[a][ax]);
        expect(d, `${name}: term ${t} grad atom ${a} axis ${ax}`).toBeLessThan(1e-8);
      }
    }
  }
}

describe('fast system vs readable terms', () => {
  const sdfFiles = readdirSync(SDF_DIR).filter(f => f.endsWith('.sdf')).sort();
  const suite = parse_mmd(readFileSync(join(SUITE_DIR, 'MMFF94.mmd'), 'utf-8'));

  for (const f of sdfFiles) {
    const mol = parse_sdf(readFileSync(join(SDF_DIR, f), 'utf-8'));
    mol.name = f.replace('.sdf', '');
    // trpcage handled separately (larger, exercises everything incl. the
    // zwitterion charge flow — and it must not blow the runtime here
    // more than once at one amplitude).
    if (f === 'trpcage.sdf') continue;
    for (const amplitude of [0, 0.05, 0.2]) {
      it(`${mol.name} @ ±${amplitude} Å`, () => {
        check(mol.name!, perturb(assign_bci_charges(assign_atom_types(mol)), amplitude));
      });
    }
  }

  it('trpcage @ ±0.05 Å (zwitterion, 304 atoms)', () => {
    const mol = parse_sdf(readFileSync(join(SDF_DIR, 'trpcage.sdf'), 'utf-8'));
    mol.name = 'trpcage';
    check('trpcage', perturb(assign_bci_charges(assign_atom_types(mol)), 0.05));
  });

  // Charged chemistry (CUVJOS) and the gradient-regression suite set
  // (FUVDOP/FILNOD/JIYJAC) — the molecules the FD sweep pins.
  for (const code of ['FUVDOP', 'FILNOD', 'JIYJAC', 'CUVJOS']) {
    it(`suite ${code} @ ±0.1 Å`, () => {
      const m = suite.find(x => x.name === code)!;
      check(code, perturb(assign_bci_charges(assign_atom_types(m)), 0.1));
    });
  }
});