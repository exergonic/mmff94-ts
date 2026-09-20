import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { parse_sdf } from '../src/sdf.js';
import { assign_atom_types } from '../src/mmff94/assign-atom-types.js';
import { calc_energy } from '../src/mmff94/index.js';
import { optimize_lbfgs } from '../src/optimize/l-bfgs.js';

// The α-imino carbocation of a five-membered ring (C₄H₆N⁺; the
// "1-pyrrolinium" drawn with the charge on the α-carbon) — found
// through the WebMO differential of 2026-09-20, where the submitted
// cation reached the engine with its charge dropped and its
// three-coordinate carbon typed as an alkyl carbon, so the conjugated
// C–N–C came back bent.
//
// MMFF94 defines no carbocation type. A three-coordinate carbon
// carrying +1 is an sp² center with an empty p orbital, so it takes
// the vinylic class (2) — the convention RDKit's MMFF94 typing
// follows. The charge must be visible to the typer (the SDF's M CHG
// block) or the atom falls back to type 1 and there is nothing in the
// potential to keep the α-carbon planar.
//
// Reference minimum: RDKit 2026.03, MMFF94, 10.6240 kcal/mol, flat;
// this library from the same fixture: 10.6244. The 761-suite has no
// carbocations and cannot arbitrate this class, so RDKit is the
// reference and the fixture is pinned in INTENTIONALLY_UNREFERENCED.
const SDF = readFileSync('tests/fixtures/sdf/pyrrolinium-carbocation.sdf', 'utf-8');
// Atom order (0-based): 0 C⁺ · 1 N · 2 C=N · 3 CH₂ · 4 CH₂ · 5–10 H.

function geom(a: { x: number; y: number; z: number }[], i: number, j: number) {
  return Math.hypot(a[i].x - a[j].x, a[i].y - a[j].y, a[i].z - a[j].z);
}
function ang(a: { x: number; y: number; z: number }[], i: number, j: number, k: number) {
  const ax = a[i].x - a[j].x, ay = a[i].y - a[j].y, az = a[i].z - a[j].z;
  const bx = a[k].x - a[j].x, by = a[k].y - a[j].y, bz = a[k].z - a[j].z;
  const na = Math.hypot(ax, ay, az), nb = Math.hypot(bx, by, bz);
  return (Math.acos(Math.max(-1, Math.min(1, (ax * bx + ay * by + az * bz) / (na * nb)))) * 180) / Math.PI;
}
function dih(a: { x: number; y: number; z: number }[], i: number, j: number, k: number, l: number) {
  const sub = (p: number, q: number) => [a[p].x - a[q].x, a[p].y - a[q].y, a[p].z - a[q].z];
  const b0 = sub(i, j), b1 = sub(k, j), b2 = sub(l, k);
  const n1 = Math.hypot(...b1);
  const bh = b1.map((v) => v / n1);
  const dot = (u: number[], v: number[]) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const v = b0.map((x, m) => x - dot(b0, bh) * bh[m]);
  const w = b2.map((x, m) => x - dot(b2, bh) * bh[m]);
  const cr = [v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2], v[0] * w[1] - v[1] * w[0]];
  return (Math.atan2(dot(cr, bh), dot(v, w)) * 180) / Math.PI;
}

describe('carbocation typing (three-coordinate C⁺ → vinylic class 2, 2026-09-20)', () => {
  it('types the charged α-carbon as 2 and the imine nitrogen as 9', () => {
    const t = assign_atom_types(parse_sdf(SDF));
    expect(t.atom_types).toEqual([2, 9, 3, 1, 1, 5, 5, 5, 5, 5, 5]);
  });

  it('reads the formal charge from the SDF M CHG block', () => {
    const mol = parse_sdf(SDF);
    expect(mol.atoms[0].formal_charge).toBe(1);
    for (let i = 1; i < mol.atoms.length; i++) {
      expect(mol.atoms[i].formal_charge ?? 0).toBe(0);
    }
  });

  it('without the charge the same carbon falls back to the alkyl type 1', () => {
    // The pre-fix behaviour, kept as a contrast: the rule is the charge,
    // not the connectivity — a neutral three-coordinate carbon is not a
    // carbocation and must not be typed as one.
    const mol = parse_sdf(SDF);
    for (const a of mol.atoms) a.formal_charge = 0;
    expect(assign_atom_types(mol).atom_types[0]).toBe(1);
  });

  it('minimises the conjugated C–N–C flat, at the RDKit minimum', () => {
    const mol = parse_sdf(SDF);
    const start = calc_energy(mol).total;
    expect(start).toBeCloseTo(21.1495, 1);           // the ETKDG start, unrelaxed
    const res = optimize_lbfgs(mol, { gradient_tolerance: 0.05 });
    expect(res.energy.total).toBeCloseTo(10.6244, 2); // RDKit 2026.03: 10.6240
    const a = res.molecule.atoms;
    // Both sp² centres stay trigonal planar …
    expect(ang(a, 4, 0, 1) + ang(a, 4, 0, 5) + ang(a, 1, 0, 5)).toBeCloseTo(360, 0);
    expect(ang(a, 1, 2, 3) + ang(a, 1, 2, 6) + ang(a, 3, 2, 6)).toBeCloseTo(360, 0);
    // … the C⁺–N=C torsion collapses to the conjugated plane …
    expect(Math.abs(dih(a, 4, 0, 1, 2))).toBeLessThan(2);
    // … and the C⁺–N bond contracts to the delocalised 1.366 Å
    // (RDKit 1.366; the alkyl-typed minimum keeps 1.474 Å).
    expect(geom(a, 0, 1)).toBeCloseTo(1.366, 2);
  });
});
