import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { parse_sdf } from '../src/sdf.js';
import { assign_atom_types } from '../src/mmff94/assign-atom-types.js';
import { angle_parameters } from '../src/mmff94/parameters/parameter-classes.js';
import { torsion_terms } from '../src/mmff94/energy/torsion.js';
import { make_class_context } from '../src/mmff94/parameters/parameter-classes.js';
import { optimize_lbfgs } from '../src/optimize/l-bfgs.js';
import { calc_energy, calc_gradient } from '../src/mmff94/index.js';
import type { TypedMolecule } from '../src/types.js';

// Methylenetriphenylphosphorane (Ph₃P=CH₂, PubChem CID 137960) — the
// 2026-08-12 dogfooding report (the molecule rendered H-less and
// asymmetric on the local path). The H-less bug was the VALENCE
// embedder's multi-ring ring-walk crash; the asymmetry was THIS
// typing gap: our P=C branch typed the 4-σ P as 75 (the crd-2 ylide
// type), so the empirical angle protocol emitted θ₀ = 94.9° at P and
// the potential's minimum became an asymmetric tripod (measured
// three-way vs OpenBabel and Tinker, 2026-08-12: with P=25 all three
// engines converge to a symmetric C–P–C ≈ 104°/114° tripod; with 75
// both our L-BFGS AND Tinker find P–C(aryl) spreads of 0.07–0.13 Å
// and C–P–C 93.7–125.3°). The rule is now crd-aware: P doubly
// bonded to C/N with 4 σ types 25 (phosphonium/phosphorane), with
// 2 σ types 75 (H–P=N–OCH₃ — the 2026-08-10 phosphine-imide fix
// stays). The suite cannot arbitrate: zero P=C/P=N in the 761.
const SDF = readFileSync('tests/fixtures/sdf/methylenetriphenylphosphorane.sdf', 'utf-8');

function typed(): TypedMolecule {
  return assign_atom_types(parse_sdf(SDF));
}

describe('methylenetriphenylphosphorane typing (P=C → 25, crd-aware, 2026-08-12)', () => {
  it('types the 4-σ ylide P as 25 and the CH₂ as 3, matching OpenBabel\'s SDF typing', () => {
    // OB's canonical SDF typing (measured): P=25, ylide C=3, aryl C=37,
    // H=5 — identical except the P, which we now match. The earlier
    // "OB types the ylide C as 60" reading was a SMILES-parse artifact
    // (the SDF with explicit H's types it 3, same as us).
    const t = typed();
    expect(t.atom_types[1]).toBe(25);  // P — the tetracoordinate PO4 family
    expect(t.atom_types[0]).toBe(3);   // ylide CH₂ — the generic sp² C (C=O covers C=P)
    expect(t.atom_types[2]).toBe(37);  // aryl C
    expect(t.atom_types[20]).toBe(5);  // H
    // The old branch gave 75 here — the crd-2 ylide type on a
    // 4-coordinate center (the asymmetric-tripod bug).
  });

  it('resolves the C–P–C angles to Tinker\'s refs (99.158°/107.124°), not the crd-2 94.9°', () => {
    // The mechanism of the render fix: type 25 is crd 4, so the
    // empirical θ₀ protocol emits the tetracoordinate-P refs —
    // Tinker's analyze prints EXACTLY 99.158 (KB 1.072) for
    // C(ylide)-P-C(aryl) and 107.124 (KB 0.947) for C(aryl)-P-C(aryl)
    // on its 25/3 leg; our resolutions match to 6 decimals. The old
    // 75 (crd 2) emitted 94.9° — forcing four substituents onto a
    // 2-coordinate geometry (the asymmetric minimum).
    const t = typed();
    const adj: number[][] = Array.from({ length: t.atoms.length }, () => []);
    for (const b of t.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
    const ctx = make_class_context(t, adj);
    const aryl = angle_parameters(ctx, 0, 1, 2);   // C(ylide)-P-C(aryl)
    const arylAryl = angle_parameters(ctx, 2, 1, 8); // C(aryl)-P-C(aryl)
    expect(aryl.theta0).toBeCloseTo(99.158, 3);
    expect(arylAryl.theta0).toBeCloseTo(107.124, 3);
  });

  it('zeroes the P=C torsion — the par row 0-0-3-25-0, matching Tinker and OpenBabel', () => {
    // With P=25 the dihedral types (5,3,25,37) hit the par file's
    // zero row 0-0-3-25-0 (the same zero-row family as the
    // phosphine-imide case's 0-0-9-25-0): the 4-σ ylide's P=C
    // rotation is table-zero in all three engines (Tinker's and
    // OpenBabel's start-geometry torsion = 0.0055 — the ring rows
    // only). The C94 row 0-0-3-75-0 (V2 = 19.0, OB's mmfftor.par
    // "C94 0:*-3-75-*") now applies only to central (3,75) — the
    // 2-σ ylide P.
    const t = typed();
    const adj: number[][] = Array.from({ length: t.atoms.length }, () => []);
    for (const b of t.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
    const ctx = make_class_context(t, adj);
    const tr = torsion_terms(ctx, t, 20, 0, 1, 2); // H-C=P-C(aryl)
    expect(tr?.v1).toBe(0);
    expect(tr?.v2).toBe(0);
    expect(tr?.v3).toBe(0);
  });

  it('minimizes to a symmetric tripod at P (the app-visible regression)', () => {
    // The user-visible guarantee: from the PubChem/Avogadro geometry
    // (with the same symmetry-breaking kick the Valence embedder
    // applies — the start is an exact-symmetry point where L-BFGS
    // crawls), the MMFF94 minimum is a symmetric P tripod: P–C(aryl)
    // spread ≤ 0.02 Å, C–P–C within 12°, P=C ≈ 1.77 (measured
    // 0.0007/10.7°/1.7716 — the P=C matches Tinker's minimum
    // 1.7716 exactly). With the old 75 typing the minimum was
    // asymmetric (spread 0.067 Å, angles spread 21°).
    const t = typed();
    for (let i = 0; i < t.atoms.length; i++) {
      const seed = (i * 2654435761) >>> 0;
      t.atoms[i].x += 0.1 * ((seed >>> 13) % 7 - 3) / 3;
      t.atoms[i].y += 0.1 * ((seed >>> 7) % 7 - 3) / 3;
      t.atoms[i].z += 0.1 * ((seed >>> 3) % 7 - 3) / 3;
    }
    const res = optimize_lbfgs(
      t,
      (m) => ({ energy: calc_energy(m), gradient: calc_gradient(m) }),
      { gradient_tolerance: 0.05, max_iterations: 250 },
    );
    const a = res.molecule.atoms;
    const d = (i: number, j: number) => Math.hypot(
      a[i].x - a[j].x, a[i].y - a[j].y, a[i].z - a[j].z,
    );
    const ang = (i: number, j: number, k: number) => {
      const ax = a[i].x - a[j].x, ay = a[i].y - a[j].y, az = a[i].z - a[j].z;
      const bx = a[k].x - a[j].x, by = a[k].y - a[j].y, bz = a[k].z - a[j].z;
      const na = Math.hypot(ax, ay, az), nb = Math.hypot(bx, by, bz);
      return (Math.acos(Math.max(-1, Math.min(1, (ax * bx + ay * by + az * bz) / (na * nb)))) * 180) / Math.PI;
    };
    const pcAryl = [2, 8, 14].map((i) => d(1, i));
    const spread = Math.max(...pcAryl) - Math.min(...pcAryl);
    expect(spread).toBeLessThan(0.02);
    const cpc = [ang(0, 1, 2), ang(0, 1, 8), ang(0, 1, 14), ang(2, 1, 8), ang(2, 1, 14), ang(8, 1, 14)];
    expect(Math.max(...cpc) - Math.min(...cpc)).toBeLessThan(12);
    expect(d(1, 0)).toBeGreaterThan(1.7);
    expect(d(1, 0)).toBeLessThan(1.8);
  });
});
