/**
 * The compiled fast path for MMFF94 energy + gradient.
 *
 * PHILOSOPHY. The readable term functions under src/mmff94/{energy,gradient}/
 * are the ground truth: chemistry vocabulary, why-comments, validated term by
 * term against BatchMin/TINKER/FD. This module is their machine-code twin for
 * the optimizer hot loop, written under a deliberately relaxed ethos: the domain
 * here is linear algebra and floating-point behavior, not chemistry, so the
 * tricks of the trade (flat typed arrays, hoisted resolution, zero per-call
 * allocation, scalar-expanded derivatives) are warranted and expected. What is
 * NOT relaxed is the project's truth discipline:
 *
 *   - EVERY parameter is resolved through the SAME helpers the readable terms
 *     use (bond_parameters, angle_parameters, stretch_bend_angle_terms,
 *     torsion_terms, oop_force_constant, nonbonded_context_for) — once, at
 *     build time. The fast loop only replays the arithmetic.
 *   - tests/fast-system.test.ts differentially checks this kernel against the
 *     readable terms: energies bit-for-bit, gradients to <=1e-8 absolute
 *     (the angle/dihedral/OOP derivative expansions below reassociate the
 *     reference's axis-loop arithmetic, which moves last-ULP rounding only),
 *     across fixtures, suite molecules, and randomly perturbed geometries.
 *
 * STRUCTURE. One FastSystem per prepared molecule (WeakMap cache — coordinates
 * may move between evaluations; chemistry may not, the optimizer's own
 * contract). Build resolves every interaction once into parallel typed arrays;
 * evaluate() walks those arrays with scalars only — no Vec3 tuples, no filter/
 * map, no Map lookups, nothing for the GC to see. Enumeration order matches the
 * readable terms' loops exactly (bonds in file order, angles/stretch-bends by
 * ascending center then adjacency order, torsions by central-bond order then
 * substituent adjacency order, OOP centers ascending, nonbonded pairs in
 * nonbonded-context order), so energy summation order — and therefore the
 * energy bits themselves — are unchanged.
 */

import type { TypedMolecule, EnergyComponents } from '../types.js';
import {
  class_context_for,
  bond_parameters,
  angle_parameters,
} from '../mmff94/parameters/parameter-classes.js';
import type { ClassContext } from '../mmff94/parameters/parameter-classes.js';
import { empirical_bond_parameters } from '../mmff94/parameters/empirical.js';
import {
  stretch_bend_angle_terms,
} from '../mmff94/energy/stretch-bend.js';
import { torsion_terms } from '../mmff94/energy/torsion.js';
import { oop_force_constant } from '../mmff94/energy/out-of-plane.js';
import { nonbonded_context_for } from '../mmff94/nonbonded-context.js';
import { assign_bci_charges } from '../mmff94/charges.js';

// ── constants — identical expressions to the readable terms ──────────
const BOND_UNIT = 143.9325;
const CS = -2.0;
const ANGLE_UNIT = 143.9325 * (Math.PI / 180) ** 2; // the published 0.043844, exact form
const LINEAR_UNIT = 143.9325;
const CB = -0.4 * (Math.PI / 180);
const SB_UNIT = 2.51210;
const OOP_UNIT = 143.9325 * (Math.PI / 180) ** 2;
const ELEC_UNIT = 332.0716;
const SCALE_1_4 = 0.75;
const ELEC_BUFFER = 0.05;
const RAD_PER_DEG = Math.PI / 180.0;

export interface FastSystem {
  readonly n_atoms: number;
  /**
   * Energy (+ optionally gradient) at `coords` (flat x,y,z per atom).
   * Zero allocations. After the call, `components` holds the per-term
   * breakdown and `total` the sum. `grad` (length 3·n_atoms) receives
   * dE/dx when non-null.
   *
   * `term_mask` (optional, 7 entries: bond, angle, strbnd, torsion,
   * oop, vdw, elec — a 0 skips that term entirely) exists for the
   * differential tests; the optimizer never passes it.
   */
  evaluate(coords: Float64Array, grad: Float64Array | null, term_mask?: Uint8Array): void;
  /** Valid immediately after evaluate(). Reused storage — copy if kept. */
  readonly components: EnergyComponents;
  readonly total: number;
}

// Resolved interaction tables (parallel arrays; counts in _n fields).
interface FastTables {
  // bonds: a1, a2, k_b, r0
  bond_n: number; bond_a1: Int32Array; bond_a2: Int32Array;
  bond_k: Float64Array; bond_r0: Float64Array;
  // angles (eq. 3/4): i-j-k, k_a, theta0(deg), linear flag
  ang_n: number; ang_i: Int32Array; ang_j: Int32Array; ang_k: Int32Array;
  ang_ka: Float64Array; ang_t0: Float64Array; ang_lin: Uint8Array;
  // stretch-bend (skips linear): i-j-k, k_ij, k_kj, r0_ij, r0_kj, theta0(deg)
  sb_n: number; sb_i: Int32Array; sb_j: Int32Array; sb_k: Int32Array;
  sb_kij: Float64Array; sb_kkj: Float64Array;
  sb_rij0: Float64Array; sb_rkj0: Float64Array; sb_t0: Float64Array;
  // torsions: i-j-k-l, V1..V3
  tor_n: number; tor_i: Int32Array; tor_j: Int32Array; tor_k: Int32Array; tor_l: Int32Array;
  tor_v1: Float64Array; tor_v2: Float64Array; tor_v3: Float64Array;
  // out-of-plane centers: j with substituents a,c,d and constant k
  oop_n: number; oop_j: Int32Array; oop_a: Int32Array; oop_c: Int32Array; oop_d: Int32Array;
  oop_k: Float64Array;
}

const systems = new WeakMap<TypedMolecule, FastSystemImpl>();

/** The fast system for a prepared molecule (built on first use). */
export function create_fast_system(prepared: TypedMolecule): FastSystem {
  let sys = systems.get(prepared);
  if (!sys) {
    sys = new FastSystemImpl(prepared);
    systems.set(prepared, sys);
  }
  return sys;
}

class FastSystemImpl implements FastSystem {
  readonly n_atoms: number;
  readonly components: EnergyComponents;
  total = 0;

  private t: FastTables;
  // nonbonded context (shared topology cache) + charge vector
  private nb: ReturnType<typeof nonbonded_context_for>;
  private q: Float64Array;

  constructor(prepared: TypedMolecule) {
    this.n_atoms = prepared.atoms.length;
    this.components = {
      total: 0, bond_stretch: 0, angle_bend: 0, stretch_bend: 0,
      torsion: 0, van_der_waals: 0, electrostatic: 0, out_of_plane: 0,
    };
    this.t = build_tables(prepared);
    this.nb = nonbonded_context_for(prepared);
    // BCI charges on demand, mirroring the readable term functions: a
    // bare typed molecule (no assign_bci_charges call) must not NaN the
    // electrostatic kernel with an empty charge array.
    this.q = Float64Array.from(
      prepared.partial_charges ?? assign_bci_charges(prepared).partial_charges!,
    );
  }

  evaluate(coords: Float64Array, grad: Float64Array | null, term_mask?: Uint8Array): void {
    const c = this.components;
    if (grad) grad.fill(0);
    const m = term_mask;
    c.bond_stretch   = (!m || m[0]) ? this.bonds(coords, grad) : 0;
    c.angle_bend     = (!m || m[1]) ? this.angles(coords, grad) : 0;
    c.stretch_bend   = (!m || m[2]) ? this.stretch_bend(coords, grad) : 0;
    c.torsion        = (!m || m[3]) ? this.torsions(coords, grad) : 0;
    c.out_of_plane   = (!m || m[4]) ? this.oop(coords, grad) : 0;
    c.van_der_waals  = (!m || m[5]) ? this.vdw(coords, grad) : 0;
    c.electrostatic  = (!m || m[6]) ? this.elec(coords, grad) : 0;
    c.total = c.bond_stretch + c.angle_bend + c.stretch_bend +
              c.torsion + c.van_der_waals + c.electrostatic + c.out_of_plane;
    this.total = c.total;
  }

  // ── bonds (eq. 2) ────────────────────────────────────────────────
  private bonds(x: Float64Array, g: Float64Array | null): number {
    const t = this.t;
    let e = 0;
    for (let p = 0; p < t.bond_n; p++) {
      const a = t.bond_a1[p], b = t.bond_a2[p];
      const dx = x[3*a] - x[3*b], dy = x[3*a+1] - x[3*b+1], dz = x[3*a+2] - x[3*b+2];
      const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
      const dr = r - t.bond_r0[p];
      const k = t.bond_k[p];
      const harmonic = BOND_UNIT * (0.5 * k) * dr * dr;
      const anharmonic = 1.0 + CS * dr + (7.0 / 12.0) * CS * CS * dr * dr;
      e += harmonic * anharmonic;

      if (g && r >= 1e-15) {
        // dE/dr = UNIT·k·Δr·[1 + cs·Δr + 7/12·cs²·Δr²]
        //         + UNIT·(k/2)·Δr²·[cs + 7/6·cs²·Δr];  dΔr/dx = ±r̂
        const d_anharmonic = CS + (7.0 / 6.0) * CS * CS * dr;
        const dE_dr = BOND_UNIT * k * dr * anharmonic +
                      BOND_UNIT * (k / 2.0) * dr * dr * d_anharmonic;
        const ux = dx / r, uy = dy / r, uz = dz / r;
        g[3*a] += dE_dr * ux; g[3*a+1] += dE_dr * uy; g[3*a+2] += dE_dr * uz;
        g[3*b] -= dE_dr * ux; g[3*b+1] -= dE_dr * uy; g[3*b+2] -= dE_dr * uz;
      }
    }
    return e;
  }

  // ── angles (eq. 3, or eq. 4 at lin-flagged centers) ──────────────
  private angles(x: Float64Array, g: Float64Array | null): number {
    const t = this.t;
    let e = 0;
    for (let p = 0; p < t.ang_n; p++) {
      const i = t.ang_i[p], j = t.ang_j[p], k = t.ang_k[p];
      const vix = x[3*i] - x[3*j], viy = x[3*i+1] - x[3*j+1], viz = x[3*i+2] - x[3*j+2];
      const vkx = x[3*k] - x[3*j], vky = x[3*k+1] - x[3*j+1], vkz = x[3*k+2] - x[3*j+2];
      const li = norm3(vix, viy, viz), lk = norm3(vkx, vky, vkz);
      // mirror vec_normalize's degenerate guard (zero vector → zero dot → π/2)
      const uix = li < 1e-15 ? 0 : vix / li, uiy = li < 1e-15 ? 0 : viy / li, uiz = li < 1e-15 ? 0 : viz / li;
      const ukx = lk < 1e-15 ? 0 : vkx / lk, uky = lk < 1e-15 ? 0 : vky / lk, ukz = lk < 1e-15 ? 0 : vkz / lk;
      let dot = uix*ukx + uiy*uky + uiz*ukz;
      if (dot > 1) dot = 1; else if (dot < -1) dot = -1;
      const theta_rad = Math.acos(dot);

      if (t.ang_lin[p]) {
        // eq. (4): E = LINEAR_UNIT·k·(1 + cos θ);  dE/dθ = −LINEAR_UNIT·k·sin θ
        const ka = t.ang_ka[p];
        e += LINEAR_UNIT * ka * (1.0 + Math.cos(theta_rad));
        if (g) {
          const dE_dth = -LINEAR_UNIT * ka * Math.sin(theta_rad);
          add_angle_grad(g, i, j, k, li, lk, uix, uiy, uiz, ukx, uky, ukz, dot, dE_dth, x);
        }
      } else {
        const theta_deg = theta_rad * (180.0 / Math.PI);
        const dth = theta_deg - t.ang_t0[p];
        const ka = t.ang_ka[p];
        const harmonic = ANGLE_UNIT * (0.5 * ka) * dth * dth;
        const anharmonic = 1.0 + CB * dth;
        e += harmonic * anharmonic;
        if (g) {
          const dE_ddeg = ANGLE_UNIT * ka * dth * anharmonic +
                          ANGLE_UNIT * (0.5 * ka) * dth * dth * CB;
          const dE_dth = dE_ddeg / RAD_PER_DEG;
          add_angle_grad(g, i, j, k, li, lk, uix, uiy, uiz, ukx, uky, ukz, dot, dE_dth, x);
        }
      }
    }
    return e;
  }

  // ── stretch-bend (eq. 5) ─────────────────────────────────────────
  private stretch_bend(x: Float64Array, g: Float64Array | null): number {
    const t = this.t;
    let e = 0;
    for (let p = 0; p < t.sb_n; p++) {
      const i = t.sb_i[p], j = t.sb_j[p], k = t.sb_k[p];
      const vix = x[3*i] - x[3*j], viy = x[3*i+1] - x[3*j+1], viz = x[3*i+2] - x[3*j+2];
      const vkx = x[3*k] - x[3*j], vky = x[3*k+1] - x[3*j+1], vkz = x[3*k+2] - x[3*j+2];
      const lij = norm3(vix, viy, viz), lkj = norm3(vkx, vky, vkz);
      const uix = vix / lij, uiy = viy / lij, uiz = viz / lij;
      const ukx = vkx / lkj, uky = vky / lkj, ukz = vkz / lkj;
      let dot = uix*ukx + uiy*uky + uiz*ukz;
      if (dot > 1) dot = 1; else if (dot < -1) dot = -1;
      const theta_rad = Math.acos(dot);
      const theta_deg = theta_rad * (180.0 / Math.PI);

      const dr_ij = lij - t.sb_rij0[p];
      const dr_kj = lkj - t.sb_rkj0[p];
      const dth = theta_deg - t.sb_t0[p];
      const kij = t.sb_kij[p], kkj = t.sb_kkj[p];
      e += SB_UNIT * (kij * dr_ij + kkj * dr_kj) * dth;

      if (g) {
        // dE/dx = SB·[k_ij·dΔr_ij/dx·Δθ + k_kj·dΔr_kj/dx·Δθ
        //             + (k_ij·Δr_ij + k_kj·Δr_kj)·dΔθ_deg/dx]
        const dE_ddth_deg = SB_UNIT * (kij * dr_ij + kkj * dr_kj);
        // bond-side pieces: dΔr/dx = ±r̂ (skip degenerate zero-length sides)
        if (lij >= 1e-15) {
          const f = SB_UNIT * kij * dth;
          g[3*i] += f * uix; g[3*i+1] += f * uiy; g[3*i+2] += f * uiz;
          g[3*j] -= f * uix; g[3*j+1] -= f * uiy; g[3*j+2] -= f * uiz;
        }
        if (lkj >= 1e-15) {
          const f = SB_UNIT * kkj * dth;
          g[3*k] += f * ukx; g[3*k+1] += f * uky; g[3*k+2] += f * ukz;
          g[3*j] -= f * ukx; g[3*j+1] -= f * uky; g[3*j+2] -= f * ukz;
        }
        // angle-side piece (degree-based): closed-form dθ/dx / RAD_PER_DEG
        add_angle_grad_deg(g, i, j, k, lij, lkj, uix, uiy, uiz, ukx, uky, ukz, dot, dE_ddth_deg, x);
      }
    }
    return e;
  }

  // ── torsions (eq. 7) ─────────────────────────────────────────────
  private torsions(x: Float64Array, g: Float64Array | null): number {
    const t = this.t;
    let e = 0;
    for (let p = 0; p < t.tor_n; p++) {
      const i = t.tor_i[p], j = t.tor_j[p], k = t.tor_k[p], l = t.tor_l[p];
      const v1x = x[3*i] - x[3*j], v1y = x[3*i+1] - x[3*j+1], v1z = x[3*i+2] - x[3*j+2];
      const v2x = x[3*k] - x[3*j], v2y = x[3*k+1] - x[3*j+1], v2z = x[3*k+2] - x[3*j+2];
      const v3x = x[3*l] - x[3*k], v3y = x[3*l+1] - x[3*k+1], v3z = x[3*l+2] - x[3*k+2];
      // n1 = v1×v2, n2 = v3×v2 (same handedness as dihedral_angle)
      const n1x = v1y*v2z - v1z*v2y, n1y = v1z*v2x - v1x*v2z, n1z = v1x*v2y - v1y*v2x;
      const n2x = v3y*v2z - v3z*v2y, n2y = v3z*v2x - v3x*v2z, n2z = v3x*v2y - v3y*v2x;
      const l1 = norm3(n1x, n1y, n1z), l2 = norm3(n2x, n2y, n2z), lv = norm3(v2x, v2y, v2z);
      let tau = 0;
      let s = 0, cc = 0;
      let ok = true;
      if (l1 < 1e-15 || l2 < 1e-15 || lv < 1e-15) {
        ok = false; // mirror dihedral_derivatives' degenerate zero-gradient guard
      } else {
        const h1x = n1x/l1, h1y = n1y/l1, h1z = n1z/l1;
        const h2x = n2x/l2, h2y = n2y/l2, h2z = n2z/l2;
        const hvx = v2x/lv, hvy = v2y/lv, hvz = v2z/lv;
        cc = h1x*h2x + h1y*h2y + h1z*h2z;
        // s = (n̂1×n̂2)·v̂2
        const cx = h1y*h2z - h1z*h2y, cy = h1z*h2x - h1x*h2z, cz = h1x*h2y - h1y*h2x;
        s = cx*hvx + cy*hvy + cz*hvz;
        tau = Math.atan2(s, cc);
      }
      const v1 = t.tor_v1[p], v2 = t.tor_v2[p], v3 = t.tor_v3[p];
      const cos1 = Math.cos(tau), cos2 = Math.cos(2.0*tau), cos3 = Math.cos(3.0*tau);
      e += 0.5 * v1 * (1.0 + cos1) +
           0.5 * v2 * (1.0 - cos2) +
           0.5 * v3 * (1.0 + cos3);

      if (g && ok) {
        const dE_dtau =
          -0.5 * v1 * Math.sin(tau) +
          v2 * Math.sin(2.0*tau) -
          1.5 * v3 * Math.sin(3.0*tau);
        add_dihedral_grad(g, i, j, k, l, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z,
                          n1x, n1y, n1z, n2x, n2y, n2z, l1, l2, lv, s, cc, dE_dtau);
      }
    }
    return e;
  }

  // ── out-of-plane (eq. 6) ─────────────────────────────────────────
  private oop(x: Float64Array, g: Float64Array | null): number {
    const t = this.t;
    let e = 0;
    for (let p = 0; p < t.oop_n; p++) {
      const j = t.oop_j[p], a = t.oop_a[p], c = t.oop_c[p], d = t.oop_d[p];
      const k = t.oop_k[p];
      const jx = x[3*j], jy = x[3*j+1], jz = x[3*j+2];
      // the three Wilson angles, same turns as the readable term
      const xa = wilson_deg(x[3*d]-jx, x[3*d+1]-jy, x[3*d+2]-jz,
                            x[3*c]-jx, x[3*c+1]-jy, x[3*c+2]-jz,
                            x[3*a]-jx, x[3*a+1]-jy, x[3*a+2]-jz);
      const xc = wilson_deg(x[3*a]-jx, x[3*a+1]-jy, x[3*a+2]-jz,
                            x[3*d]-jx, x[3*d+1]-jy, x[3*d+2]-jz,
                            x[3*c]-jx, x[3*c+1]-jy, x[3*c+2]-jz);
      const xd = wilson_deg(x[3*a]-jx, x[3*a+1]-jy, x[3*a+2]-jz,
                            x[3*c]-jx, x[3*c+1]-jy, x[3*c+2]-jz,
                            x[3*d]-jx, x[3*d+1]-jy, x[3*d+2]-jz);
      e += OOP_UNIT * (k / 2.0) * (xa*xa + xc*xc + xd*xd);

      if (g) {
        // dE/dx = OOP_UNIT·k·χ·dχ_deg/dx per turn (χ in degrees here)
        const dE = OOP_UNIT * k;
        add_oop_grad(g, d, j, c, a, xa * dE, x);
        add_oop_grad(g, a, j, d, c, xc * dE, x);
        add_oop_grad(g, a, j, c, d, xd * dE, x);
      }
    }
    return e;
  }

  // ── van der Waals (buffered 14-7, eq. 8) ─────────────────────────
  private vdw(x: Float64Array, g: Float64Array | null): number {
    const nb = this.nb;
    let e = 0;
    for (let p = 0; p < nb.n_pairs; p++) {
      const eps = nb.pair_epsilon_ij[p];
      if (isNaN(eps)) continue;
      const i = nb.pair_i[p], j = nb.pair_j[p];
      const dx = x[3*i] - x[3*j], dy = x[3*i+1] - x[3*j+1], dz = x[3*i+2] - x[3*j+2];
      const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
      const a = nb.pair_vdw_a[p], b = nb.pair_vdw_b[p];
      const C = nb.pair_vdw_C[p], D = nb.pair_vdw_D[p];
      const r7 = Math.pow(r, 7);
      const f_rep = Math.pow(a / (r + b), 7);
      const f_att = C / (r7 + D) - 2;
      e += eps * f_rep * f_att;

      if (g) {
        const f_rep_p = -7.0 * Math.pow(a, 7) / Math.pow(r + b, 8);
        const f_att_p = -7.0 * C * Math.pow(r, 6) / Math.pow(r7 + D, 2);
        const dE_dr = eps * (f_rep_p * f_att + f_rep * f_att_p);
        if (r >= 1e-15) {
          const ux = dx / r, uy = dy / r, uz = dz / r;
          g[3*i] += dE_dr * ux; g[3*i+1] += dE_dr * uy; g[3*i+2] += dE_dr * uz;
          g[3*j] -= dE_dr * ux; g[3*j+1] -= dE_dr * uy; g[3*j+2] -= dE_dr * uz;
        }
      }
    }
    return e;
  }

  // ── electrostatics (part III eq. 6, 1-4 ×0.75) ───────────────────
  private elec(x: Float64Array, g: Float64Array | null): number {
    const nb = this.nb;
    const q = this.q;
    let e = 0;
    for (let p = 0; p < nb.n_pairs; p++) {
      const i = nb.pair_i[p], j = nb.pair_j[p];
      const qi = q[i], qj = q[j];
      if (qi === 0 || qj === 0) continue;
      const dx = x[3*i] - x[3*j], dy = x[3*i+1] - x[3*j+1], dz = x[3*i+2] - x[3*j+2];
      const rb = Math.sqrt(dx*dx + dy*dy + dz*dz) + ELEC_BUFFER;
      let pair_e = (ELEC_UNIT * qi * qj) / rb;
      if (nb.pair_is_14[p]) pair_e *= SCALE_1_4;
      e += pair_e;

      if (g) {
        let dE_dr = -ELEC_UNIT * qi * qj / (rb * rb);
        if (nb.pair_is_14[p]) dE_dr *= SCALE_1_4;
        const r = rb - ELEC_BUFFER;
        if (r >= 1e-15) {
          const ux = dx / r, uy = dy / r, uz = dz / r;
          g[3*i] += dE_dr * ux; g[3*i+1] += dE_dr * uy; g[3*i+2] += dE_dr * uz;
          g[3*j] -= dE_dr * ux; g[3*j+1] -= dE_dr * uy; g[3*j+2] -= dE_dr * uz;
        }
      }
    }
    return e;
  }
}

// ── build: resolve every interaction once ────────────────────────────
function build_tables(mol: TypedMolecule): FastTables {
  const n = mol.atoms.length;
  const ctx: ClassContext = (() => {
    const adj: number[][] = Array.from({ length: n }, () => []);
    for (const b of mol.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }
    return class_context_for(mol, adj);
  })();
  // adjacency in the same construction order the readable terms use
  const adj: number[][] = Array.from({ length: n }, () => []);
  for (const b of mol.bonds) { adj[b.atom1].push(b.atom2); adj[b.atom2].push(b.atom1); }

  // bonds (file order; empirical-rule fallback mirrors the readable term)
  const nb_bonds = mol.bonds.length;
  const b_a1 = new Int32Array(nb_bonds), b_a2 = new Int32Array(nb_bonds);
  const b_k = new Float64Array(nb_bonds), b_r0 = new Float64Array(nb_bonds);
  let bn = 0;
  for (const bond of mol.bonds) {
    let params = bond_parameters(ctx, bond.atom1, bond.atom2);
    if (!params) params = empirical_bond_parameters(mol.atoms[bond.atom1], mol.atoms[bond.atom2]);
    if (!params) continue;
    b_a1[bn] = bond.atom1; b_a2[bn] = bond.atom2;
    b_k[bn] = params.k_b; b_r0[bn] = params.r0;
    bn++;
  }

  // count angles first (j ascending, adjacency-order pairs)
  let n_ang = 0;
  for (let j = 0; j < n; j++) {
    const nbrs = adj[j];
    n_ang += nbrs.length * (nbrs.length - 1) / 2;
  }
  const a_i = new Int32Array(n_ang), a_j = new Int32Array(n_ang), a_k = new Int32Array(n_ang);
  const a_ka = new Float64Array(n_ang), a_t0 = new Float64Array(n_ang), a_lin = new Uint8Array(n_ang);
  const s_i = new Int32Array(n_ang), s_j = new Int32Array(n_ang), s_k = new Int32Array(n_ang);
  const s_kij = new Float64Array(n_ang), s_kkj = new Float64Array(n_ang);
  const s_rij0 = new Float64Array(n_ang), s_rkj0 = new Float64Array(n_ang), s_t0 = new Float64Array(n_ang);
  let an = 0, sn = 0;
  for (let j = 0; j < n; j++) {
    const nbrs = adj[j];
    for (let ii = 0; ii < nbrs.length; ii++) {
      for (let kk = ii + 1; kk < nbrs.length; kk++) {
        const i = nbrs[ii], k = nbrs[kk];
        const ap = angle_parameters(ctx, i, j, k);
        a_i[an] = i; a_j[an] = j; a_k[an] = k;
        a_ka[an] = ap.k_a; a_t0[an] = ap.theta0; a_lin[an] = ap.linear ? 1 : 0;
        an++;
        const sb = stretch_bend_angle_terms(ctx, mol, i, j, k);
        if (sb && !sb.linear) {
          s_i[sn] = i; s_j[sn] = j; s_k[sn] = k;
          s_kij[sn] = sb.k_ij; s_kkj[sn] = sb.k_kj;
          s_rij0[sn] = sb.r0_ij; s_rkj0[sn] = sb.r0_kj; s_t0[sn] = sb.theta0;
          sn++;
        }
      }
    }
  }

  // torsions (central-bond order, then substituent adjacency order)
  const n_tor_cap = mol.bonds.length * 9; // crude cap; grown via plain arrays then copied
  const t_i: number[] = [], t_j: number[] = [], t_k: number[] = [], t_l: number[] = [];
  const t_v1: number[] = [], t_v2: number[] = [], t_v3: number[] = [];
  void n_tor_cap;
  for (const bond of mol.bonds) {
    const j = bond.atom1, k = bond.atom2;
    const ins = adj[j].filter(nn => nn !== k);
    const lns = adj[k].filter(nn => nn !== j);
    if (ins.length === 0 || lns.length === 0) continue;
    for (const i of ins) {
      for (const l of lns) {
        if (l === i) continue; // 3-ring closure — BatchMin's d == a rule
        const terms = torsion_terms(ctx, mol, i, j, k, l);
        if (!terms) continue;
        t_i.push(i); t_j.push(j); t_k.push(k); t_l.push(l);
        t_v1.push(terms.v1); t_v2.push(terms.v2); t_v3.push(terms.v3);
      }
    }
  }

  // out-of-plane centers (exactly three neighbors)
  const oop_j: number[] = [], oop_a: number[] = [], oop_c: number[] = [], oop_d: number[] = [];
  const oop_k: number[] = [];
  for (let j = 0; j < n; j++) {
    const nbrs = adj[j];
    if (nbrs.length !== 3) continue;
    const kc = oop_force_constant(mol, j, nbrs[0], nbrs[1], nbrs[2]);
    if (kc === undefined) continue;
    oop_j.push(j); oop_a.push(nbrs[0]); oop_c.push(nbrs[1]); oop_d.push(nbrs[2]);
    oop_k.push(kc);
  }

  return {
    bond_n: bn, bond_a1: b_a1, bond_a2: b_a2, bond_k: b_k, bond_r0: b_r0,
    ang_n: an, ang_i: a_i, ang_j: a_j, ang_k: a_k,
    ang_ka: a_ka, ang_t0: a_t0, ang_lin: a_lin,
    sb_n: sn, sb_i: s_i, sb_j: s_j, sb_k: s_k,
    sb_kij: s_kij, sb_kkj: s_kkj, sb_rij0: s_rij0, sb_rkj0: s_rkj0, sb_t0: s_t0,
    tor_n: t_i.length,
    tor_i: Int32Array.from(t_i), tor_j: Int32Array.from(t_j),
    tor_k: Int32Array.from(t_k), tor_l: Int32Array.from(t_l),
    tor_v1: Float64Array.from(t_v1), tor_v2: Float64Array.from(t_v2), tor_v3: Float64Array.from(t_v3),
    oop_n: oop_j.length,
    oop_j: Int32Array.from(oop_j), oop_a: Int32Array.from(oop_a),
    oop_c: Int32Array.from(oop_c), oop_d: Int32Array.from(oop_d),
    oop_k: Float64Array.from(oop_k),
  };
}

// ── scalar geometry kernel helpers (module-scope, monomorphic) ───────

/** Euclidean length of (x,y,z). */
function norm3(x: number, y: number, z: number): number {
  return Math.sqrt(x*x + y*y + z*z);
}

/**
 * Wilson χ in DEGREES for plane atoms (i,j) + out-of-plane l — the flat
 * transcription of wilson_oop_angle (unit vectors, n̂ = û_ji × û_jk,
 * χ = asin(clamp(n̂·û_jl))). Guard mirrors: zero-length or collinear → 0.
 */
function wilson_deg(
  vix: number, viy: number, viz: number,
  vkx: number, vky: number, vkz: number,
  vlx: number, vly: number, vlz: number,
): number {
  const li = norm3(vix, viy, viz), lk = norm3(vkx, vky, vkz), ll = norm3(vlx, vly, vlz);
  if (li < 1e-15 || lk < 1e-15 || ll < 1e-15) return 0;
  const ux = vix/li, uy = viy/li, uz = viz/li;
  const wx = vkx/lk, wy = vky/lk, wz = vkz/lk;
  const nx = uy*wz - uz*wy, ny = uz*wx - ux*wz, nz = ux*wy - uy*wx;
  const ln = norm3(nx, ny, nz);
  if (ln < 1e-15) return 0;
  let s = (nx/ln)*(vlx/ll) + (ny/ln)*(vly/ll) + (nz/ln)*(vlz/ll);
  if (s > 1) s = 1; else if (s < -1) s = -1;
  return Math.asin(s) * (180.0 / Math.PI);
}

/**
 * Angle-gradient accumulator (RADIANS form): adds dE/dθ · dθ/dx for the
 * angle i-j-k. Closed-form transcription of angle_derivatives:
 *   dθ/dx_i[α] = −(û_jk[α] − û_ji[α]·cos)/(sinθ·|v_ji|)
 *   dθ/dx_k[α] = −(û_ji[α] − û_jk[α]·cos)/(sinθ·|v_jk|)
 *   dθ/dx_j[α] = −Σ (the reference builds j's piece the same way)
 * The reference falls back to a forward-difference of θ when sin θ < 1e-12
 * (cusps at 0°/180°); the scalar FD fallback below mirrors that.
 */
function add_angle_grad(
  g: Float64Array, i: number, j: number, k: number,
  li: number, lk: number,
  uix: number, uiy: number, uiz: number,
  ukx: number, uky: number, ukz: number,
  cos_t: number, dE_dtheta: number,
  x: Float64Array,
): void {
  add_angle_grad_deg(g, i, j, k, li, lk,
                     uix, uiy, uiz, ukx, uky, ukz, cos_t, dE_dtheta * RAD_PER_DEG, x);
}

/** Degree-form variant: `dE_dtheta_deg` multiplies dθ_deg/dx (= dθ/dx / RAD_PER_DEG). */
function add_angle_grad_deg(
  g: Float64Array, i: number, j: number, k: number,
  li: number, lk: number,
  uix: number, uiy: number, uiz: number,
  ukx: number, uky: number, ukz: number,
  cos_t: number, dE_dtheta_deg: number,
  x: Float64Array,
): void {
  const sin_t = Math.sqrt(Math.max(0.0, 1.0 - cos_t * cos_t));
  if (li < 1e-15 || lk < 1e-15) return; // zero-length arm: reference returns zeros via guards
  if (sin_t < 1e-12) {
    fd_angle_grad(g, i, j, k, dE_dtheta_deg, x);
    return;
  }
  // dθ/dx_i[α] = −(û_jk[α] − û_ji[α]·cos)/(sin·|v_ji|) — the closed form
  // of the reference's unit_vec_deriv chain; j's piece is minus the sum.
  const si = -dE_dtheta_deg / (sin_t * li * RAD_PER_DEG);
  const sk = -dE_dtheta_deg / (sin_t * lk * RAD_PER_DEG);
  const gix = si * (ukx - uix * cos_t);
  const giy = si * (uky - uiy * cos_t);
  const giz = si * (ukz - uiz * cos_t);
  const gkx = sk * (uix - ukx * cos_t);
  const gky = sk * (uiy - uky * cos_t);
  const gkz = sk * (uiz - ukz * cos_t);
  g[3*i] += gix; g[3*i+1] += giy; g[3*i+2] += giz;
  g[3*k] += gkx; g[3*k+1] += gky; g[3*k+2] += gkz;
  g[3*j] -= gix + gkx; g[3*j+1] -= giy + gky; g[3*j+2] -= giz + gkz;
}

// Scratch for the FD fallbacks — evaluate() never re-enters, so shared
// module-level buffers are safe and keep the fallback allocation-free.
const FD_D = new Float64Array(9); // 3 atoms × 3 axes of d(angle)/dx (radians)

/** Forward-difference fallback for the angle gradient (mirrors the
 *  readable term's angle_derivatives_fd: cusped at exactly 0°/180°). */
function fd_angle_grad(
  g: Float64Array, i: number, j: number, k: number,
  dE_dtheta_deg: number, x: Float64Array,
): void {
  const DELTA = 1e-6;
  const th = (px: number, py: number, pz: number) => {
    const ax = px - x[3*j], ay = py - x[3*j+1], az = pz - x[3*j+2];
    const bx = x[3*k] - x[3*j], by = x[3*k+1] - x[3*j+1], bz = x[3*k+2] - x[3*j+2];
    const la = norm3(ax, ay, az), lb = norm3(bx, by, bz);
    if (la < 1e-15 || lb < 1e-15) return Math.acos(0);
    let dd = (ax/la)*(bx/lb) + (ay/lb===0?0:ay/lb)*0; // placeholder, replaced below
    void dd;
    const ux = ax/la, uy = ay/la, uz = az/la;
    const vx = bx/lb, vy = by/lb, vz = bz/lb;
    dd = ux*vx + uy*vy + uz*vz;
    if (dd > 1) dd = 1; else if (dd < -1) dd = -1;
    return Math.acos(dd);
  };
  const atoms = [i, j, k];
  for (let n = 0; n < 3; n++) {
    const a = atoms[n];
    for (let axis = 0; axis < 3; axis++) {
      const orig = x[3*a + axis];
      x[3*a + axis] = orig + DELTA;
      // θ(i-j-k) with j as center: pass i's position; k read live
      const t1 = th(x[3*i], x[3*i+1], x[3*i+2]);
      x[3*a + axis] = orig;
      // forward difference of the angle itself (radians per Å)
      FD_D[n*3 + axis] = (t1 - th(x[3*i], x[3*i+1], x[3*i+2])) / DELTA;
    }
  }
  // NOTE: the lambda above recomputes BASE theta after restore, giving the
  // forward difference (moved − base) — matches the reference's convention.
  const scale = dE_dtheta_deg / RAD_PER_DEG;
  for (let n = 0; n < 3; n++) {
    const a = atoms[n];
    g[3*a]   += scale * FD_D[n*3];
    g[3*a+1] += scale * FD_D[n*3 + 1];
    g[3*a+2] += scale * FD_D[n*3 + 2];
  }
}

/**
 * Dihedral-gradient accumulator: scalar expansion of the readable term's
 * dihedral_derivatives — same construction (n1 = v1×v2, n2 = v3×v2,
 * τ = atan2(s,c)), same per-axis contribution algebra, evaluated with
 * scalars instead of tuple allocations. Degenerate geometries (any length
 * < 1e-15) were already screened by the caller (zero contribution).
 */
function add_dihedral_grad(
  g: Float64Array, i: number, j: number, k: number, l: number,
  v1x: number, v1y: number, v1z: number,
  v2x: number, v2y: number, v2z: number,
  v3x: number, v3y: number, v3z: number,
  n1x: number, n1y: number, n1z: number,
  n2x: number, n2y: number, n2z: number,
  len1: number, len2: number, lenv: number,
  s: number, c: number, dE_dtau: number,
): void {
  // hats
  const h1x = n1x/len1, h1y = n1y/len1, h1z = n1z/len1;
  const h2x = n2x/len2, h2y = n2y/len2, h2z = n2z/len2;
  const hvx = v2x/lenv, hvy = v2y/lenv, hvz = v2z/lenv;
  // n̂1×n̂2 (constant across axes)
  const nx2x = h1y*h2z - h1z*h2y, nx2y = h1z*h2x - h1x*h2z, nx2z = h1x*h2y - h1y*h2x;

  // unit_vec_deriv(v, dv): ((dv − û(û·dv))/|v|) — inlined per use below.

  // per-axis dτ/dx for the current axis (the loop variable IS the axis)
  let dix = 0, djx = 0, dkx = 0, dlx = 0;

  for (let ax = 0; ax < 3; ax++) {
    // basis-vector cross products with v2 (sparse): e×v2 and the dv patterns
    // e1×v2 = (0,-v2z,v2y); e2×v2 = (v2z,0,-v2x); e3×v2 = (-v2y,v2x,0)
    let ev2x = 0, ev2y = 0, ev2z = 0;
    if (ax === 0) { ev2y = -v2z; ev2z = v2y; }
    else if (ax === 1) { ev2x = v2z; ev2z = -v2x; }
    else { ev2x = -v2y; ev2y = v2x; }

    // ── atom i: dv1 = e, dv2 = 0, dv3 = 0 → dn1 = e×v2, dn2 = 0
    {
      // dn̂1 = (dn1 − n̂1(n̂1·dn1))/|n1|
      const d = h1x*ev2x + h1y*ev2y + h1z*ev2z;
      const u1x = (ev2x - h1x*d)/len1, u1y = (ev2y - h1y*d)/len1, u1z = (ev2z - h1z*d)/len1;
      const dc = u1x*h2x + u1y*h2y + u1z*h2z;
      // ds = (dn̂1×n̂2)·v̂2
      const cx = u1y*h2z - u1z*h2y, cy = u1z*h2x - u1x*h2z, cz = u1x*h2y - u1y*h2x;
      const ds = cx*hvx + cy*hvy + cz*hvz;
      const dtau = c*ds - s*dc;
      dix = dtau;
    }
    // ── atom j: dv1 = −e, dv2 = −e, dv3 = 0 → dn1 = (−e)×v2 + v1×(−e), dn2 = v3×(−e)
    {
      // v×w = (vy wz − vz wy, vz wx − vx wz, vx wy − vy wx)
      let mvx = 0, mvy = 0, mvz = 0; // the −e vector
      if (ax === 0) { mvx = -1; } else if (ax === 1) { mvy = -1; } else { mvz = -1; }
      const dn1x = mvy*v2z - mvz*v2y + v1y*mvz - v1z*mvy;
      const dn1y = mvz*v2x - mvx*v2z + v1z*mvx - v1x*mvz;
      const dn1z = mvx*v2y - mvy*v2x + v1x*mvy - v1y*mvx;
      const dn2x = v3y*mvz - v3z*mvy;
      const dn2y = v3z*mvx - v3x*mvz;
      const dn2z = v3x*mvy - v3y*mvx;
      // dn̂1 = unit_deriv(n1, dn1); dn̂2 = unit_deriv(n2, dn2); dv̂2 = unit_deriv(v2, dv2=mv)
      const d1 = h1x*dn1x + h1y*dn1y + h1z*dn1z;
      const un1x = (dn1x - h1x*d1)/len1, un1y = (dn1y - h1y*d1)/len1, un1z = (dn1z - h1z*d1)/len1;
      const d2 = h2x*dn2x + h2y*dn2y + h2z*dn2z;
      const un2x = (dn2x - h2x*d2)/len2, un2y = (dn2y - h2y*d2)/len2, un2z = (dn2z - h2z*d2)/len2;
      const dv = mvx*hvx + mvy*hvy + mvz*hvz;
      const uvx = (mvx - hvx*dv)/lenv, uvy = (mvy - hvy*dv)/lenv, uvz = (mvz - hvz*dv)/lenv;
      const dc = un1x*h2x + un1y*h2y + un1z*h2z + h1x*un2x + h1y*un2y + h1z*un2z;
      const cx = un1y*h2z - un1z*h2y, cy = un1z*h2x - un1x*h2z, cz = un1x*h2y - un1y*h2x;
      const cx2 = h1y*un2z - h1z*un2y, cy2 = h1z*un2x - h1x*un2z, cz2 = h1x*un2y - h1y*un2x;
      const ds = (cx + cx2)*hvx + (cy + cy2)*hvy + (cz + cz2)*hvz + nx2x*uvx + nx2y*uvy + nx2z*uvz;
      const dtau = c*ds - s*dc;
      djx = dtau;
    }
    // ── atom k: dv1 = 0, dv2 = e, dv3 = −e → dn1 = v1×e, dn2 = (−e)×v2 + v3×e
    {
      let ex = 0, ey = 0, ez = 0;
      if (ax === 0) ex = 1; else if (ax === 1) ey = 1; else ez = 1;
      const dn1x = v1y*ez - v1z*ey, dn1y = v1z*ex - v1x*ez, dn1z = v1x*ey - v1y*ex;
      // (−e)×v2 + v3×e — both crosses needed: k moves BOTH v2 (+e) and v3 (−e)
      const dn2x = ez*v2y - ey*v2z + v3y*ez - v3z*ey;
      const dn2y = ex*v2z - ez*v2x + v3z*ex - v3x*ez;
      const dn2z = ey*v2x - ex*v2y + v3x*ey - v3y*ex;
      const d1 = h1x*dn1x + h1y*dn1y + h1z*dn1z;
      const un1x = (dn1x - h1x*d1)/len1, un1y = (dn1y - h1y*d1)/len1, un1z = (dn1z - h1z*d1)/len1;
      const d2 = h2x*dn2x + h2y*dn2y + h2z*dn2z;
      const un2x = (dn2x - h2x*d2)/len2, un2y = (dn2y - h2y*d2)/len2, un2z = (dn2z - h2z*d2)/len2;
      const dv = ex*hvx + ey*hvy + ez*hvz;
      const uvx = (ex - hvx*dv)/lenv, uvy = (ey - hvy*dv)/lenv, uvz = (ez - hvz*dv)/lenv;
      const dc = un1x*h2x + un1y*h2y + un1z*h2z + h1x*un2x + h1y*un2y + h1z*un2z;
      const cx = un1y*h2z - un1z*h2y, cy = un1z*h2x - un1x*h2z, cz = un1x*h2y - un1y*h2x;
      const cx2 = h1y*un2z - h1z*un2y, cy2 = h1z*un2x - h1x*un2z, cz2 = h1x*un2y - h1y*un2x;
      const ds = (cx + cx2)*hvx + (cy + cy2)*hvy + (cz + cz2)*hvz + nx2x*uvx + nx2y*uvy + nx2z*uvz;
      const dtau = c*ds - s*dc;
      dkx = dtau;
    }
    // ── atom l: dv3 = e → dn2 = e×v2
    {
      let dn2x = 0, dn2y = 0, dn2z = 0;
      if (ax === 0) { dn2y = -v2z; dn2z = v2y; }
      else if (ax === 1) { dn2x = v2z; dn2z = -v2x; }
      else { dn2x = -v2y; dn2y = v2x; }
      const d2 = h2x*dn2x + h2y*dn2y + h2z*dn2z;
      const un2x = (dn2x - h2x*d2)/len2, un2y = (dn2y - h2y*d2)/len2, un2z = (dn2z - h2z*d2)/len2;
      const dc = h1x*un2x + h1y*un2y + h1z*un2z;
      const cx2 = h1y*un2z - h1z*un2y, cy2 = h1z*un2x - h1x*un2z, cz2 = h1x*un2y - h1y*un2x;
      const ds = cx2*hvx + cy2*hvy + cz2*hvz;
      const dtau = c*ds - s*dc;
      dlx = dtau;
    }

    // commit this axis's pieces (scaled by dE/dτ)
    g[3*i + ax] += dE_dtau * dix;
    g[3*j + ax] += dE_dtau * djx;
    g[3*k + ax] += dE_dtau * dkx;
    g[3*l + ax] += dE_dtau * dlx;
  }
}

/**
 * Out-of-plane gradient accumulator for one Wilson turn (atoms i-j center-
 * j-k plane, l out of plane): scalar expansion of oop_angle_derivatives,
 * scaled by dE/dχ_deg (which folds in the /RAD_PER_DEG degree conversion).
 */
function add_oop_grad(
  g: Float64Array, i: number, j: number, k: number, l: number,
  dE_dchi_deg: number, x: Float64Array,
): void {
  const vix = x[3*i] - x[3*j], viy = x[3*i+1] - x[3*j+1], viz = x[3*i+2] - x[3*j+2];
  const vkx = x[3*k] - x[3*j], vky = x[3*k+1] - x[3*j+1], vkz = x[3*k+2] - x[3*j+2];
  const vlx = x[3*l] - x[3*j], vly = x[3*l+1] - x[3*j+1], vlz = x[3*l+2] - x[3*j+2];
  const li = norm3(vix, viy, viz), lk = norm3(vkx, vky, vkz), ll = norm3(vlx, vly, vlz);
  if (li < 1e-15 || lk < 1e-15 || ll < 1e-15) return;
  const uix = vix/li, uiy = viy/li, uiz = viz/li;
  const ukx = vkx/lk, uky = vky/lk, ukz = vkz/lk;
  const ulx = vlx/ll, uly = vly/ll, ulz = vlz/ll;
  const nx = uiy*ukz - uiz*uky, ny = uiz*ukx - uix*ukz, nz = uix*uky - uiy*ukx;
  const ln = norm3(nx, ny, nz);
  if (ln < 1e-15) return; // collinear: χ ≡ 0, zero gradient (mirrors the reference)
  const hx = nx/ln, hy = ny/ln, hz = nz/ln;
  const sc = hx*ulx + hy*uly + hz*ulz;
  const cos_chi = Math.sqrt(Math.max(0.0, 1.0 - sc*sc));
  if (cos_chi < 1e-12) {
    fd_oop_grad(g, i, j, k, l, dE_dchi_deg, x);
    return;
  }

  // per-axis contribution(): dn = du_ji×u_jk + u_ji×du_jk; dn̂ = ud(n, dn);
  // ds = dn̂·u_jl + n̂·du_jl; dχ = ds/cosχ. Each atom branch computes the
  // dχ/dx for THIS axis only (the loop variable is the Cartesian axis).
  for (let ax = 0; ax < 3; ax++) {
    let ri = 0, rj = 0, rk = 0, rl = 0;
    // atom i: du_ji = ud(v_ji, e) — only u_ji moves, so dn = du_ji×u_jk
    {
      const ex = ax === 0 ? 1 : 0, ey = ax === 1 ? 1 : 0, ez = ax === 2 ? 1 : 0;
      const d = uix*ex + uiy*ey + uiz*ez;
      const bx = (ex - uix*d)/li, by = (ey - uiy*d)/li, bz = (ez - uiz*d)/li;
      const dnx = by*ukz - bz*uky;
      const dny = bz*ukx - bx*ukz;
      const dnz = bx*uky - by*ukx;
      const dh = hx*dnx + hy*dny + hz*dnz;
      const hnx = (dnx - hx*dh)/ln, hny = (dny - hy*dh)/ln, hnz = (dnz - hz*dh)/ln;
      ri = ((hnx*ulx + hny*uly + hnz*ulz) / cos_chi);
    }
    // atom k: du_jk = ud(v_jk, e)
    {
      const ex = ax === 0 ? 1 : 0, ey = ax === 1 ? 1 : 0, ez = ax === 2 ? 1 : 0;
      const d = ukx*ex + uky*ey + ukz*ez;
      const bx = (ex - ukx*d)/lk, by = (ey - uky*d)/lk, bz = (ez - ukz*d)/lk;
      const dnx = uiy*bz - uiz*by;
      const dny = uiz*bx - uix*bz;
      const dnz = uix*by - uiy*bx;
      const dh = hx*dnx + hy*dny + hz*dnz;
      const hnx = (dnx - hx*dh)/ln, hny = (dny - hy*dh)/ln, hnz = (dnz - hz*dh)/ln;
      rk = ((hnx*ulx + hny*uly + hnz*ulz) / cos_chi);
    }
    // atom l: du_jl = ud(v_jl, e)
    {
      const ex = ax === 0 ? 1 : 0, ey = ax === 1 ? 1 : 0, ez = ax === 2 ? 1 : 0;
      const d = ulx*ex + uly*ey + ulz*ez;
      const bx = (ex - ulx*d)/ll, by = (ey - uly*d)/ll, bz = (ez - ulz*d)/ll;
      rl = ((hx*bx + hy*by + hz*bz) / cos_chi);
    }
    // atom j: du_ji = ud(v_ji, −e), du_jk = ud(v_jk, −e), du_jl = ud(v_jl, −e)
    {
      const ex = ax === 0 ? -1 : 0, ey = ax === 1 ? -1 : 0, ez = ax === 2 ? -1 : 0;
      const di = uix*ex + uiy*ey + uiz*ez;
      const bix = (ex - uix*di)/li, biy = (ey - uiy*di)/li, biz = (ez - uiz*di)/li;
      const dk = ukx*ex + uky*ey + ukz*ez;
      const bkx = (ex - ukx*dk)/lk, bky = (ey - uky*dk)/lk, bkz = (ez - ukz*dk)/lk;
      const dl = ulx*ex + uly*ey + ulz*ez;
      const blx = (ex - ulx*dl)/ll, bly = (ey - uly*dl)/ll, blz = (ez - ulz*dl)/ll;
      const dnx = biy*ukz - biz*uky + uiy*bkz - uiz*bky;
      const dny = biz*ukx - bix*ukz + uiz*bkx - uix*bkz;
      const dnz = bix*uky - biy*ukx + uix*bky - uiy*bkx;
      const dh = hx*dnx + hy*dny + hz*dnz;
      const hnx = (dnx - hx*dh)/ln, hny = (dny - hy*dh)/ln, hnz = (dnz - hz*dh)/ln;
      rj = ((hnx*ulx + hny*uly + hnz*ulz + hx*blx + hy*bly + hz*blz) / cos_chi);
    }
    const f = dE_dchi_deg / RAD_PER_DEG;
    g[3*i + ax] += f * ri;
    g[3*j + ax] += f * rj;
    g[3*k + ax] += f * rk;
    g[3*l + ax] += f * rl;
  }
}

/** Forward-difference fallback for the OOP gradient near χ = ±90°. */
function fd_oop_grad(
  g: Float64Array, i: number, j: number, k: number, l: number,
  dE_dchi_deg: number, x: Float64Array,
): void {
  const DELTA = 1e-6;
  const atoms = [i, j, k, l];
  for (let n = 0; n < 4; n++) {
    const a = atoms[n];
    for (let axis = 0; axis < 3; axis++) {
      const orig = x[3*a + axis];
      x[3*a + axis] = orig + DELTA;
      const moved = wilson_deg(
        x[3*i]-x[3*j], x[3*i+1]-x[3*j+1], x[3*i+2]-x[3*j+2],
        x[3*k]-x[3*j], x[3*k+1]-x[3*j+1], x[3*k+2]-x[3*j+2],
        x[3*l]-x[3*j], x[3*l+1]-x[3*j+1], x[3*l+2]-x[3*j+2],
      );
      x[3*a + axis] = orig;
      const base = wilson_deg(
        x[3*i]-x[3*j], x[3*i+1]-x[3*j+1], x[3*i+2]-x[3*j+2],
        x[3*k]-x[3*j], x[3*k+1]-x[3*j+1], x[3*k+2]-x[3*j+2],
        x[3*l]-x[3*j], x[3*l+1]-x[3*j+1], x[3*l+2]-x[3*j+2],
      );
      g[3*a + axis] += dE_dchi_deg * ((moved - base) / DELTA) / RAD_PER_DEG;
    }
  }
}
