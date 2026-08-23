/**
 * Analytical derivatives of the geometric quantities the energy terms
 * depend on: bond length, bond angle, dihedral angle, and the Wilson
 * out-of-plane angle.
 *
 * Every derivative here is the chain rule through the SAME computational
 * path the energy term uses (see src/utils/vector.ts): the angle helpers
 * differentiate the normalized bond vectors exactly as
 * angle_in_radians() / dihedral_angle() / wilson_oop_angle() build them.
 * That is what makes the analytical gradient agree with the finite
 * difference of the energy term by construction — a derivative of a
 * differently-ordered expression would still be mathematically equal,
 * but floating-point rounding would show up as test noise.
 *
 * All angles are returned in RADIANS per Ångström. The energy terms
 * convert their degree-based constants themselves.
 */

import {
  Vec3,
  vec_sub,
  vec_dot,
  vec_cross,
  vec_length,
} from '../../utils/vector.js';

/** Derivative of the unit vector u = v/|v| with respect to a displacement dv of v. */
export function unit_vec_deriv(v: Vec3, dv: Vec3): Vec3 {
  const len = vec_length(v);
  if (len < 1e-15) return [0, 0, 0];
  // d(v/|v|) = (I − ûûᵀ)/|v| · dv
  const u = [v[0] / len, v[1] / len, v[2] / len];
  const dot = u[0] * dv[0] + u[1] * dv[1] + u[2] * dv[2];
  return [(dv[0] - u[0] * dot) / len, (dv[1] - u[1] * dot) / len, (dv[2] - u[2] * dot) / len];
}

/**
 * Derivative of the bond length r = |pos_a − pos_b|.
 * Returns dr/dx for each endpoint (same shape as distance()).
 *
 * Zero-length guard: r = 0 is a genuine cusp (the one-sided limit of
 * dr/dx is direction-dependent), so there is no correct finite value.
 * The zero return keeps the optimizer NaN-free; the energy term stays
 * finite-but-huge there, so a collapsed pair stalls at high energy
 * rather than exploding. Degenerate input only — never reached on any
 * validated geometry.
 */
export function bond_length_derivatives(
  pos_a: Vec3,
  pos_b: Vec3,
): { d_dx_a: Vec3; d_dx_b: Vec3 } {
  const ab = vec_sub(pos_a, pos_b);
  const r = vec_length(ab);
  if (r < 1e-15) return { d_dx_a: [0, 0, 0], d_dx_b: [0, 0, 0] };
  // dr/dx_a = (pos_a − pos_b)/r, dr/dx_b = −dr/dx_a
  const u = [ab[0] / r, ab[1] / r, ab[2] / r];
  return {
    d_dx_a: [u[0], u[1], u[2]],
    d_dx_b: [-u[0], -u[1], -u[2]],
  };
}

/**
 * Finite-difference fallback for the angle derivative at near-linear
 * geometries. The analytical form dθ/dx = −(1/sin θ)·dcos/dx is a 0/0
 * limit at exactly 0° and 180°, but the true derivative is finite —
 * the angle deficit is first-order in a transverse displacement.
 *
 * The angle is CUSPED in the transverse coordinate (θ = π − |δ|/r, a
 * cone like |x|), so a central difference at the exact cusp averages
 * the two one-sided slopes to zero. Forward differences of the angle
 * itself (δ = 1e-6 Å, the suite's FD convention) give the correct
 * directional derivative — and on each side of the cusp the angle is
 * smooth, so the forward difference is accurate to second order.
 */
function angle_derivatives_fd(
  pos_j: Vec3,
  pos_i: Vec3,
  pos_k: Vec3,
): { d_dx_i: Vec3; d_dx_j: Vec3; d_dx_k: Vec3 } {
  const DELTA = 1e-6;
  const theta_of = (pi: Vec3, pj: Vec3, pk: Vec3): number => {
    const v1 = vec_sub(pi, pj);
    const v2 = vec_sub(pk, pj);
    const l1 = vec_length(v1);
    const l2 = vec_length(v2);
    const u1: Vec3 = [v1[0] / l1, v1[1] / l1, v1[2] / l1];
    const u2: Vec3 = [v2[0] / l2, v2[1] / l2, v2[2] / l2];
    const c = Math.max(-1.0, Math.min(1.0, vec_dot(u1, u2)));
    return Math.acos(c);
  };
  const theta0 = theta_of(pos_i, pos_j, pos_k);
  const positions: [Vec3, Vec3, Vec3] = [pos_i, pos_j, pos_k];
  const d: [Vec3, Vec3, Vec3] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  const axes: [Vec3, Vec3, Vec3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let n = 0; n < 3; n++) {
    for (let a = 0; a < 3; a++) {
      const e = axes[a];
      const moved = positions.map((p, idx) =>
        idx === n
          ? [p[0] + DELTA * e[0], p[1] + DELTA * e[1], p[2] + DELTA * e[2]]
          : p,
      ) as [Vec3, Vec3, Vec3];
      d[n][a] = (theta_of(moved[0], moved[1], moved[2]) - theta0) / DELTA;
    }
  }
  return { d_dx_i: d[0], d_dx_j: d[1], d_dx_k: d[2] };
}

/**
 * Derivative of the bond angle θ (radians) at the central atom j,
 * for the angle i−j−k. Mirrors angle_in_radians(): unit vectors
 * (i−j) and (k−j), then θ = acos(û₁·û₂).
 *
 * dθ/dx = −(1/sin θ) · d(û₁·û₂)/dx
 */
export function angle_derivatives(
  pos_j: Vec3,
  pos_i: Vec3,
  pos_k: Vec3,
): { d_dx_i: Vec3; d_dx_j: Vec3; d_dx_k: Vec3 } {
  const v_ji = vec_sub(pos_i, pos_j);
  const v_jk = vec_sub(pos_k, pos_j);
  const len_ji = vec_length(v_ji);
  const len_jk = vec_length(v_jk);
  const u_ji: Vec3 = [v_ji[0] / len_ji, v_ji[1] / len_ji, v_ji[2] / len_ji];
  const u_jk: Vec3 = [v_jk[0] / len_jk, v_jk[1] / len_jk, v_jk[2] / len_jk];

  const cos_theta = vec_dot(u_ji, u_jk);
  const sin_theta = Math.sqrt(Math.max(0.0, 1.0 - cos_theta * cos_theta));
  if (sin_theta < 1e-12) {
    // Linear or degenerate angle: the 1/sin θ form is a 0/0 limit,
    // not a zero derivative — the true force bends the angle (the
    // energy terms use the cosine form for linear centers; non-linear
    // centers at exactly 0°/180° — e.g. a sketcher's linear water —
    // get the finite-difference limit here).
    return angle_derivatives_fd(pos_j, pos_i, pos_k);
  }

  // d(û₁·û₂)/dx_i = û₂ · dû₁/dx_i (û₂ fixed), and û₁ depends only on i and j.
  // d(û₁·û₂)/dx_k = û₁ · dû₂/dx_k. For j, both unit vectors move.
  const dcos_dx_i: Vec3 = [0, 0, 0];
  const dcos_dx_k: Vec3 = [0, 0, 0];
  const dcos_dx_j: Vec3 = [0, 0, 0];
  const axes: [Vec3, Vec3, Vec3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let a = 0; a < 3; a++) {
    const e = axes[a];
    const du_ji = unit_vec_deriv(v_ji, e);
    const du_jk = unit_vec_deriv(v_jk, e);
    dcos_dx_i[a] = vec_dot(du_ji, u_jk);
    dcos_dx_k[a] = vec_dot(u_ji, du_jk);
    dcos_dx_j[a] = -dcos_dx_i[a] - dcos_dx_k[a];
  }

  const scale = -1.0 / sin_theta;
  return {
    d_dx_i: [scale * dcos_dx_i[0], scale * dcos_dx_i[1], scale * dcos_dx_i[2]],
    d_dx_j: [scale * dcos_dx_j[0], scale * dcos_dx_j[1], scale * dcos_dx_j[2]],
    d_dx_k: [scale * dcos_dx_k[0], scale * dcos_dx_k[1], scale * dcos_dx_k[2]],
  };
}

/**
 * Derivative of the dihedral angle τ (radians) for i−j−k−l.
 * Mirrors dihedral_angle() exactly: v1 = i−j, v2 = k−j, v3 = l−k,
 * normals n1 = v1×v2 and n2 = v3×v2 (same handedness), then
 * τ = atan2(s, c) with s = (n̂₁×n̂₂)·v̂₂ and c = n̂₁·n̂₂.
 *
 * dτ = c·ds − s·dc, because s² + c² = 1 for unit vectors.
 */
export function dihedral_derivatives(
  pos_i: Vec3,
  pos_j: Vec3,
  pos_k: Vec3,
  pos_l: Vec3,
): { d_dx_i: Vec3; d_dx_j: Vec3; d_dx_k: Vec3; d_dx_l: Vec3 } {
  const v1 = vec_sub(pos_i, pos_j);
  const v2 = vec_sub(pos_k, pos_j);
  const v3 = vec_sub(pos_l, pos_k);

  const n1 = vec_cross(v1, v2);
  const n2 = vec_cross(v3, v2);

  const len1 = vec_length(n1);
  const len2 = vec_length(n2);
  const len_v2 = vec_length(v2);
  if (len1 < 1e-15 || len2 < 1e-15 || len_v2 < 1e-15) {
    return { d_dx_i: [0, 0, 0], d_dx_j: [0, 0, 0], d_dx_k: [0, 0, 0], d_dx_l: [0, 0, 0] };
  }

  const n1_hat: Vec3 = [n1[0] / len1, n1[1] / len1, n1[2] / len1];
  const n2_hat: Vec3 = [n2[0] / len2, n2[1] / len2, n2[2] / len2];
  const v2_hat: Vec3 = [v2[0] / len_v2, v2[1] / len_v2, v2[2] / len_v2];

  const c = vec_dot(n1_hat, n2_hat);
  const n1x2 = vec_cross(n1_hat, n2_hat);
  const s = vec_dot(n1x2, v2_hat);

  // How each vector moves when atom X is displaced along axis α:
  //   v1 = i − j  →  +δ on i, −δ on j
  //   v2 = k − j  →  +δ on k, −δ on j
  //   v3 = l − k  →  +δ on l, −δ on k
  const axes: [Vec3, Vec3, Vec3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

  function contribution(dv1: Vec3, dv2: Vec3, dv3: Vec3): number {
    // dn1 = dv1×v2 + v1×dv2, dn2 = dv3×v2 + v3×dv2
    const dn1 = vec_add3(vec_cross(dv1, v2), vec_cross(v1, dv2));
    const dn2 = vec_add3(vec_cross(dv3, v2), vec_cross(v3, dv2));
    const dn1_hat = unit_vec_deriv(n1, dn1);
    const dn2_hat = unit_vec_deriv(n2, dn2);
    const dv2_hat = unit_vec_deriv(v2, dv2);
    // dc = dn̂₁·n̂₂ + n̂₁·dn̂₂
    const dc = vec_dot(dn1_hat, n2_hat) + vec_dot(n1_hat, dn2_hat);
    // ds = (dn̂₁×n̂₂ + n̂₁×dn̂₂)·v̂₂ + (n̂₁×n̂₂)·dv̂₂
    const ds =
      vec_dot(vec_add3(vec_cross(dn1_hat, n2_hat), vec_cross(n1_hat, dn2_hat)), v2_hat) +
      vec_dot(n1x2, dv2_hat);
    // dτ = c·ds − s·dc
    return c * ds - s * dc;
  }

  const d_dx_i: Vec3 = [0, 0, 0];
  const d_dx_j: Vec3 = [0, 0, 0];
  const d_dx_k: Vec3 = [0, 0, 0];
  const d_dx_l: Vec3 = [0, 0, 0];
  for (let a = 0; a < 3; a++) {
    const e = axes[a];
    d_dx_i[a] = contribution(e, [0, 0, 0], [0, 0, 0]);
    d_dx_j[a] = contribution([-e[0], -e[1], -e[2]], [-e[0], -e[1], -e[2]], [0, 0, 0]);
    d_dx_k[a] = contribution([0, 0, 0], e, [-e[0], -e[1], -e[2]]);
    d_dx_l[a] = contribution([0, 0, 0], [0, 0, 0], e);
  }

  return { d_dx_i, d_dx_j, d_dx_k, d_dx_l };
}

/**
 * Finite-difference fallback for the Wilson out-of-plane angle
 * derivative at χ near ±90° (the substituent perpendicular to the
 * reference plane). The analytical form dχ/dx = (1/cos χ)·ds/dx is a
 * 0/0 limit there. Like the bond angle at collinearity, χ is cusped
 * in the transverse coordinate (χ = 90° − |δ|/r), so forward
 * differences (δ = 1e-6 Å) give the correct directional derivative —
 * the true force pushes the perpendicular substituent back toward the
 * plane.
 */
function oop_angle_derivatives_fd(
  pos_i: Vec3,
  pos_j: Vec3,
  pos_k: Vec3,
  pos_l: Vec3,
): { d_dx_i: Vec3; d_dx_j: Vec3; d_dx_k: Vec3; d_dx_l: Vec3 } {
  const DELTA = 1e-6;
  const chi_of = (pi: Vec3, pj: Vec3, pk: Vec3, pl: Vec3): number => {
    const unit = (p: Vec3, q: Vec3): Vec3 => {
      const v = vec_sub(p, q);
      const l = vec_length(v);
      return [v[0] / l, v[1] / l, v[2] / l];
    };
    const u_ji = unit(pi, pj);
    const u_jk = unit(pk, pj);
    const u_jl = unit(pl, pj);
    const n = vec_cross(u_ji, u_jk);
    const len = vec_length(n);
    if (len < 1e-15) return 0; // i and k collinear with j — χ = 0 everywhere
    const n_hat: Vec3 = [n[0] / len, n[1] / len, n[2] / len];
    const s = Math.max(-1.0, Math.min(1.0, vec_dot(n_hat, u_jl)));
    return Math.asin(s);
  };
  const chi0 = chi_of(pos_i, pos_j, pos_k, pos_l);
  const positions: [Vec3, Vec3, Vec3, Vec3] = [pos_i, pos_j, pos_k, pos_l];
  const d: [Vec3, Vec3, Vec3, Vec3] = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]];
  const axes: [Vec3, Vec3, Vec3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let n = 0; n < 4; n++) {
    for (let a = 0; a < 3; a++) {
      const e = axes[a];
      const moved = positions.map((p, idx) =>
        idx === n
          ? [p[0] + DELTA * e[0], p[1] + DELTA * e[1], p[2] + DELTA * e[2]]
          : p,
      ) as [Vec3, Vec3, Vec3, Vec3];
      d[n][a] =
        (chi_of(moved[0], moved[1], moved[2], moved[3]) - chi0) / DELTA;
    }
  }
  return { d_dx_i: d[0], d_dx_j: d[1], d_dx_k: d[2], d_dx_l: d[3] };
}

/**
 * Derivative of the Wilson out-of-plane angle χ (radians) for the
 * tri-coordinate center j. Mirrors wilson_oop_angle() exactly: unit
 * vectors (i−j), (k−j), (l−j); normal n = û_ji × û_jk; then
 * sin χ = n̂·û_jl, χ = asin(sin χ).
 *
 * dχ/dx = (1/√(1−s²)) · ds/dx
 */
export function oop_angle_derivatives(
  pos_i: Vec3,
  pos_j: Vec3,
  pos_k: Vec3,
  pos_l: Vec3,
): { d_dx_i: Vec3; d_dx_j: Vec3; d_dx_k: Vec3; d_dx_l: Vec3 } {
  const v_ji = vec_sub(pos_i, pos_j);
  const v_jk = vec_sub(pos_k, pos_j);
  const v_jl = vec_sub(pos_l, pos_j);

  const len_ji = vec_length(v_ji);
  const len_jk = vec_length(v_jk);
  const len_jl = vec_length(v_jl);
  if (len_ji < 1e-15 || len_jk < 1e-15 || len_jl < 1e-15) {
    return { d_dx_i: [0, 0, 0], d_dx_j: [0, 0, 0], d_dx_k: [0, 0, 0], d_dx_l: [0, 0, 0] };
  }

  const u_ji: Vec3 = [v_ji[0] / len_ji, v_ji[1] / len_ji, v_ji[2] / len_ji];
  const u_jk: Vec3 = [v_jk[0] / len_jk, v_jk[1] / len_jk, v_jk[2] / len_jk];
  const u_jl: Vec3 = [v_jl[0] / len_jl, v_jl[1] / len_jl, v_jl[2] / len_jl];

  const n = vec_cross(u_ji, u_jk);
  const len_n = vec_length(n);
  if (len_n < 1e-15) {
    // i and k collinear with j: every point is "in the plane", χ = 0
    // everywhere (the energy term returns 0 for this geometry too).
    return { d_dx_i: [0, 0, 0], d_dx_j: [0, 0, 0], d_dx_k: [0, 0, 0], d_dx_l: [0, 0, 0] };
  }
  const n_hat: Vec3 = [n[0] / len_n, n[1] / len_n, n[2] / len_n];

  const sin_chi = vec_dot(n_hat, u_jl);
  const cos_chi = Math.sqrt(Math.max(0.0, 1.0 - sin_chi * sin_chi));
  if (cos_chi < 1e-12) {
    // χ at ±90°: the substituent is perpendicular to the reference
    // plane — the 1/cos χ form is a 0/0 limit, not a zero derivative
    // (the true force pushes the substituent back toward the plane).
    return oop_angle_derivatives_fd(pos_i, pos_j, pos_k, pos_l);
  }

  // ds/dx = dn̂·û_jl + n̂·dû_jl, with n̂ = n/|n| built from unit vectors
  // u_ji × u_jk (mirroring wilson_oop_angle's normalization order).
  const axes: [Vec3, Vec3, Vec3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

  function contribution(du_ji: Vec3, du_jk: Vec3, du_jl: Vec3): number {
    const dn = vec_add3(vec_cross(du_ji, u_jk), vec_cross(u_ji, du_jk));
    const dn_hat = unit_vec_deriv(n, dn);
    return vec_dot(dn_hat, u_jl) + vec_dot(n_hat, du_jl);
  }

  const d_dx_i: Vec3 = [0, 0, 0];
  const d_dx_j: Vec3 = [0, 0, 0];
  const d_dx_k: Vec3 = [0, 0, 0];
  const d_dx_l: Vec3 = [0, 0, 0];
  for (let a = 0; a < 3; a++) {
    const e = axes[a];
    // v_ji = i − j: moves with i (+e) and j (−e); same for v_jk with k.
    const du_ji_i = unit_vec_deriv(v_ji, e);
    const du_ji_j = unit_vec_deriv(v_ji, [-e[0], -e[1], -e[2]]);
    const du_jk_k = unit_vec_deriv(v_jk, e);
    const du_jk_j = unit_vec_deriv(v_jk, [-e[0], -e[1], -e[2]]);
    const du_jl_l = unit_vec_deriv(v_jl, e);
    const du_jl_j = unit_vec_deriv(v_jl, [-e[0], -e[1], -e[2]]);
    d_dx_i[a] = contribution(du_ji_i, [0, 0, 0], [0, 0, 0]);
    d_dx_k[a] = contribution([0, 0, 0], du_jk_k, [0, 0, 0]);
    d_dx_l[a] = contribution([0, 0, 0], [0, 0, 0], du_jl_l);
    d_dx_j[a] = contribution(du_ji_j, du_jk_j, du_jl_j);
  }

  const scale = 1.0 / cos_chi;
  return {
    d_dx_i: [scale * d_dx_i[0], scale * d_dx_i[1], scale * d_dx_i[2]],
    d_dx_j: [scale * d_dx_j[0], scale * d_dx_j[1], scale * d_dx_j[2]],
    d_dx_k: [scale * d_dx_k[0], scale * d_dx_k[1], scale * d_dx_k[2]],
    d_dx_l: [scale * d_dx_l[0], scale * d_dx_l[1], scale * d_dx_l[2]],
  };
}

/** a + b + c — small local helper to keep the chain-rule lines readable. */
function vec_add3(a: Vec3, b: Vec3): Vec3 {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

/** Export for the gradient terms' unit conventions. */
export const RAD_PER_DEG = Math.PI / 180.0;
