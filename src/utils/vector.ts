/**
 * 3D vector math — pure functions, no classes.
 *
 * These are the geometric primitives used by every energy term
 * and gradient function. They operate on plain [x, y, z] tuples.
 * No dependencies, no objects, no surprises.
 */

export type Vec3 = [number, number, number];

export function vec_add(a: Vec3, b: Vec3): Vec3 {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

export function vec_sub(a: Vec3, b: Vec3): Vec3 {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

export function vec_scale(v: Vec3, s: number): Vec3 {
  return [v[0] * s, v[1] * s, v[2] * s];
}

export function vec_dot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export function vec_cross(a: Vec3, b: Vec3): Vec3 {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

export function vec_length(v: Vec3): number {
  return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

export function vec_normalize(v: Vec3): Vec3 {
  const len = vec_length(v);
  if (len < 1e-15) return [0, 0, 0];
  return [v[0] / len, v[1] / len, v[2] / len];
}

/**
 * Distance between two points.
 */
export function distance(a: Vec3, b: Vec3): number {
  return vec_length(vec_sub(a, b));
}

/**
 * Angle in radians between vectors a→b and c→b (with common origin at b).
 */
export function angle_in_radians(b: Vec3, a: Vec3, c: Vec3): number {
  const ba = vec_normalize(vec_sub(a, b));
  const bc = vec_normalize(vec_sub(c, b));
  const dot = vec_dot(ba, bc);
  // Clamp to [-1, 1] to avoid NaN from floating-point rounding
  return Math.acos(Math.max(-1, Math.min(1, dot)));
}

/**
 * Dihedral (torsion) angle in radians for four consecutive atoms i−j−k−l.
 *
 * Returns the angle between the plane (i, j, k) and the plane (j, k, l).
 * Sign follows the right-hand rule about the j→k axis.
 * τ = 0 when i−j and k−l are eclipsed (cis), τ = π when staggered (trans).
 */
export function dihedral_angle(
  i: Vec3, j: Vec3, k: Vec3, l: Vec3
): number {
  const v1 = vec_sub(i, j);
  const v2 = vec_sub(k, j);
  const v3 = vec_sub(l, k);

  const n1 = vec_cross(v1, v2);
  const n2 = vec_cross(v2, v3);

  const n1_norm = vec_normalize(n1);
  const n2_norm = vec_normalize(n2);

  const cos_angle = vec_dot(n1_norm, n2_norm);
  const sin_angle = vec_dot(vec_cross(n1_norm, n2_norm), vec_normalize(v2));

  return Math.atan2(sin_angle, cos_angle);
}

/**
 * Rotate a vector v about an axis by angle theta (radians).
 * Uses Rodrigues' rotation formula.
 */
export function rotate_around_axis(
  v: Vec3, axis: Vec3, theta: number
): Vec3 {
  const a = vec_normalize(axis);
  const cos_t = Math.cos(theta);
  const sin_t = Math.sin(theta);
  const dot = vec_dot(v, a);

  const term1 = vec_scale(v, cos_t);
  const term2 = vec_scale(vec_cross(a, v), sin_t);
  const term3 = vec_scale(a, dot * (1 - cos_t));

  return vec_add(vec_add(term1, term2), term3);
}
