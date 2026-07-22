/**
 * Unit tests for 3D vector math.
 *
 * Tests every function in src/utils/vector.ts.
 * These are the primitives every energy term depends on,
 * so they must be rigorously correct.
 */

import { describe, it, expect } from 'vitest';
import {
  vec_add, vec_sub, vec_scale,
  vec_dot, vec_cross, vec_length, vec_normalize,
  distance, angle_in_radians, dihedral_angle, rotate_around_axis,
} from '../src/utils/vector';

describe('vector math', () => {
  it('vec_add adds two vectors elementwise', () => {
    expect(vec_add([1, 2, 3], [4, 5, 6])).toEqual([5, 7, 9]);
  });

  it('vec_sub subtracts two vectors elementwise', () => {
    expect(vec_sub([4, 5, 6], [1, 2, 3])).toEqual([3, 3, 3]);
  });

  it('vec_scale multiplies each component by scalar', () => {
    expect(vec_scale([1, 2, 3], 2)).toEqual([2, 4, 6]);
  });

  it('vec_dot returns correct dot product', () => {
    expect(vec_dot([1, 0, 0], [1, 0, 0])).toBeCloseTo(1);
    expect(vec_dot([1, 0, 0], [0, 1, 0])).toBeCloseTo(0);
  });

  it('vec_cross returns correct cross product', () => {
    const cross = vec_cross([1, 0, 0], [0, 1, 0]);
    expect(cross[0]).toBeCloseTo(0);
    expect(cross[1]).toBeCloseTo(0);
    expect(cross[2]).toBeCloseTo(1);
  });

  it('vec_length returns correct magnitude', () => {
    expect(vec_length([3, 4, 0])).toBeCloseTo(5);
  });

  it('vec_normalize returns unit vector', () => {
    const n = vec_normalize([3, 4, 0]);
    expect(vec_length(n)).toBeCloseTo(1);
    expect(n[0]).toBeCloseTo(0.6);
    expect(n[1]).toBeCloseTo(0.8);
  });

  it('vec_normalize of zero vector returns zero vector', () => {
    expect(vec_normalize([0, 0, 0])).toEqual([0, 0, 0]);
  });

  it('distance returns Euclidean distance', () => {
    expect(distance([0, 0, 0], [3, 4, 0])).toBeCloseTo(5);
  });

  it('angle_in_radians gives 90° for orthogonal vectors', () => {
    const angle = angle_in_radians([0, 0, 0], [1, 0, 0], [0, 1, 0]);
    expect(angle).toBeCloseTo(Math.PI / 2);
  });

  it('angle_in_radians gives 0 for parallel vectors', () => {
    const angle = angle_in_radians([0, 0, 0], [1, 0, 0], [2, 0, 0]);
    expect(angle).toBeCloseTo(0);
  });

  it('dihedral_angle gives 0° for cis (eclipsed) conformation', () => {
    const tau = dihedral_angle(
      [0, 0, 0],   // i
      [1, 0, 0],   // j
      [2, 0, 0],   // k
      [3, 0, 0],   // l
    );
    // All atoms colinear: dihedral is undefined (or 0 by convention)
    expect(Math.abs(tau)).toBeCloseTo(0, 1);
  });

  it('rotate_around_axis rotates a vector 90° about the z-axis', () => {
    const rotated = rotate_around_axis([1, 0, 0], [0, 0, 1], Math.PI / 2);
    expect(rotated[0]).toBeCloseTo(0, 5);
    expect(rotated[1]).toBeCloseTo(1, 5);
    expect(rotated[2]).toBeCloseTo(0, 5);
  });
});
