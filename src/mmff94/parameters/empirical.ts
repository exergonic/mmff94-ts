// Part V empirical-rule parameter generation (Halgren 1996, J. Comput.
// Chem. 17, 616-641): the designed fallback when a class-scoped lookup
// misses — the single home for the "derive parameters on the fly from
// the spec" machinery: the bond rules (eqs. 18-19), the angle rules
// (the θ₀ protocol + eq. 20), and the torsion rules (pp. 631-632, with
// the Table X U_i/V_i/W constants). The eq. 17 BCI fallback stays
// inside charges.ts — one term of the charge sum, inseparable from the
// q⁰/α sharing machinery.
//
//   eqs. (18)-(19) — the bond rules:
//     r0(ij) = r_i + r_j − c·|χ_i − χ_j|^1.4
//       r_i/r_j  covalent single-bond radii; c = 0.050 if either atom
//                is H, else 0.085
//       χ        Pauling-scale electronegativities
//     k_b = k_ref·(r_ref/r0)^6   (k_ref/r_ref: Table V)
//   the θ₀ protocol + eq. (20) — the angle rules:
//     θ₀ from the central atom's coordination/valency/ring class;
//     k_a = β·Z_i·C_j·Z_k / (θ₀²·(r̃_ij+r̃_jk)·exp(2D)) with
//       D = ((r̃_ij − r̃_jk)²/(r̃_ij + r̃_jk))² — BOTH numerator and
//           denominator squared (the paper's form, verified visually;
//           Tinker's kbond.f agrees; OpenBabel's GetAngleParam omits
//           the numerator square — a deviation),
//       β = 1.75 (×0.85 in 4-rings, ×0.05 in 3-rings), Z/C: Table VI.
//   the torsion rules (pp. 631-632) — the V2/V3 barrier heights from
//     the Table X U_i/V_i constants and the W values (2 for O, 8 for
//     S), by the central bond's aromaticity, the sp2/sp3 character of
//     j and k, and their mltb/pilp flags; rules (a)/(e)/(f)/(g-case-1)
//     SKIP the torsion (linear centers, unsaturated sp2 pairs, and
//     specific mltb/pilp combinations).
//
// The bond form here is MEASURED against the reference, not the
// paper's literal eq. (18): the paper adds a δ = 0.008 Å shrinkage and
// the mltb hybridization-index / BOij radius reductions, but the
// reference's only empirical bond (OHMW1's O–H, "Empirical rule bond
// parameters: 0 4 5" in the CCL MMFF94.empirical_rule_parameters file)
// matches the PLAIN form — r0 = 0.72 + 0.33 − 0.050·1.30^1.4 = 0.9778 —
// to 3.5e-6, and Tinker's own transcription (kbond.f) implements
// exactly this plain form (its bl(i) = rad0a + rad0b − cst·|χa−χb|^1.4,
// no δ, no corrections). The published tables are used as-is: r(O) =
// 0.72 (the CCL errata and Tinker's mmffcovrad agree).
//
// The angle rules never fire in the suite (the par's wildcards always
// resolve, validated 747/747); they are pinned against the spec by
// tests/angle-empirical.test.ts.
import type { Atom } from '../../types.js';
import type { AtomTypeProperties } from './atom-type-properties.js';

/** Covalent single-bond radii (Å) — the CCL eq. (18) list. */
export const COVALENT_RADII: Record<string, number> = {
  H: 0.33, Li: 1.34, Be: 0.90, B: 0.81, C: 0.77, N: 0.73, O: 0.72,
  F: 0.74, Na: 1.54, Mg: 1.30, Al: 1.22, Si: 1.15, P: 1.09, S: 1.03,
  Cl: 1.01, K: 1.96, Ca: 1.74, Sc: 1.44, Ti: 1.36, Cu: 1.38, Zn: 1.31,
  Ga: 1.19, Ge: 1.20, As: 1.20, Se: 1.16, Br: 1.15, Rb: 2.11, Sr: 1.92,
  Y: 1.62, Zr: 1.48, Ag: 1.53, Cd: 1.48, In: 1.46, Sn: 1.40, Sb: 1.41,
  Te: 1.35, I: 1.33,
};

/** Pauling-scale (Allred-Rochow) electronegativities — the CCL eq. (16) list. */
export const ELECTRONEGATIVITY: Record<string, number> = {
  H: 2.2, Li: 0.97, Be: 1.47, B: 2.01, C: 2.5, N: 3.07, O: 3.5,
  F: 4.1, Na: 1.01, Mg: 1.23, Al: 1.47, Si: 1.74, P: 2.06, S: 2.44,
  Cl: 2.83, K: 0.91, Ca: 1.04, Sc: 1.3, Ti: 1.5, V: 1.6, Cr: 1.6,
  Mn: 1.5, Fe: 1.8, Co: 1.8, Ni: 1.8, Cu: 1.9, Zn: 1.6, Ga: 1.82,
  Ge: 2.02, As: 2.2, Se: 2.48, Br: 2.74, Rb: 0.89, Sr: 0.99, Y: 1.3,
  Zr: 1.4, Nb: 1.6, Mo: 1.8, Tc: 1.9, Ru: 2.2, Rh: 2.2, Pd: 2.2,
  Ag: 1.9, Cd: 1.7, In: 1.49, Sn: 1.72, Sb: 1.82, Te: 2.01, I: 2.21,
};

/**
 * Table V — reference r0/k_b for the empirical force-constant rule,
 * keyed by the element pair (atomic numbers, ascending). Mechanical
 * transcription of the paper's table (identical to mmffbndk.par).
 */
export const BNDK_REF: Record<string, { r0_ref: number; k_ref: number }> = {
  '1-6': { r0_ref: 1.084, k_ref: 5.15 }, '1-7': { r0_ref: 1.001, k_ref: 7.35 },
  '1-8': { r0_ref: 0.947, k_ref: 9.1 }, '1-9': { r0_ref: 0.92, k_ref: 10.6 },
  '1-14': { r0_ref: 1.48, k_ref: 2.3 }, '1-15': { r0_ref: 1.415, k_ref: 2.95 },
  '1-16': { r0_ref: 1.326, k_ref: 4.3 }, '1-17': { r0_ref: 1.28, k_ref: 4.3 },
  '1-35': { r0_ref: 1.41, k_ref: 4.2 }, '1-53': { r0_ref: 1.6, k_ref: 2.7 },
  '6-6': { r0_ref: 1.512, k_ref: 3.8 }, '6-7': { r0_ref: 1.439, k_ref: 4.55 },
  '6-8': { r0_ref: 1.393, k_ref: 5.4 }, '6-9': { r0_ref: 1.353, k_ref: 6.2 },
  '6-14': { r0_ref: 1.86, k_ref: 2.6 }, '6-15': { r0_ref: 1.84, k_ref: 2.7 },
  '6-16': { r0_ref: 1.812, k_ref: 2.85 }, '6-17': { r0_ref: 1.781, k_ref: 2.75 },
  '6-35': { r0_ref: 1.94, k_ref: 2.6 }, '6-53': { r0_ref: 2.16, k_ref: 1.4 },
  '7-7': { r0_ref: 1.283, k_ref: 6.0 }, '7-8': { r0_ref: 1.333, k_ref: 5.9 },
  '7-9': { r0_ref: 1.36, k_ref: 5.9 }, '7-14': { r0_ref: 1.74, k_ref: 3.7 },
  '7-15': { r0_ref: 1.65, k_ref: 4.8 }, '7-16': { r0_ref: 1.674, k_ref: 3.5 },
  '7-17': { r0_ref: 1.75, k_ref: 3.5 }, '7-35': { r0_ref: 1.9, k_ref: 2.9 },
  '7-53': { r0_ref: 2.1, k_ref: 1.6 },
  '8-8': { r0_ref: 1.48, k_ref: 3.6 }, '8-9': { r0_ref: 1.42, k_ref: 4.6 },
  '8-14': { r0_ref: 1.63, k_ref: 5.2 }, '8-15': { r0_ref: 1.66, k_ref: 4.7 },
  '8-16': { r0_ref: 1.47, k_ref: 9.9 }, '8-17': { r0_ref: 1.7, k_ref: 4.1 },
  '8-35': { r0_ref: 1.85, k_ref: 3.4 }, '8-53': { r0_ref: 2.05, k_ref: 1.6 },
  '9-14': { r0_ref: 1.57, k_ref: 6.4 }, '9-15': { r0_ref: 1.54, k_ref: 7.1 },
  '9-16': { r0_ref: 1.55, k_ref: 6.9 },
  '14-14': { r0_ref: 2.32, k_ref: 1.3 }, '14-15': { r0_ref: 2.25, k_ref: 1.5 },
  '14-16': { r0_ref: 2.15, k_ref: 2.0 }, '14-17': { r0_ref: 2.02, k_ref: 3.1 },
  '14-35': { r0_ref: 2.19, k_ref: 2.1 }, '14-53': { r0_ref: 2.44, k_ref: 1.5 },
  '15-15': { r0_ref: 2.21, k_ref: 1.7 }, '15-16': { r0_ref: 2.1, k_ref: 2.4 },
  '15-17': { r0_ref: 2.03, k_ref: 3.0 }, '15-35': { r0_ref: 2.21, k_ref: 2.0 },
  '15-53': { r0_ref: 2.47, k_ref: 1.4 },
  '16-16': { r0_ref: 2.052, k_ref: 2.5 }, '16-17': { r0_ref: 2.04, k_ref: 2.9 },
  '16-35': { r0_ref: 2.24, k_ref: 1.9 }, '16-53': { r0_ref: 2.4, k_ref: 1.7 },
  '17-17': { r0_ref: 1.99, k_ref: 3.5 },
  '35-35': { r0_ref: 2.28, k_ref: 2.4 },
  '53-53': { r0_ref: 2.67, k_ref: 1.6 },
};

/** Element → atomic number (the Table V keys and the θ₀ Z > 10 rule). */
export const ELEMENT_ATOMIC_NUMBER: Record<string, number> = {
  H: 1, Li: 3, Be: 4, B: 5, C: 6, N: 7, O: 8, F: 9, Na: 11, Mg: 12,
  Al: 13, Si: 14, P: 15, S: 16, Cl: 17, K: 19, Ca: 20, Sc: 21, Ti: 22,
  V: 23, Cr: 24, Mn: 25, Fe: 26, Co: 27, Ni: 28, Cu: 29, Zn: 30, Ga: 31,
  Ge: 32, As: 33, Se: 34, Br: 35, Rb: 37, Sr: 38, Y: 39, Zr: 40, Nb: 41,
  Mo: 42, Tc: 43, Ru: 44, Rh: 45, Pd: 46, Ag: 47, Cd: 48, In: 49, Sn: 50,
  Sb: 51, Te: 52, I: 53,
};

// MMFF part V, table VI: Z_i — the eq. (20) atomic parameter for the
// peripheral atoms (NOT the atomic number; the paper's footnote: the
// current version does not use B and As, so they are absent here, and
// Si's Z is assumed equal to phosphorus's).
const ELEMENT_Z: Record<string, number> = {
  H: 1.395, C: 2.494, N: 2.711, O: 3.045, F: 2.847,
  Si: 2.350, P: 2.350, S: 2.980, Cl: 2.909, Br: 3.017, I: 3.086,
};

// MMFF part V, table VI: C_i — the central atom's constant. No H/F/Br/I
// entries (0 by the ?? fallback) — an H-centered empirical angle gets
// k_a = 0, the spec.
const ELEMENT_C: Record<string, number> = {
  B: 0.704, C: 1.016, N: 1.113, O: 1.337,
  Si: 0.811, P: 1.068, S: 1.249, Cl: 1.078, As: 0.825,
};

/** Element → periodic-table row (the empirical torsion rules): 0 = Z
 *  ≤ 2, 1 = 3–10, 2 = 11–18, 3 = 19–36, 4 = 37–54. */
export const ELEMENT_ROW: Record<string, number> = {
  H: 0, B: 1, C: 1, N: 1, O: 1, F: 1,
  Si: 2, P: 2, S: 2, Cl: 2, Br: 3, I: 4,
};

// MMFF part V, table X: U_i — the torsion V2 driving constant.
const ELEMENT_U: Record<string, number> = {
  C: 2.0, N: 2.0, O: 2.0, Si: 1.25, P: 1.25, S: 1.25,
};
// MMFF part V, table X: V_i — the torsion V3 driving constant.
const ELEMENT_V: Record<string, number> = {
  C: 2.12, N: 1.5, O: 0.2, Si: 1.22, P: 2.4, S: 0.49,
};

/**
 * Eq. (18): the empirical reference bond length for i-j — the plain
 * Schomaker-Stevenson/Blom-Haaland form (no δ, no corrections; see the
 * header note). Returns undefined when the element data is missing
 * (radii or χ unknown).
 */
export function empirical_bond_length(a1: Atom, a2: Atom): number | undefined {
  const r0a = COVALENT_RADII[a1.element];
  const r0b = COVALENT_RADII[a2.element];
  const Xa = ELECTRONEGATIVITY[a1.element];
  const Xb = ELECTRONEGATIVITY[a2.element];
  if (r0a === undefined || r0b === undefined || Xa === undefined || Xb === undefined) {
    return undefined;
  }

  const c = a1.element === 'H' || a2.element === 'H' ? 0.05 : 0.085;
  return r0a + r0b - c * Math.pow(Math.abs(Xa - Xb), 1.4);
}

/**
 * Eqs. (18)-(19): the empirical bond parameters for i-j — the k_b from
 * Table V scaled by the inverse sixth power of the length ratio.
 * Undefined when the element pair is outside Table V (the paper's
 * Badger's-rule fallback is not implemented — the suite never needs
 * it, and its constants come from a cited external paper).
 */
export function empirical_bond_parameters(
  a1: Atom,
  a2: Atom,
): { k_b: number; r0: number } | undefined {
  const z1 = ELEMENT_ATOMIC_NUMBER[a1.element];
  const z2 = ELEMENT_ATOMIC_NUMBER[a2.element];
  if (z1 === undefined || z2 === undefined) return undefined;
  const key = z1 < z2 ? `${z1}-${z2}` : `${z2}-${z1}`;
  const ref = BNDK_REF[key];
  if (!ref) return undefined;

  const r0 = empirical_bond_length(a1, a2);
  if (r0 === undefined) return undefined;

  const rr = ref.r0_ref / r0;
  const rr6 = rr * rr * rr * rr * rr * rr;
  return { k_b: ref.k_ref * rr6, r0 };
}

/**
 * The θ₀ default protocol (part V, p. 626): the reference angle when
 * the lookup chain misses entirely — 120° initial, then the central
 * atom's coordination number / valency / linearity / ring class.
 */
export function empirical_theta0(
  prop: Pick<AtomTypeProperties, 'crd' | 'val' | 'mltb' | 'lin'> | undefined,
  element: string,
  cls: number,
): number {
  let theta0 = 120.0;
  if (prop) {
    if (prop.crd === 4) theta0 = 109.45;
    // The 105°/95°/180° branches live INSIDE the crd = 2 case (the
    // paper's structure) — an ungated atomic-number check would
    // override the tricoordinate 92°/107° rule for heavy centers.
    if (prop.crd === 2) {
      if (element === 'O') theta0 = 105.0;
      if ((ELEMENT_ATOMIC_NUMBER[element] ?? 0) > 10) theta0 = 95.0;
      if (prop.lin) theta0 = 180.0;
    }
    if (prop.crd === 3 && prop.val === 3 && !prop.mltb) {
      theta0 = element === 'N' ? 107.0 : 92.0;
    }
  }
  // Small-ring angles are forced to the ring geometry.
  if (cls === 3 || cls === 5 || cls === 6) theta0 = 60.0;
  if (cls === 4 || cls === 7 || cls === 8) theta0 = 90.0;
  return theta0;
}

/**
 * Eq. (20): the angle-bending force constant when the lookup misses —
 * β·Z_i·C_j·Z_k / (θ₀²·(r̃_ij+r̃_jk)·exp(2D)). The r̃'s are the
 * REFERENCE bond lengths (the par's r0, or the eq. (18) value for a
 * par-less bond — the caller resolves that chain); D is the squared
 * ratio (both numerator and denominator squared).
 */
export function empirical_ka(
  i: Atom,
  j: Atom,
  k: Atom,
  r0ab: number,
  r0bc: number,
  theta0: number,
  cls: number,
): number {
  const Za = ELEMENT_Z[i.element] ?? 0;
  const Cb = ELEMENT_C[j.element] ?? 0; // 0 for H — spec
  const Zc = ELEMENT_Z[k.element] ?? 0;

  // Ring strain scales the constant down: ×0.85 for 4-rings, ×0.05
  // for 3-rings.
  let beta = 1.75;
  if (cls === 4 || cls === 7 || cls === 8) beta *= 0.85;
  if (cls === 3 || cls === 5 || cls === 6) beta *= 0.05;

  const rr = r0ab + r0bc;
  const D = ((r0ab - r0bc) * (r0ab - r0bc)) / (rr * rr);
  const rad2 = (Math.PI / 180) * (Math.PI / 180);
  return (beta * Za * Cb * Zc * Math.exp(-2 * D)) / (rr * theta0 * theta0 * rad2);
}

/**
 * The empirical torsion barrier heights for i-j-k-l (part V,
 * pp. 631-632 — the extension of the part IV torsion parameterization;
 * the U_i/V_i constants are Table X), or skip = true when the rules
 * say this torsion does not exist (linear centers, and the rule
 * (e)/(f)/(g) exclusions). Rules (b)-(h) always set something. The
 * graph queries — the j-k bond order and aromaticity — are gathered
 * by the caller; the properties and elements of j and k arrive as
 * arguments.
 */
export function empirical_torsion(
  pj: Pick<AtomTypeProperties, 'lin' | 'pilp' | 'val' | 'mltb' | 'crd'> | undefined,
  pk: Pick<AtomTypeProperties, 'lin' | 'pilp' | 'val' | 'mltb' | 'crd'> | undefined,
  ej: string,
  ek: string,
  order_jk: number,
  aromatic_jk: boolean,
): { v1: number; v2: number; v3: number; skip: boolean } {
  const ub = ELEMENT_U[ej] ?? 0;
  const uc = ELEMENT_U[ek] ?? 0;
  const vb = ELEMENT_V[ej] ?? 0;
  const vc = ELEMENT_V[ek] ?? 0;
  const v = { v1: 0, v2: 0, v3: 0, skip: false };
  let found = false;

  // rule (a): linear centers carry no torsion
  if (pj?.lin || pk?.lin) { v.skip = true; return v; }

  // rules (b)/(c): V2 from the U parameters
  if (aromatic_jk) {
    // rule (b): aromatic central bond — π = 0.5 without lone pairs
    // on either atom, 0.3 with; β = 3 for the val 3/4 combination,
    // 6 otherwise
    const pi_bc = (!pj?.pilp && !pk?.pilp) ? 0.5 : 0.3;
    const beta =
      (pj?.val === 3 && pk?.val === 4) || (pj?.val === 4 && pk?.val === 3) ? 3.0 : 6.0;
    v.v2 = beta * pi_bc * Math.sqrt(ub * uc);
    found = true;
  } else if (order_jk === 2) {
    // rule (c): the j-k bond has a formal bond order of 2 — π = 1.0
    // only for the full double bond (mltb 2 on both), else 0.4 (e.g.
    // the formal C=N of a guanidinium resonance structure, mltb 2/1).
    // The paper's condition is on the j-k bond order; OpenBabel's
    // GetTorsionParam checks the i-j bond instead and Tinker's ktors
    // checks only the mltb flags — both readings give π = 1.0 for the
    // central single of butadiene, which the paper's rule (g) case
    // (5) assigns π = 0.15.
    const pi_bc = pj?.mltb === 2 && pk?.mltb === 2 ? 1.0 : 0.4;
    v.v2 = 6.0 * pi_bc * Math.sqrt(ub * uc);
    found = true;
  }

  // rule (d): both sp3 → V3
  if (!found && pj?.crd === 4 && pk?.crd === 4) {
    v.v3 = Math.sqrt(vb * vc) / 9.0;
    found = true;
  }

  // rules (e)/(f): sp3/sp2 mixed — no torsion for the unsaturated
  // sp2 combinations listed (a saturated sp3-sp2 bond falls through
  // to the remaining rules)
  const unsaturated_sp2 = (
    p: Pick<AtomTypeProperties, 'crd' | 'val' | 'mltb'> | undefined,
  ): boolean => {
    if (!p) return false;
    if (p.crd === 3) return p.val === 4 || p.val === 34 || p.mltb !== 0;
    if (p.crd === 2) return p.val === 3 || p.mltb !== 0;
    return false;
  };
  if (!found && pj?.crd === 4 && pk?.crd !== 4) {
    if (unsaturated_sp2(pk)) { v.skip = true; return v; }
  } else if (!found && pk?.crd === 4 && pj?.crd !== 4) {
    if (unsaturated_sp2(pj)) { v.skip = true; return v; }
  }

  // rule (g): order-1 central bond between mltb/pilp-carrying types
  // → V2 (the π value depends on which side carries what)
  if (!found) {
    const central_single = order_jk === 1 && !aromatic_jk;
    if (
      central_single &&
      ((pj?.mltb && pk?.mltb) || (pj?.mltb && pk?.pilp) || (pk?.mltb && pj?.pilp))
    ) {
      if (pj?.pilp && pk?.pilp) { v.skip = true; return v; } // case (1)
      let pi_bc = 0.15;
      if (pj?.pilp && pk?.mltb) { // case (2)
        // π = 0.5 when the PILP atom (j) is itself a strongly
        // delocalized single-bond former (mltb(J) = 1 — e.g. the
        // amide N type 10) — the paper and Tinker's ktors check j's
        // own mltb; OpenBabel checks k's, which is the wrong atom.
        if (pj.mltb === 1) pi_bc = 0.5;
        else if (ELEMENT_ROW[ej] === 1 && ELEMENT_ROW[ek] === 1) pi_bc = 0.3;
        else pi_bc = 0.15;
        found = true;
      }
      if (pk?.pilp && pj?.mltb) { // case (3)
        // Indices interchanged from case (2): now k carries the lone
        // pair, so the 0.5 test reads k's own mltb.
        if (pk.mltb === 1) pi_bc = 0.5;
        else if (ELEMENT_ROW[ej] === 1 && ELEMENT_ROW[ek] === 1) pi_bc = 0.3;
        else pi_bc = 0.15;
        found = true;
      }
      if (!found && (pj?.mltb === 1 || pk?.mltb === 1) && (ej !== 'C' || ek !== 'C')) {
        pi_bc = 0.4;
        found = true;
      }
      if (!found) pi_bc = 0.15;
      v.v2 = 6.0 * pi_bc * Math.sqrt(ub * uc);
      found = true;
    }
  }

  // rule (h): O/S central pair → negative V2 (W = 2 for O, 8 for S);
  // otherwise V3 from the V parameters over N_bc — the number of
  // torsion interactions about the j-k bond, (crd(J)−1)·(crd(K)−1)
  // (eq. (22), rule (d)'s definition; Tinker's ktors agrees, OpenBabel
  // uses crd·crd). The crd ?? 2 fallback keeps N_bc ≥ 1 for unknown
  // types — a torsion needs j and k to each carry a substituent.
  if (!found) {
    const o_s = (e: string) => e === 'O' || e === 'S';
    if (o_s(ej) && o_s(ek)) {
      v.v2 = -Math.sqrt((ej === 'O' ? 2.0 : 8.0) * (ek === 'O' ? 2.0 : 8.0));
    } else {
      v.v3 = Math.sqrt(vb * vc) / (((pj?.crd ?? 2) - 1) * ((pk?.crd ?? 2) - 1));
    }
  }
  return v;
}
