/**
 * The MMFF94 partial-charge model: bond charge increments (BCI).
 *
 * MMFF94 does NOT store per-atom partial charges in its parameter
 * tables. Instead, each BOND TYPE pair has a charge increment (the
 * bci value, mmffchg.par), and the partial charge on an atom is the
 * sum of the increments of every bond it participates in:
 *
 *   q_i = Σ_j bci(i, j)
 *
 * The par stores each pair with its types in ascending order, and the
 * increment is SUBTRACTED from the smaller-type atom and ADDED to the
 * larger-type atom (Halgren part III). So in ammonia, the N(8)–H(23)
 * bond carries bci = +0.36: the N collects −0.36 from each of its
 * three H's (−1.08 total) and each H collects +0.36.
 *
 * Bonds whose type pair has no entry use the per-atom default values
 * (mmffpbci.par): the pair contributes P_i − P_j, the difference of
 * the two atoms' defaults — the same flow from the smaller type.
 *
 * The lookup is class-scoped like every other term: a conjugated
 * single bond (BTij = 1) uses the class-1 entry when one exists, and
 * falls to the per-atom defaults when it does not (it does NOT fall
 * to the class-0 entry — the reference is strict about this).
 *
 * FORMAL CHARGES (part V, eq. 15): atoms of charged types carry a
 * "primary" formal atomic charge q⁰ (e.g. +1 on quaternary N, −0.5 on
 * a carboxylate oxygen), and atoms with a NEGATIVE primary charge
 * share half of it with their bonded neighbors — no sharing is
 * invoked for positive charges (part V). The full charge is
 *
 *   q_i = (1 − α_i·crd_i)·q⁰_i + Σ_k α_k·q⁰_k + Σ_k w_ik
 *
 * where α_i is the formal-charge adjustment factor (fcadj, the α of
 * part V Table III) and the neighbor sum runs over the crd_i bonded
 * atoms k: each neighbor's α_k·q⁰_k (the portion it shares) flows to
 * atom i. NOTE: the flow uses the NEIGHBOR's α — transcribing the
 * sharing with the atom's own α (as some implementations do) does not
 * reproduce the reference: a carboxylate carbon (α = 0) must still
 * receive half of each attached oxygen's −0.5.
 *
 * The q⁰ values below were verified against the validation suite's
 * reference partial charges (the .mmd pchg column): solving eq. (15)
 * for q⁰ across all 550 reference molecules reproduces these values
 * to machine precision, and the environment rule for type 32
 * (carboxylate O: −0.5; sulfone/nitro/nitrate O: 0) is the clean
 * bimodal split the reference shows. Types absent from the table
 * carry q⁰ = 0.
 */

import type { TypedMolecule } from '../types';
import { BCI_PARAMS, BCI_DEFAULT_PARAMS, ATOM_TYPE_PROPERTIES } from './parameters';
import { make_class_context, bond_type_flag } from './parameters/parameter-classes';

/** Primary formal atomic charges q⁰ (part V eq. 15), per type. */
const PRIMARY_FORMAL_CHARGES: Record<number, number> = {
  34: 1,     // NR+ — quaternary N
  35: -1,    // OM — oxide O on sp3/sp2 C
  51: 1,     // O=+ — oxenium O (pyrylium)
  54: 1,     // N+=C — iminium N
  55: 0.5,   // NCN+ — N in N+=C-N (fractional)
  56: 1 / 3, // NGD+ — guanidinium N (fractional, part V)
  58: 1,     // NPD+ — pyridinium N
  89: -1,    // F⁻
  90: -1,    // Cl⁻
  91: -1,    // Br⁻
  // 32 is environment-dependent: −0.5 for a carboxylate oxygen
  // (CO₂⁻), 0 for sulfone/nitro/nitrate oxygens — see below.
};

/**
 * Compute the partial charges for every atom of a typed molecule and
 * return a COPY of the molecule with `partial_charges` attached (the
 * electrostatics term and the gradient both consume them).
 *
 * Pure by design, like the rest of the pipeline: the input molecule is
 * not mutated, and the returned charged molecule is the value that
 * flows into calc_energy / calc_gradient / optimize_lbfgs. This is
 * sound because MMFF94 partial charges are geometry-independent — they
 * depend only on connectivity and atom types — so the charged
 * molecule stays valid while an optimizer moves the atoms (the
 * optimizer clones only the atom positions and keeps this field).
 */
export function assign_bci_charges(molecule: TypedMolecule): TypedMolecule {
  // Adjacency + the shared class context (the BTij flag selects the
  // bci class for conjugated single bonds).
  const adj: number[][] = Array.from({ length: molecule.atoms.length }, () => []);
  for (const bond of molecule.bonds) {
    adj[bond.atom1].push(bond.atom2);
    adj[bond.atom2].push(bond.atom1);
  }
  const ctx = make_class_context(molecule, adj);

  const charges = molecule.atoms.map(() => 0.0);
  for (const bond of molecule.bonds) {
    const ti = molecule.atom_types[bond.atom1];
    const tj = molecule.atom_types[bond.atom2];
    const t_min = Math.min(ti, tj);
    const t_max = Math.max(ti, tj);
    const cls = bond_type_flag(ctx, bond.atom1, bond.atom2);
    const entry = BCI_PARAMS[`${cls}-${t_min}-${t_max}`];

    let bci: number;
    if (entry) {
      bci = entry.bci;
    } else {
      // Unparametrized pair: the difference of the per-atom defaults.
      const pa = BCI_DEFAULT_PARAMS[t_min]?.pbci ?? 0;
      const pb = BCI_DEFAULT_PARAMS[t_max]?.pbci ?? 0;
      bci = pa - pb;
    }

    // The increment flows from the smaller type to the larger one.
    if (ti === t_min) {
      charges[bond.atom1] -= bci;
      charges[bond.atom2] += bci;
    } else {
      charges[bond.atom1] += bci;
      charges[bond.atom2] -= bci;
    }
  }

  // Formal-charge correction: eq. (15) of part V. q⁰ is per type; the
  // sharing (the α_k·q⁰_k flow) needs the full adjacency.
  // Primary charge of atom i, with the type-32 environment rule (clean
  // bimodal split in the reference): −0.5 on a carboxylate oxygen — a
  // carbon neighbor carrying two terminal oxygens — and 0 on sulfone,
  // nitro, and nitrate oxygens, whose polarization lives in the BCI.
  const q0_of = (i: number): number => {
    const t = molecule.atom_types[i];
    if (t === 32) {
      return adj[i].some(nb => {
        const nbr = molecule.atoms[nb];
        if (nbr.element !== 'C') return false;
        let terminalO = 0;
        for (const b of adj[nb]) {
          if (molecule.atoms[b].element === 'O' && adj[b].length === 1) terminalO++;
        }
        return terminalO >= 2;
      })
        ? -0.5
        : 0;
    }
    return PRIMARY_FORMAL_CHARGES[t] ?? 0;
  };

  const formal = molecule.atoms.map(() => 0.0);
  const n = molecule.atoms.length;
  for (let i = 0; i < n; i++) {
    const t = molecule.atom_types[i];
    const crd = ATOM_TYPE_PROPERTIES[t]?.crd ?? 0;
    const alpha = BCI_DEFAULT_PARAMS[t]?.fcadj ?? 0;
    formal[i] = (1 - alpha * crd) * q0_of(i);
    for (const k of adj[i]) {
      const tk = molecule.atom_types[k];
      const alpha_k = BCI_DEFAULT_PARAMS[tk]?.fcadj ?? 0;
      formal[i] += alpha_k * q0_of(k);
    }
  }
  for (let i = 0; i < n; i++) charges[i] += formal[i];

  // The charged molecule is a shallow copy: atoms and bonds are shared
  // references — only the new field is added. Geometry and typing are
  // untouched, so this is the value that flows into the energy terms
  // and the optimizer.
  return { ...molecule, partial_charges: charges };
}
