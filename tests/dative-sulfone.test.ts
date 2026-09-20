import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import type { Molecule } from '../src/types.js';
import { parse_sdf } from '../src/sdf.js';
import { assign_atom_types } from '../src/mmff94/assign-atom-types.js';
import { assign_bci_charges } from '../src/mmff94/charges.js';
import { calc_energy } from '../src/mmff94/energy/total.js';
import { normalize_dative_sulfur } from '../src/mmff94/normalize-dative.js';

// One molecule, two drawings. dimethyl-sulfone.sdf is dimethyl sulfone as
// PubChem publishes it (CID 6213, explicit S=O double bonds);
// dimethyl-sulfone-dative.sdf is the same geometry with those bonds written
// as dative S(+2)–O(−) pairs — the form a SMILES round-trip or an editor can
// hand you. MMFF94 has types for the double-bonded form only, so typed
// exactly as drawn the dative sulfur lands in the thiol/sulfide class and
// the S=O parameters are never reached: 70.1491 kcal/mol against 4.8984.
// assign_atom_types() normalizes the dative form first
// (src/mmff94/normalize-dative.ts), and these tests pin the pair to one
// answer. OpenBabel normalizes it too — its reference logs for the two
// drawings are byte-identical.

function fixture(name: string): Molecule {
  return parse_sdf(readFileSync(`tests/fixtures/sdf/${name}.sdf`, 'utf-8'));
}

function typed_and_charged(name: string) {
  const typed = assign_atom_types(fixture(name));
  return { typed, charged: assign_bci_charges(typed) };
}

describe('dative S(+2)–O(−) drawings of a sulfone', () => {
  it('type identically to the hypervalent drawing', () => {
    const hyper = typed_and_charged('dimethyl-sulfone');
    const dative = typed_and_charged('dimethyl-sulfone-dative');
    expect(dative.typed.atom_types).toEqual(hyper.typed.atom_types);
    // the sulfone sulfur and both S=O oxygens, not the sulfide class 15
    expect(hyper.typed.atom_types.slice(0, 3)).toEqual([18, 32, 32]);
  });

  it('carry identical partial charges, and no net charge', () => {
    const hyper = typed_and_charged('dimethyl-sulfone');
    const dative = typed_and_charged('dimethyl-sulfone-dative');
    expect(dative.charged.partial_charges).toEqual(hyper.charged.partial_charges);
    const net = (q?: number[]) => (q ?? []).reduce((a, b) => a + b, 0);
    expect(net(hyper.charged.partial_charges)).toBeCloseTo(0, 9);
    expect(net(dative.charged.partial_charges)).toBeCloseTo(0, 9);
  });

  it('give the same energy', () => {
    const hyper = calc_energy(typed_and_charged('dimethyl-sulfone').charged);
    const dative = calc_energy(typed_and_charged('dimethyl-sulfone-dative').charged);
    expect(dative.total).toBeCloseTo(hyper.total, 9);
    expect(hyper.total).toBeCloseTo(4.8984, 3);
  });

  it('leaves a charge-separated molecule that MMFF94 does type alone', () => {
    // The Wittig ylide is P(+1)–C(−1) by design: MMFF94 has ylide types for
    // exactly that form, so the normalization must return it untouched.
    const ylide = fixture('methylenetriphenylphosphorane');
    expect(normalize_dative_sulfur(ylide)).toBe(ylide);
  });
});
