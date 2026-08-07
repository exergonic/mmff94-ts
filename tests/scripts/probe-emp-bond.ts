/**
 * Probe: our empirical bond generation for element pairs.
 * Run: npx tsx tests/scripts/probe-emp-bond.ts
 */
import { empirical_bond_length, empirical_bond_parameters } from '../../src/mmff94/parameters/empirical.js';

const mk = (element: string) => ({ element } as any);
for (const [a, b] of [['P', 'Si'], ['O', 'H'], ['S', 'H'], ['P', 'S']] as const) {
  const r = empirical_bond_length(mk(a), mk(b));
  const p = empirical_bond_parameters(mk(a), mk(b));
  console.log(`${a}-${b}: r0 ${r.toFixed(6)}  k_b ${p?.k_b ?? 'none'}`);
}
