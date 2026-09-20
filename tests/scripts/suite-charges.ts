/**
 * The per-atom reference partial charges from the validation suite's .mmd
 * `pchg` column — the shared extractor behind the charges gate
 * (tests/charges-suite.test.ts) and the validation-report generator
 * (tests/scripts/generate-validation-doc.ts).
 *
 * This is NOT the molecule parser's `atoms[].partial_charge`: that column
 * is empty across the suite, so the generator's earlier comparison ran on
 * undefined values, every deviation came out NaN, the NaN never advanced
 * the `worst` accumulator, and the report printed "Worst 0.00e+0 ()" with
 * a 100% pass rate — a missing value wearing a lab coat (fixed 2026-09-20).
 */
export function reference_charges(mmdText: string, name: string, nAtoms: number): number[] {
  const pchg: number[] = new Array(nAtoms).fill(0);
  let inMol = false;
  for (const line of mmdText.split('\n')) {
    const head = line.match(/^\s*\d+\s+\[(\w+),/);
    if (head) { inMol = head[1] === name; continue; }
    if (!inMol) continue;
    const p = line.trim().split(/\s+/);
    // atom line: bmin_type + 6 neighbor pairs + x y z label idx fchg pchg name subname serial
    if (p.length >= 20) {
      const serial = parseInt(p[p.length - 1], 10);
      if (!isNaN(serial)) pchg[serial - 1] = parseFloat(p[p.length - 4]);
    }
  }
  return pchg;
}
