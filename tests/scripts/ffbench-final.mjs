// Final: per-molecule time-to-converged with restart+jitter, vs single run.
import { readFileSync, readdirSync, writeFileSync } from "fs";
import { parse_sdf } from "./dist/sdf.js";
import { optimize_lbfgs } from "./dist/optimize/l-bfgs.js";
const SDF = "/home/bwayne/ff-bench/data/typing";
const rows = [];
let tSingle = 0, tRestart = 0, okS = 0, okR = 0;
for (const f of readdirSync(SDF).filter(x => x.endsWith(".sdf") && x !== "mol_017.sdf").sort()) {
  const name = f.replace(".sdf", "");
  // single run
  let mol = parse_sdf(readFileSync(SDF + "/" + f, "utf-8"));
  let t0 = performance.now();
  let r = optimize_lbfgs(mol, { max_iterations: 3000 });
  let ms1 = performance.now() - t0;
  tSingle += ms1; if (r.converged) okS++;
  const eSingle = r.energy.total;
  // restart+jitter rounds
  let cur = mol, conv = false, t1 = performance.now(), eFinal = r.energy.total;
  for (let round = 0; round < 6 && !conv; round++) {
    const rr = optimize_lbfgs(cur, { max_iterations: 500 });
    conv = rr.converged;
    cur = { ...rr.molecule, bonds: rr.molecule.bonds };
    eFinal = rr.energy.total;
    if (!conv) {
      let s = 42 + round;
      for (const a of cur.atoms) for (const k of ["x","y","z"]) {
        s = (s * 1103515245 + 12345) & 0x7fffffff;
        a[k] += ((s / 0x7fffffff) - 0.5) * 0.01;
      }
    }
  }
  let ms2 = performance.now() - t1;
  tRestart += ms2; if (conv) okR++;
  rows.push(`${name},${r.converged ? "OK" : "stuck"},${ms1.toFixed(0)},${conv ? "OK" : "stuck"},${ms2.toFixed(0)},${eSingle.toFixed(4)},${eFinal.toFixed(4)}`);
}
writeFileSync("/home/bwayne/ff-bench/data/ts_mmff_fast/final_timing.csv",
  "name,single_status,single_ms,restart_status,restart_ms,E_single,E_restart\n" + rows.join("\n") + "\n");
console.log("single-run:  " + okS + "/29 converged, avg " + (tSingle / 29).toFixed(0) + " ms/molecule (3000-iter cap)");
console.log("restart+jit: " + okR + "/29 converged, avg " + (tRestart / 29).toFixed(0) + " ms/molecule (6x500)");
