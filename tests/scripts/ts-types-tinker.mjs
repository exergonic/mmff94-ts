// Feed Tinker OUR atom types (the 761/761 OpenBabel-validated set,
// class -> first-original-numbering via mmff94.prm) and minimize.
// Any energy difference vs our own minimizer is then a REAL
// optimizer/force-field-implementation difference, not typing.
//
// Usage: node ts-types-tinker.mjs [mol_006 mol_008 ...]  (default: all 29)
import { readFileSync, readdirSync, mkdirSync, writeFileSync } from "fs";
import { execFileSync } from "child_process";
import { join } from "path";
import { parse_sdf } from "/home/bwayne/src/mmff94-ts/dist/sdf.js";
import { assign_atom_types } from "/home/bwayne/src/mmff94-ts/dist/mmff94/assign-atom-types.js";

const SDF_DIR = "/home/bwayne/ff-bench/data/typing";
const OUT = "/tmp/ts_types_tinker";
const PRM = "/home/bwayne/src/tinker/params/mmff94.prm";
const MINIMIZE = "/usr/local/bin/minimize";
mkdirSync(OUT, { recursive: true });

// class -> first original number (same mapping as gen-tinker-fixtures.ts)
function class_to_original(prmText) {
  const map = new Map();
  for (const line of prmText.split("\n")) {
    const m = line.match(/^atom\s+(\d+)\s+(\d+)\s/);
    if (!m) continue;
    const orig = parseInt(m[1], 10);
    const cls = parseInt(m[2], 10);
    if (!map.has(cls)) map.set(cls, orig);
  }
  return map;
}
const prmText = readFileSync(PRM, "utf-8");
const map = class_to_original(prmText);

const names = process.argv.slice(2).length
  ? process.argv.slice(2)
  : readdirSync(SDF_DIR)
      .filter((f) => f.endsWith(".sdf") && f !== "mol_017.sdf")
      .map((f) => f.replace(".sdf", ""))
      .sort();

console.log("name,ours_E,tinker_E,tinker_ms,tinker_iters,delta");
for (const name of names) {
  // Our side: parse, type, minimize with the default criterion.
  const mol = parse_sdf(readFileSync(join(SDF_DIR, `${name}.sdf`), "utf-8"));
  const typed = assign_atom_types(mol);

  // Re-minimize a copy so the comparison geometry comes from THIS run.
  const { optimize_lbfgs } = await import("/home/bwayne/src/mmff94-ts/dist/optimize/l-bfgs.js");
  const work = JSON.parse(JSON.stringify(typed));
  const r = optimize_lbfgs(work, { max_iterations: 3000 });
  const minimized = r.converged ? r.molecule : optimize_lbfgs(typed, { max_iterations: 3000 }).molecule;

  // Tinker side: minimized geometry (from our run) + OUR types.
  const neighbors = Array.from({ length: minimized.atoms.length }, () => []);
  for (const b of minimized.bonds) {
    neighbors[b.atom1].push(b.atom2 + 1);
    neighbors[b.atom2].push(b.atom1 + 1);
  }
  const lines = [`${minimized.atoms.length} ${name}`];
  let badType = false;
  for (let i = 0; i < minimized.atoms.length; i++) {
    const a = minimized.atoms[i];
    const cls = typed.atom_types[i];
    const orig = map.get(cls);
    if (orig === undefined) { badType = true; break; }
    lines.push(
      `${String(i + 1).padStart(5)}  ${a.element.padEnd(2)} ${a.x.toFixed(6).padStart(13)} ${a.y.toFixed(6).padStart(13)} ${a.z.toFixed(6).padStart(13)} ${String(orig).padStart(5)}  ${neighbors[i].sort((x, y) => x - y).join(" ")}`
    );
  }
  if (badType) {
    console.log(`${name},SKIP(no-original-for-class),,,`);
    continue;
  }
  const txyz = join(OUT, `${name}.txyz`);
  const key = join(OUT, `${name}.key`);
  writeFileSync(txyz, lines.join("\n") + "\n");
  writeFileSync(key, `parameters ${PRM}\nMMFF-PIBOND\n`);

  // Minimize from OUR converged geometry: if the force fields agree,
  // Tinker's energy change should be ~0 (already at a minimum).
  let out;
  try {
    out = execFileSync(MINIMIZE, [txyz, key, "3000"], {
      input: "0.05\n",
      timeout: 300000,
      encoding: "utf-8",
    });
  } catch (e) {
    console.log(`${name},TINKER_FAIL,,,`);
    continue;
  }
  let eFinal = null;
  for (const line of out.split("\n")) {
    if (line.includes("Final Function Value")) eFinal = parseFloat(line.split(":")[1]);
  }
  const d = eFinal === null ? NaN : eFinal - r.energy.total;
  console.log(
    `${name},${r.energy.total.toFixed(4)},${eFinal === null ? "FAIL" : eFinal.toFixed(4)},,${d >= 0 ? "+" : ""}${d.toFixed(4)}`
  );
}
