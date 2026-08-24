// Minimal-repro hunt: small neutral N/O/S molecules -> our types -> txyz
// -> tinker analyze [M]. Which functional groups trigger a spurious net
// charge in Tinker's MMFF94 charge derivation?
//
// Usage: node minrepro.mjs   (writes /tmp/minrepro/*.txyz, prints table)
import { readFileSync, writeFileSync, mkdirSync } from "fs";
import { execFileSync } from "child_process";
import { join } from "path";

const PRM = "/home/bwayne/src/tinker/params/mmff94.prm";
const MINIMIZE_DIR = "/tmp/minrepro";
mkdirSync(MINIMIZE_DIR, { recursive: true });

const { parse_smiles_or_die, assign_atom_types, assign_bci_charges } =
  await import("/home/bwayne/src/mmff94-ts/dist/index.js").catch(() => ({}));

// The lib may not export a SMILES parser; use RDKit via python for embeds,
// then parse the SDF here. Simpler: pre-embedded SDFs made by RDKit below.
const MOLECULES = {
  pyridine: "c1ccncc1",
  imidazole: "c1cnc[nH]1",
  indole: "c1ccc2[nH]ccc2c1",
  aniline: "Nc1ccccc1",
  Nmethylacetamide: "CC(=O)NC",
  methanesulfonamide: "CS(=O)(=O)N",
  dimethylsulfone: "CS(=O)(=O)C",
  urea: "NC(=O)N",
  acetophenone: "CC(=O)c1ccccc1",
};

function class_to_original(prmText) {
  const m = new Map();
  for (const line of prmText.split("\n")) {
    const parts = line.trim().split(/\s+/);
    if (parts.length >= 3 && parts[0] === "atom") {
      const orig = parseInt(parts[1], 10);
      const cls = parseInt(parts[2], 10);
      if (!m.has(cls)) m.set(cls, orig);
    }
  }
  return m;
}

const map = class_to_original(readFileSync(PRM, "utf-8"));
const { parse_sdf } = await import("/home/bwayne/src/mmff94-ts/dist/sdf.js");
const { assign_atom_types: type_atoms } = await import("/home/bwayne/src/mmff94-ts/dist/mmff94/assign-atom-types.js");
const { assign_bci_charges: charge_atoms } = await import("/home/bwayne/src/mmff94-ts/dist/mmff94/charges.js");

console.log("name,type_sum_check,tinker_net_charge");
for (const [name, smiles] of Object.entries(MOLECULES)) {
  // 3D embed via obabel (available on lenovo)
  const sdf = execFileSync(
    "/usr/bin/obabel",
    [`-:` + smiles, `-osdf`, `--gen3d`, `best`],
    { encoding: "utf-8", timeout: 60000 },
  );
  const sdfPath = join(MINIMIZE_DIR, `${name}.sdf`);
  writeFileSync(sdfPath, sdf);

  const mol = parse_sdf(sdf);
  const typedCharged = charge_atoms(type_atoms(mol));
  const sum = typedCharged.partial_charges.reduce((a, b) => a + b, 0);

  // write txyz with OUR types on the embedded geometry
  const neighbors = Array.from({ length: mol.atoms.length }, () => []);
  for (const b of mol.bonds) {
    neighbors[b.atom1].push(b.atom2 + 1);
    neighbors[b.atom2].push(b.atom1 + 1);
  }
  const lines = [`${mol.atoms.length} ${name}`];
  let ok = true;
  for (let i = 0; i < mol.atoms.length; i++) {
    const a = mol.atoms[i];
    const orig = map.get(typedCharged.atom_types[i]);
    if (orig === undefined) { ok = false; break; }
    lines.push(
      `${String(i + 1).padStart(5)}  ${a.element.padEnd(2)} ${a.x.toFixed(6).padStart(13)} ${a.y.toFixed(6).padStart(13)} ${a.z.toFixed(6).padStart(13)} ${String(orig).padStart(5)}  ${neighbors[i].sort((x, y) => x - y).join(" ")}`
    );
  }
  if (!ok) {
    console.log(`${name},SKIP(no-original),`);
    continue;
  }
  const stem = join(MINIMIZE_DIR, name);
  writeFileSync(`${stem}.txyz`, lines.join("\n") + "\n");
  // Tinker's suffix() looks for <stem>.xyz (not .txyz) — mirror the file
  writeFileSync(`${stem}.xyz`, lines.join("\n") + "\n");
  writeFileSync(`${stem}.key`, `parameters ${PRM}\nMMFF-PIBOND\n`);

  let out;
  try {
    out = execFileSync("/usr/local/bin/analyze", [stem], {
      input: "M\n", encoding: "utf-8", timeout: 120000,
    });
  } catch (e) {
    console.log(`${name},${sum.toFixed(4)},ANALYZE_FAIL`);
    continue;
  }
  let net = "?";
  for (const line of out.split("\n")) {
    if (line.includes("Total Electric")) net = line.split(":")[1].trim().split()[0];
  }
  console.log(`${name},${sum.toFixed(4)},${net}`);
}
