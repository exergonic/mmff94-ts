#!/usr/bin/env node
/**
 * WebMO engine shim for mmff94-ts.
 *
 * Runs inside a WebMO job directory as the engine executable (the contract
 * of run_mmff94ts.cgi, a close relative of WebMO's run_tinker.cgi):
 *
 *   stdin   <- input.stdin     (the calculation mode: "energy" | "optimize")
 *   stdout  -> output.out      (a Tinker-style report WebMO parses)
 *   cwd     =  the job directory
 *
 * The structure comes from the editor's own exports, staged by WebMO:
 *   input.xyz    element + x/y/z, one atom per line
 *   connections  "i j order" per bond, 1-based, orders 1..3 (kekulized)
 *
 * The final geometry of an optimization is written to input.xyz_2 in the
 * Tinker layout WebMO's parser expects (index, element, x, y, z, type).
 */
import { readFileSync, writeFileSync } from 'node:fs';
import { calc_energy, diagnose_molecule, optimize_lbfgs } from '../../dist/index.js';

const VERSION = '0.1.0-alpha.2';

function fail(msg) {
  console.error(`mmff94-ts: ${msg}`);
  process.exit(1);
}

function read_xyz(path) {
  let text;
  try {
    text = readFileSync(path, 'utf-8');
  } catch {
    fail(`cannot read ${path}`);
  }
  const atoms = [];
  for (const line of text.split('\n')) {
    const t = line.trim().split(/\s+/);
    if (t.length < 4) continue;
    const m = t[0].match(/^([A-Z][a-z]?)$/);
    if (!m) continue;
    const x = Number(t[1]), y = Number(t[2]), z = Number(t[3]);
    if (![x, y, z].every(Number.isFinite)) continue;
    atoms.push({ index: atoms.length, element: m[1], x, y, z });
  }
  if (atoms.length === 0) fail('no atoms parsed from input.xyz');
  return atoms;
}

function read_connections(path, atom_count) {
  let text;
  try {
    text = readFileSync(path, 'utf-8');
  } catch {
    fail(`cannot read ${path}`);
  }
  const bonds = [];
  for (const line of text.split('\n')) {
    const t = line.trim().split(/\s+/);
    if (t.length < 3) continue;
    const i = parseInt(t[0], 10), j = parseInt(t[1], 10), order = parseInt(t[2], 10);
    if (!Number.isInteger(i) || !Number.isInteger(j)) continue;
    if (i < 1 || j < 1 || i > atom_count || j > atom_count) continue;
    if (order < 1) continue;   // the editor lists non-bonded contacts as order 0
    if (order > 3) continue;   // kekulized structures carry orders 1..3 only
    bonds.push({ atom1: i - 1, atom2: j - 1, bond_order: order });
  }
  return bonds;
}

/**
 * The editor's per-atom formal charges: one line per atom, in input.xyz
 * order, the first column being the charge (the total is the job's net
 * charge). The typer needs them — a three-coordinate carbon carrying +1
 * is a carbocation, and MMFF94's typing only places it in the vinylic
 * sp2 class when the charge is visible. A missing file (older jobs)
 * leaves every atom neutral, exactly as before.
 */
function read_charges(path, atoms) {
  let text;
  try {
    text = readFileSync(path, 'utf-8');
  } catch {
    return 0;
  }
  const lines = text.split('\n').filter((line) => line.trim() !== '');
  let total = 0;
  for (let i = 0; i < Math.min(lines.length, atoms.length); i++) {
    const q = parseInt(lines[i].trim().split(/\s+/)[0], 10);
    if (!Number.isInteger(q) || q === 0) continue;
    atoms[i].formal_charge = q;
    total += q;
  }
  return total;
}

function banner() {
  const bar = '#' .repeat(66);
  return [
    '',
    `     ${bar}`,
    `   ${bar}####`,
    '  ###                                                                      ###',
    ' ###            mmff94-ts  ---  MMFF94 in pure TypeScript                  ###',
    ' ##                                                                          ##',
    ` ##                         Version ${VERSION}                                 ##`,
    ' ##                                                                          ##',
    ' ##       MMFF94 energy model of Halgren 1996 - energy, gradients, L-BFGS    ##',
    ' ###                                                                        ###',
    '  ###                                                                      ###',
    `   ${bar}####`,
    `     ${bar}`,
    '',
  ].join('\n');
}

/** The energy table WebMO's parser reads (Tinker's "Energy Component
 *  Breakdown"; mmff94-ts's terms map onto the Tinker names one-for-one). */
function component_block(e) {
  const rows = [
    ['Bond Stretching', e.bond_stretch],
    ['Angle Bending', e.angle_bend],
    ['Stretch-Bend', e.stretch_bend],
    ['Torsional Angle', e.torsion],
    ['Out-of-Plane Bend', e.out_of_plane],
    ['Van der Waals', e.van_der_waals],
    ['Charge-Charge', e.electrostatic],
  ];
  const lines = ['', ' Energy Component Breakdown :           Kcal/mole', ''];
  for (const [name, value] of rows) {
    lines.push(` ${name.padEnd(38)} ${value.toFixed(4).padStart(10)}`);
  }
  lines.push(` ${'Total'.padEnd(38)} ${e.total.toFixed(4).padStart(10)}`);
  return lines.join('\n');
}

/** input.xyz_2 in Tinker's layout: index, element, x, y, z, type, neighbours. */
function final_geometry(molecule) {
  const adjacency = new Map();
  for (const b of molecule.bonds) {
    if (!adjacency.has(b.atom1)) adjacency.set(b.atom1, []);
    if (!adjacency.has(b.atom2)) adjacency.set(b.atom2, []);
    adjacency.get(b.atom1).push(b.atom2 + 1);
    adjacency.get(b.atom2).push(b.atom1 + 1);
  }
  const lines = [`${String(molecule.atoms.length).padStart(6)}  WebMO Mechanics`];
  for (let i = 0; i < molecule.atoms.length; i++) {
    const a = molecule.atoms[i];
    const nbrs = (adjacency.get(i) ?? []).sort((x, y) => x - y).join(' ');
    lines.push(
      `${String(i + 1).padStart(6)}  ${a.element.padEnd(2)} ` +
      `${a.x.toFixed(6).padStart(12)} ${a.y.toFixed(6).padStart(12)} ${a.z.toFixed(6).padStart(12)}` +
      `  ${String(molecule.atom_types?.[i] ?? 0).padStart(5)}  ${nbrs}`,
    );
  }
  return lines.join('\n') + '\n';
}

/** The parameter-gap disclosure, in the report the user actually reads.
 *  diagnose_molecule() costs one extra energy evaluation; on a WebMO job
 *  that is the right trade — a silently dropped interaction would
 *  otherwise leave the results page looking complete. */
function diagnostics_block(d) {
  const empirical = Object.entries(d.empirical)
    .filter(([, n]) => n > 0)
    .map(([k, n]) => `${n} ${k.replace(/_/g, '-')}`)
    .join(', ');
  const lines = ['', ' Parameter Diagnostics :', ''];
  lines.push(
    ` ${'Coordination gaps'.padEnd(20)} ${d.atoms.length === 0 ? 'none' : `${d.atoms.length} atom(s): ` + d.atoms.map((a) => `${a.element}${a.index + 1} (type ${a.type}, ${a.coordination} neighbours)`).join(', ')}`,
  );
  lines.push(
    ` ${'Empirical rules'.padEnd(20)} ${empirical === '' ? 'none' : `${empirical} (part V, expected behaviour)`}`,
  );
  const omitted = Object.entries(d.dropped).filter(([k, n]) => k !== 'torsion' && n > 0);
  const strictDropped = omitted.reduce((a, [, n]) => a + n, 0);
  lines.push(
    ` ${'Omitted by rules'.padEnd(20)} ${d.dropped.torsion === 0 ? 'none' : `${d.dropped.torsion} torsion path(s) (linear centres / unsaturated-sp2)`}`,
  );
  lines.push(
    ` ${'Dropped (no rule)'.padEnd(20)} ${strictDropped === 0 ? 'none' : omitted.map(([k, n]) => `${n} ${k.replace(/_/g, '-')}`).join(', ')}`,
  );
  if (d.atoms.length > 0 || strictDropped > 0) {
    lines.push('');
    lines.push(' ***  PARAMETER GAPS ABOVE: the result is not fully parameterized  ***');
  }
  return lines.join('\n');
}

function main() {
  let stdin_text = '';
  try {
    stdin_text = readFileSync('input.stdin', 'utf-8');
  } catch {
    /* no stdin file: fall through to the default */
  }
  const mode = stdin_text.trim().split(/\s+/)[0] || 'energy';

  const atoms = read_xyz('input.xyz');
  const bonds = read_connections('connections', atoms.length);
  const net_charge = read_charges('charges', atoms);
  const molecule = { atoms, bonds, name: 'WebMO job' };

  const out = [banner()];
  if (net_charge !== 0) {
    out.push(` Net Formal Charge :                    ${net_charge > 0 ? '+' : ''}${net_charge} e`);
  }

  if (mode === 'energy') {
    const e = calc_energy(molecule);
    out.push(` Total Potential Energy :                ${e.total.toFixed(4)} Kcal/mole`);
    out.push(component_block(e));
  } else if (mode === 'optimize') {
    const result = optimize_lbfgs(molecule, { gradient_tolerance: 0.05 });
    const verdict = result.converged ? 'Normal Termination due to SmallGrad' : 'Termination: Criterion Not Met';
    out.push('');
    out.push(` LBFGS  --  ${verdict}`);
    out.push('');
    out.push(` Final Function Value :    ${result.energy.total.toFixed(4).padStart(14)}`);
    out.push(` Final RMS Gradient :       ${(result.final_rms_gradient ?? NaN).toFixed(4).padStart(14)}`);
    out.push(component_block(result.energy));
    writeFileSync('input.xyz_2', final_geometry(result.molecule));
  } else {
    fail(`unknown mode '${mode}' in input.stdin (expected 'energy' or 'optimize')`);
  }

  out.push(diagnostics_block(diagnose_molecule(molecule, { quiet: true })));
  out.push('');
  out.push(' mmff94-ts -- MMFF94 (Halgren J. Comput. Chem. 1996) in pure TypeScript');
  out.push('');
  process.stdout.write(out.join('\n') + '\n');
}

try {
  main();
} catch (err) {
  fail(err && err.message ? err.message : String(err));
}
