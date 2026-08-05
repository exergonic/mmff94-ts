// Print the bmin refs for the metal hydrates.
import { readFileSync } from 'fs';

function parse_fortran(s: string): number {
  s = s.trim();
  if (s.includes('D')) s = s.replace('D', 'e');
  return parseFloat(s);
}

const text = readFileSync('tests/fixtures/validation-suite/MMFF94_bmin.log', 'utf-8');
const lines = text.split('\n');
let currentCode = '';
const result = new Map<string, Record<string, number>>();
for (let i = 0; i < lines.length; i++) {
  const m = lines[i].match(/\[\s*(\w+),/);
  if (m) { currentCode = m[1]; continue; }
  if (!currentCode) continue;
  const s = lines[i].match(/^\s+Stretch\s*=\s*(\S+)/);
  if (s) {
    const e: Record<string, number> = {};
    e.stretch = parse_fortran(s[1]);
    e.bend = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
    e.torsion = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
    e.oop = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
    e.strbnd = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
    e.elec = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
    e.vdw = parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]);
    result.set(currentCode, e);
    currentCode = '';
  }
}
for (const c of ['FE2PW3', 'CU1PW1', 'QUICNA01', 'SAHSUP']) {
  console.log(c, JSON.stringify(result.get(c)));
}
// extra: print any code given on argv
const args = process.argv.slice(2);
if (args.length) {
  const { assign_atom_types } = await import('../../src/mmff94/assign-atom-types');
  const { calc_energy } = await import('../../src/mmff94/energy/total');
  const { parse_mmd } = await import('../../src/utils/mmd-parser');
  const mols = parse_mmd(readFileSync('tests/fixtures/validation-suite/MMFF94.mmd', 'utf-8'));
  for (const n of args) {
    const ref = result.get(n);
    const mol = mols.find(m => m.name === n);
    if (mol) {
      const e = calc_energy(assign_atom_types(mol));
      console.log(n, 'ref oop', ref?.oop?.toFixed(6), 'ours oop', e.out_of_plane.toFixed(6), 'd', ref ? (e.out_of_plane - ref.oop).toFixed(6) : '-');
    }
  }
}
