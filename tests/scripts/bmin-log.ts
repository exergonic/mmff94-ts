// Shared parser for the suite's BatchMin per-component log
// (MMFF94_bmin.log). Single source of truth for the reference per-term
// energies — used by the validation tests and the documentation
// generator, so the numbers they print cannot drift apart.
import { readFileSync } from 'fs';
import { join } from 'path';

export interface BminComponentEnergies {
  stretch: number;
  bend: number;
  torsion: number;
  oop: number;
  strbnd: number;
  elec: number;
  vdw: number;
  total: number;
}

function parse_fortran(s: string): number {
  s = s.trim();
  if (s.includes('D')) s = s.replace('D', 'e');
  return parseFloat(s);
}

export function parse_bmin_log(text: string): Map<string, BminComponentEnergies> {
  const result = new Map<string, BminComponentEnergies>();
  let currentCode = '';
  const lines = text.split('\n');

  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(/\[\s*(\w+),/);
    if (m) { currentCode = m[1]; continue; }
    if (!currentCode) continue;

    const s = lines[i].match(/^\s+Stretch\s*=\s*(\S+)/);
    if (s) {
      const e: BminComponentEnergies = {
        stretch: parse_fortran(s[1]),
        bend: parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]),
        torsion: parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]),
        oop: parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]),
        strbnd: parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]),
        elec: parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]),
        vdw: parse_fortran(lines[++i].match(/=\s*(\S+)/)![1]),
      };
      const t = lines[++i].match(/Total Energy\s*=\s*([-0-9.]+)/);
      // NOTE: the log's "Total Energy" line is printed 2-decimal
      // rounded (Fortran F9.2 — e.g. -51.42700), while the component
      // lines are 5-decimal. Do NOT use `total` at fine precision;
      // compare totals against MMFF94.energies instead (5 decimals).
      if (t) e.total = parseFloat(t[1]);
      result.set(currentCode, e);
      currentCode = '';
    }
  }
  return result;
}

export function load_bmin_log(suiteDir = 'tests/fixtures/validation-suite'): Map<string, BminComponentEnergies> {
  return parse_bmin_log(readFileSync(join(suiteDir, 'MMFF94_bmin.log'), 'utf-8'));
}
