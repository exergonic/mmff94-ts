/**
 * The parameter tables claim to be generated — every file in
 * src/mmff94/parameters/ says "Do not edit by hand". This gate proves it:
 * it runs the extractor into a temp directory and diffs every table.
 *
 * It passes with a printed skip when the OpenBabel .par sources are absent
 * (temp_ob/data is a local working copy, not part of the repository), so
 * it is a real gate on machines that carry the data and inert elsewhere.
 *
 * If this fails, either the committed tables were hand-edited (fold the
 * change into scripts/extract-mmff94-par.py) or the extractor changed.
 */
import { describe, it, expect } from 'vitest';
import { spawnSync } from 'child_process';
import { join } from 'path';

describe('parameter tables are generated, not hand-edited', () => {
  it('the extractor reproduces every committed table byte-for-byte', () => {
    const r = spawnSync('bash', ['scripts/check-parameter-tables.sh'], {
      cwd: join(__dirname, '..'),
      shell: true,
      encoding: 'utf-8',
    });
    const output = `${r.stdout ?? ''}${r.stderr ?? ''}`.trim();
    console.log(output);
    expect(r.status, output).toBe(0);
  }, 180_000);
});
