# WebMO integration

Runs mmff94-ts as a WebMO computation engine. WebMO's editor supplies the
structure (ordered bonds included) and mmff94-ts performs its own atom
typing, charges, and energy evaluation — which makes WebMO a front-end test
harness for the engine itself.

## Files here

- `mmff94ts-run.mjs` — the shim. Runs inside a WebMO job directory; reads
  `input.xyz` (elements + coordinates), `connections` (bond list with orders),
  and `input.stdin` (mode: `energy` or `optimize`); writes a Tinker-style
  report to stdout and, for optimizations, the final geometry to
  `input.xyz_2`.
- `mmff94ts` — the executable WebMO invokes. A one-line wrapper that pins an
  absolute node path for Apache's minimal environment.

## WebMO-side files

The engine pack lives in the WebMO installation tree (cgi-bin/webmo and
webmo/), adapted from WebMO's Tinker pack: `interfaces/mmff94ts.int`,
`interfaces/mmff94ts.tmpl`, `run_mmff94ts.cgi`, `parse_mmff94ts.cgi`,
`mmff94ts.cgi`, `mmff94ts.html`, `javascript/mmff94ts.js`. For a new
installation, copy WebMO's tinker pack files and rename tinker -> mmff94ts,
then apply the corresponding edits (no MMFF94 preparatory retyping is
performed: the shim consumes the editor's exports directly).

## Deploy / update

From the repository root:

```
npm run build
tar czf mmff94ts-deploy.tgz dist integrations package.json
```

Extract into the deployment directory (any path; `mmff94ts.int` points to
its `bin/` subdirectory via `mmff94tsBinDir`), then make the wrapper
executable. `package.json` must ship alongside `dist/` — it carries
`"type": "module"`, without which Node parses the build as CommonJS.

## Validation

At the same start geometry, mmff94-ts and Tinker agree to the precision
WebMO displays: nitrobenzene optimization gives 49.4720 kcal/mol here vs
49.4719 with Tinker 26.2, and the single-point energies agree to
0.0001 kcal/mol.
