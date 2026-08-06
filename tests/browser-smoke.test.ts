/**
 * Browser smoke test — the built dist/ loaded in a real browser.
 *
 * This is the regression guard for the browser target: dist/index.js is
 * plain ESM and must load in a browser WITHOUT a bundler. The historical
 * bug was extension-less relative imports emitted by tsc (browsers reject
 * them with a module-resolution error). The test rebuilds dist with tsc,
 * serves it over HTTP, loads it in headless Chromium (playwright driving
 * the system Edge via the msedge channel — no browser download), runs the
 * full pipeline in-page, and asserts the results are EXACTLY the numbers
 * the Node path computes (same V8 engine, same code), which are in turn
 * pinned against the reference logs.
 *
 * On a machine without Edge, swap the launch channel to 'chromium' and
 * run `npx playwright install chromium` once.
 */
import { describe, it, expect, beforeAll, afterAll } from 'vitest';
import { readFileSync } from 'fs';
import { join } from 'path';
import { spawnSync } from 'child_process';
import http from 'http';
import { chromium, type Browser, type Page } from 'playwright';

type Dist = typeof import('../dist/index.js');
let mmff: Dist;

// Fixture total energies from the obenergy reference logs (the same
// values reference-comparison.test.ts pins in Node).
const REFS: Record<string, number> = {
  benzene: 16.22697,
  pyridine: 15.5234,
  pyrrole: 3.2868,
  nicotine: 30.2543,
};

/** The pipeline under test, run in Node on the built dist. */
function inNode(sdfText: string, dist: Dist) {
  const mol = dist.parse_sdf(sdfText);
  const typed = dist.assign_bci_charges(dist.assign_atom_types(mol));
  const e = dist.calc_energy(typed);
  const g = dist.calc_gradient(typed);
  return {
    total: e.total,
    terms: [
      e.bond_stretch, e.angle_bend, e.stretch_bend, e.torsion,
      e.van_der_waals, e.electrostatic, e.out_of_plane,
    ],
    maxG: Math.max(...g.flat().map(Math.abs)),
    gradLen: g.length,
  };
}

/** The same pipeline inside the page, importing /dist/index.js over HTTP. */
function inPage(page: Page, sdfText: string) {
  return page.evaluate(async (sdf) => {
    // The historical regression: this dynamic import fails if dist
    // contains extension-less relative imports — browsers reject them
    // at resolution time. Built via new Function so vitest's file
    // transform does not rewrite the import into its SSR helper (which
    // only exists in the Node test process, not in the page).
    const dynamicImport = new Function('u', 'return import(u)') as (
      u: string,
    ) => Promise<typeof import('../../dist/index.js')>;
    const mmff = await dynamicImport('/dist/index.js');
    const mol = mmff.parse_sdf(sdf);
    const typed = mmff.assign_bci_charges(mmff.assign_atom_types(mol));
    const e = mmff.calc_energy(typed);
    const g = mmff.calc_gradient(typed);
    return {
      total: e.total,
      terms: [
        e.bond_stretch, e.angle_bend, e.stretch_bend, e.torsion,
        e.van_der_waals, e.electrostatic, e.out_of_plane,
      ],
      maxG: Math.max(...g.flat().map(Math.abs)),
      gradLen: g.length,
    };
  }, sdfText);
}

const MIME: Record<string, string> = {
  js: 'text/javascript', mjs: 'text/javascript', json: 'application/json',
  map: 'application/json', html: 'text/html', css: 'text/css',
};

/** A minimal static server for the dist directory (module MIME matters:
 *  .js must be text/javascript or the browser refuses the import). */
function serve(): http.Server {
  const distRoot = join(__dirname, '..', 'dist');
  return http.createServer((req, res) => {
    const url = new URL(req.url ?? '/', 'http://localhost');
    if (url.pathname === '/') {
      res.writeHead(200, { 'content-type': 'text/html' });
      res.end('<!doctype html><html><body>mmff94-ts browser smoke</body></html>');
      return;
    }
    const rel = url.pathname.replace(/^\/dist\//, '');
    const file = join(distRoot, rel);
    if (!file.startsWith(distRoot)) {
      res.writeHead(403);
      res.end();
      return;
    }
    try {
      const data = readFileSync(file);
      const ext = file.split('.').pop() ?? '';
      res.writeHead(200, { 'content-type': MIME[ext] ?? 'application/octet-stream' });
      res.end(data);
    } catch {
      res.writeHead(404);
      res.end();
    }
  });
}

describe('browser smoke — the built dist/ in headless Chromium', () => {
  let server: http.Server;
  let baseUrl: string;
  let browser: Browser;
  let page: Page;

  beforeAll(async () => {
    // Build dist fresh — the artifact under test must be current.
    const built = spawnSync('npx', ['tsc'], { cwd: join(__dirname, '..'), shell: true, encoding: 'utf-8' });
    if (built.status !== 0) {
      throw new Error(`tsc failed:\n${built.stdout}\n${built.stderr}`);
    }
    mmff = await import('../dist/index.js');

    server = serve();
    await new Promise<void>((resolve) => server.listen(0, '127.0.0.1', resolve));
    baseUrl = `http://127.0.0.1:${(server.address() as { port: number }).port}`;

    browser = await chromium.launch({ channel: 'msedge', headless: true });
    page = await browser.newPage();
    await page.goto(baseUrl + '/');
  }, 180_000);

  afterAll(async () => {
    await browser?.close();
    server?.close();
  });

  for (const [name, ref] of Object.entries(REFS)) {
    it(`dist loads in the browser and matches the reference for ${name}`, async () => {
      const sdf = readFileSync(join(__dirname, 'fixtures', 'sdf', `${name}.sdf`), 'utf-8');

      // dist-in-Node vs the reference log: absolute correctness.
      const node = inNode(sdf, mmff);
      expect(Math.abs(node.total - ref)).toBeLessThan(0.001);

      // dist-in-browser vs dist-in-Node: the same code, but Node 26's
      // V8 and Edge's V8 may differ in the last ULP of a Math.*
      // transcendental, so the comparison is 1e-7 absolute — a real
      // regression is orders of magnitude larger than that (the
      // historical dist bug was a hard load failure).
      const pageResult = await inPage(page, sdf);
      expect(Math.abs(pageResult.total - node.total)).toBeLessThan(1e-7);
      for (let t = 0; t < 7; t++) {
        expect(Math.abs(pageResult.terms[t] - node.terms[t])).toBeLessThan(1e-7);
      }
      expect(Math.abs(pageResult.maxG - node.maxG)).toBeLessThan(1e-7);
      expect(pageResult.gradLen).toBe(node.gradLen);
    }, 60_000);
  }
});
