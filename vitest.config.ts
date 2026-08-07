import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    include: ['tests/**/*.test.ts'],
    // The gradient FD suite is ~19,000 evaluations by design (~1 min
    // sequential) — a 2-minute budget keeps it honest without flaking.
    // The teardownTimeout covers the worker's task-update RPC, which
    // the long synchronous FD loops can starve.
    testTimeout: 120_000,
    teardownTimeout: 120_000,
    coverage: {
      provider: 'v8',
      include: ['src/**/*.ts'],
    },
  },
});
