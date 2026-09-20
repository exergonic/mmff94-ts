#!/usr/bin/env bash
# Prove that the committed parameter tables are exactly what the extractor
# produces — the "generated, do not edit by hand" guarantee for
# src/mmff94/parameters/*.ts. The hand-folded sulfinate/metal bridges live
# in the extractor's BRIDGES block, so they are covered too.
#
# Skips (exit 0, with a note) when the OpenBabel .par sources are absent:
# temp_ob/data is a local working copy, not part of the repository.
set -u
cd "$(dirname "$0")/.."

# Only the tables the extractor writes; empirical.ts, index.ts and
# parameter-classes.ts are hand-written and live in the same directory.
GENERATED="bond.ts angle.ts stretch-bend.ts torsion.ts van-der-waals.ts bci.ts
out-of-plane.ts atom-types.ts lookup.ts atom-type-properties.ts
default-stretch-bend.ts"

PAR_DIR="${MMFF94_PAR_DIR:-temp_ob/data}"
if [ ! -f "$PAR_DIR/mmffbond.par" ]; then
  echo "check-parameter-tables: $PAR_DIR has no .par sources - skipping (local-only data)."
  exit 0
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
# The extractor may be a native interpreter that cannot read MSYS paths.
if command -v cygpath > /dev/null 2>&1; then
  TMP_FOR_PYTHON="$(cygpath -w "$TMP")"
else
  TMP_FOR_PYTHON="$TMP"
fi

if ! MMFF94_PAR_OUT="$TMP_FOR_PYTHON" python3 scripts/extract-mmff94-par.py > "$TMP/extract.log" 2>&1; then
  echo "check-parameter-tables: the extractor FAILED:"
  cat "$TMP/extract.log"
  exit 1
fi

fail=0
for b in $GENERATED; do
  f="src/mmff94/parameters/$b"
  if [ ! -f "$TMP/$b" ]; then
    echo "  MISSING from the extractor's output: $b"
    fail=1
    continue
  fi
  if ! diff -q "$f" "$TMP/$b" > /dev/null; then
    echo "  DRIFT in $b (committed vs regenerated):"
    diff -u "$f" "$TMP/$b" | head -30
    fail=1
  fi
done

if [ "$fail" -ne 0 ]; then
  echo "check-parameter-tables: FAIL - the committed tables are not what the extractor produces."
  exit 1
fi
echo "check-parameter-tables: OK - all 11 generated tables match the extractor's output."
