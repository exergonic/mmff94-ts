#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
PATH="/c/Program Files/OpenBabel-3.1.1:$PATH"
# Native (C:/) form on purpose: obabel is a native binary and this shell does
# not translate MSYS paths in the environment, so /c/Users/... reaches it as a
# literal string and the MMFF94 parameter files are not found.
BABEL_DATADIR="C:/Users/mccan/AppData/Roaming/OpenBabel-3.1.1/data"
export BABEL_DATADIR

input="$1"
base="$(basename "$input" .sdf)"
base="${base%.SDF}"
ref_dir="$script_dir/../references"
mkdir -p "$ref_dir"
output="$ref_dir/${base}.mmff94.log"

obabel "$input" -otxt --ff mmff94 --energy --log 2>&1 | tr -d '\r' > "$output"
