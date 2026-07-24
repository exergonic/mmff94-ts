#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
PATH="/c/Program Files/OpenBabel-3.1.1:$PATH"
BABEL_DATADIR="/c/Users/mccan/AppData/Roaming/OpenBabel-3.1.1/data"
export BABEL_DATADIR

input="$1"
base="$(basename "$input" .sdf)"
base="${base%.SDF}"
ref_dir="$script_dir/../references"
mkdir -p "$ref_dir"
output="$ref_dir/${base}.log"

obabel "$input" -otxt --ff mmff94 --energy --log > "$output" 2>&1
