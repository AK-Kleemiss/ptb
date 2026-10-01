#!/usr/bin/env bash
# Usage: bash smoke_yb_defaults.sh /absolute/path/to/ptb /scratch/output [arch]
set -euo pipefail
exe=$1
out=$2
arch_name=${3:-}
for q in 0 2 3; do
  n=$((24-q))
  u=$((n%2))
  work="$out/q$q"
  mkdir -p "$work"
  printf '1\nYb release smoke\nYb 0 0 0\n' > "$work/atom.xyz"
  if [ -n "$arch_name" ]; then
    (cd "$work" && arch "-$arch_name" "$exe" atom.xyz -chrg "$q" -uhf "$u" -stda -denmat ptb.denmat > ptb.log 2>&1)
  else
    (cd "$work" && "$exe" atom.xyz -chrg "$q" -uhf "$u" -stda -denmat ptb.denmat > ptb.log 2>&1)
  fi
  test -s "$work/ptb.denmat"
  test -s "$work/wfn.xtb"
  echo "Yb q$q passed"
done
