#!/usr/bin/env bash
# Run in a clean directory containing source/ on Linux x86_64 with oneAPI.
set -eo pipefail
test "$(uname -s)" = Linux
test "$(uname -m)" = x86_64
export MKLROOT=${MKLROOT:-/opt/intel/oneapi/mkl/latest}
mkdir -p build
make -C source
file build/ptb
file build/ptb | grep -q 'statically linked'
ldd build/ptb || true # ldd exits 1 for a fully static executable
sha256sum build/ptb
