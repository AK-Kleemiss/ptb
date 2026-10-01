#!/usr/bin/env bash
# Run in a clean directory containing source/ on a Mac with both Homebrew GCC 16 toolchains.
set -euo pipefail
test "$(uname -s)" = Darwin
for arch_name in arm64 x86_64; do
  if [ "$arch_name" = arm64 ]; then
    prefix=/opt/homebrew
  else
    prefix=/usr/local
  fi
  fc="$prefix/bin/gfortran-16"
  cc="$prefix/bin/gcc-16"
  gcc_lib="$prefix/opt/gcc/lib/gcc/current"
  test -x "$fc" && test -x "$cc"
  if [ -e "src_$arch_name" ]; then
    echo "src_$arch_name already exists; use a clean build directory" >&2
    exit 1
  fi
  mkdir -p "src_$arch_name" "build_$arch_name" "staticlibs_$arch_name"
  cp -R source/. "src_$arch_name/"
  for lib in libgomp libgfortran libquadmath libgcc libgcc_eh; do
    if [ -f "$gcc_lib/$lib.a" ]; then
      ln -s "$gcc_lib/$lib.a" "staticlibs_$arch_name/$lib.a"
    fi
  done
  abs_static="$(pwd)/staticlibs_$arch_name"
  if [ "$arch_name" = x86_64 ]; then
    arch_cmd=(arch -x86_64)
  else
    arch_cmd=(arch -arm64)
  fi
  "${arch_cmd[@]}" make -C "src_$arch_name" NAME=LINUX COMPILER=gfortran \
    PROG="../build_$arch_name/ptb" FC="$fc" CC="$cc" \
    FFLAGS="-O2 -ffree-line-length-none -fopenmp" \
    CCFLAGS="-O2 -std=gnu17 -DLINUX" \
    LINKER="$fc -O2 -static-libgcc -static-libgfortran -L$abs_static" \
    LIBS="$abs_static/libgomp.a -framework Accelerate"
  file "build_$arch_name/ptb"
  otool -L "build_$arch_name/ptb"
  if otool -L "build_$arch_name/ptb" | grep -E '/(opt/homebrew|usr/local)/'; then
    echo "Homebrew runtime dependency in $arch_name build" >&2
    exit 1
  fi
done
mkdir -p build
lipo -create build_arm64/ptb build_x86_64/ptb -output build/ptb_macos
file build/ptb_macos
otool -L build/ptb_macos
shasum -a 256 build/ptb_macos
