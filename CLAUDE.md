# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PTB is a density-matrix tight-binding (P-TB) quantum chemistry method (Grimme, Müller, Hansen, J. Chem. Phys. 2023) implemented in Fortran. It computes electronic structure properties (populations, dipoles, polarizabilities, Wiberg bond orders, etc.) from a polarized valence double-zeta basis for elements up to Z=86.

This is a fork (`AK-Kleemiss/ptb`) of the upstream `grimme-lab/ptb` repository, modified for integration into NoSpherA2 without altering the original PTB model — see the version banner in `source/main.f90`'s `head` subroutine. Keep upstream-compatibility in mind: prefer additive/behind-flag changes over edits that alter core PTB physics.

## Build

```
cd source
make
```

Produces `../build/ptb`. The Makefile auto-detects OS (`Windows_NT`/`Linux`/`Darwin` via `$(OS)`/`uname -s`) and selects a compiler/linker/library set accordingly:

- **Windows**: Intel `ifx`/`icx`, links against Intel MKL (`$(MKLROOT)`), static link, OpenMP (`-Qiopenmp`).
- **Linux**: defaults to `ifx` + MKL; set `COMPILER=gfortran` to use `gfortran`/`gcc` + system LAPACK/BLAS instead (`make COMPILER=gfortran`).
- **Mac**: unsupported (no Intel Fortran compiler available).

`make clean` removes `*.o`, `*.mod`, `*genmod*`, and the built binary.

There is no test suite or linter configured in this repo — validation is done by running the binary against reference/DFT data (see Runtime inputs below).

### Adding a new source file
Object files must be added explicitly to `OBJS1` (or `MODULES` if it's a module needed early) in `source/Makefile`; there's no glob-based source discovery. Respect the manual dependency ordering rules near the bottom of the Makefile (e.g. `main.o: pgtb.o`, `purification.o: purification_settings.o symmetry_i.o norms.o metrics.o`) — Fortran module compilation order matters.

## Running

```
ptb <geometry-file> [options]
```

Geometry input (first positional arg) accepts Turbomole `coord` or `.xyz` format. Two data files are required at runtime:
- `.atompara` — empirical element parameters, default `~/.atompara`, override with `-par <path>`
- `.basis_vDZP` — vDZP basis set, default `~/.basis_vDZP`, override with `-bas <path>`

Both example files exist at the repo root (`.atompara`, `.basis_vDZP`), along with `.ecp` and `.auxbasis_vDZP`.

Other runtime file-based inputs (read from CWD, not CLI flags): `.CHRG` (total charge), `.UHF` (# unpaired electrons), `.EFIELD` (external field), `.RPBE`. Key CLI flags: `-chrg`, `-uhf`, `-stda`, `-test`, `-nogtb`, `-purify`, `-json [file]`, `-version`, `-help`.

## Architecture

### Entry point and flow
`source/main.f90` is a single large procedural driver (not modular by design — this mirrors upstream Grimme-group code style). Flow: parse CLI args → read `.atompara`/basis files → build AO/core basis (`bascom`, `cbascom`) → compute overlap integrals (`sint`) → optionally read DFT reference data (`rdtm`, only if not `-stda`/`.RPBE`) → run the PTB single-point calculation via `pgtb` (in `pgtb.f90`) → post-process (population analysis, Wiberg bond orders, dipoles, polarizabilities) → optionally emit a `.data` fit-file (used by the parameter-fitting workflow) and/or JSON output (`json_output.f90`).

Large chunks of `main.f90` (gradient/dipole-gradient/Raman finite-difference branches, D4-only dispersion path, energy-only path) are `include`d from separate files (`grad.f90`, `dipgrad.f90`, `polgrad.f90`, `testout.f90`) rather than being subroutines — these share `main`'s local variable scope directly via Fortran `include`.

### Shared-state modules (old-school Fortran globals)
Several modules hold plain public module-level arrays instead of derived types, and are `use`d throughout for shared state:
- `bascom` — AO basis bookkeeping (shell counts, exponents/contractions, shell↔AO maps)
- `cbascom` — analogous core-basis bookkeeping
- `parcom` — empirical TB parameters read from `.atompara` (indexed `(shell, element Z)`)
- `com` — general per-element constants (atomic mass, valence electrons, EN, hardness, avg. CN)
- `mocom` — reference MO/density-matrix data used for parameter fitting and MO matching
- `thresholds`, `metrics`, `norms` — numerical thresholds and matrix-norm helpers used by the purification code

### Newer object-oriented modules
`purification_settings.f90` and `purification.f90` follow a more modern Fortran style (derived types with type-bound procedures, e.g. `tPurificationSet`, `tMetricSet`, `tChempotSet`). This is the density-matrix purification engine — an alternative to explicit diagonalization for computing the density matrix from the Fock matrix, with 4 selectable schemes (`mcweeny`, iterative sign via Padé or Newton-Schulz, or plain diagonalization) and a selectable S-metric power (`inv`, `inv_sqrt`, `sqrt`). Enabled via the `-purify` CLI flag, which allocates a `tPurificationSet` (`pur`) threaded through into `pgtb`. Settings can be read from a `.PUR` control file via `%settings()` (currently commented out in `main.f90`).

### Vendored third-party libraries
Two large files are single-file amalgamations of external open-source libraries, vendored wholesale rather than as separate modules per file — do not try to split them apart, and check the corresponding upstream project before making non-trivial changes:
- `dftd4.f90` (~12.9k lines): bundles `dftd4`, `multicharge`, and `mctc-lib` (many `mctc_*`/`multicharge_*`/`dftd4_*` modules) for DFT-D4 dispersion and EEQ charge models. `dftd4.patch` documents a manual patch applied on top of the vendored source (self-consistent lattice allocation fix).
- `la.f90` + `blas_level1/2/3.f90` + `lapack_eig.f90`: a linear-algebra shim module (`gtb_la`, `gtb_la_level1/2/3`, `gtb_lapack_eig`) wrapping BLAS/LAPACK calls, compiled early via the `SOURCES` list in the Makefile since later modules depend on it.

### GPU/accelerator code (currently disabled)
`cuda_context.f90` (module `cuda_`, using `accel_lib`) and commented-out `-cuda` CLI handling in `main.f90` are experimental CUDA offload scaffolding. Not included in the Makefile `OBJS`/`SOURCES` lists and not currently built — treat as inactive/WIP unless asked to wire it back in.

### Naming/style conventions to preserve
- Fixed-form `.f` files (older, F77-style) coexist with free-form `.f90` files (newer). Match the style of whichever file you're editing.
- `real(wp)` with `wp => real64` from `iso_fortran_env` is the standard real kind in newer modules; older `.f` code uses `real*8`/`real*4` directly.
- Element-indexed arrays are conventionally sized `(86)` or `(10,86)` (10 = max shells per atom, 86 = max Z supported); don't change these bounds casually since many modules assume them.
