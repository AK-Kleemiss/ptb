# Lanthanide density fitting workflow

This directory contains the reproducible workflow for generating and fitting
La-Lu PTB density parameters against ORCA `wB97X-3c` density cubes.

The workflow is intentionally split into small steps:

1. Verify the current PTB lanthanide baseline.
2. Generate seed structures and a manifest.
3. Write ORCA inputs and SLURM scripts.
4. Run ORCA externally on the cluster.
5. Run PTB with `-denmat`, project the PTB density onto each ORCA cube grid,
   and compare real-space density residuals.
6. Optimize selected `.atompara` parameters with evolutionary/global search
   followed by local refinement.

## Quick start

```bash
python -m ptb_lnf.cli verify-baseline --repo-root ../..
python -m ptb_lnf.cli generate-structures --out work/structures --manifest work/manifest.jsonl
python -m ptb_lnf.cli write-orca --manifest work/manifest.jsonl --out work/orca
python -m ptb_lnf.cli write-slurm --orca-root work/orca --out work/orca_array.slurm
```

## MACE-Osaka structure generation

Create a dedicated venv, install CUDA-enabled PyTorch plus MACE/ASE, download
the official `mace-osaka26-small.model`, and generate MACE-relaxed La-Lu
training structures:

```bash
python setup_mace_osaka.py --per-element 4 --device cuda
```

For the intended production set of 150 geometries per element, use:

```bash
python setup_mace_osaka.py --training-set --device cuda
```

This produces 2250 relaxed structures across La-Lu. For chunked or exploratory
runs, keep using `--per-element N` with a smaller value; an explicit
`--per-element` value takes precedence over `--training-set`.

This creates:

- `.venv-mace-osaka/` with CUDA PyTorch, `mace-torch>=0.3.12`, ASE, NumPy, and SciPy
- `models/mace-osaka26-small.model`, downloaded from the official
  `qiqb-osaka/mace-osaka26` v0.0.1 release and SHA256-verified
- `work/mace_structures/` and `work/mace_manifest.jsonl`

The script defaults to the PyTorch CUDA 12.8 wheel index, which is appropriate
for recent NVIDIA drivers. Use `--cpu-torch --device cpu` only when GPU use is
not desired, or override `--torch-index-url` if your cluster/image requires a
different PyTorch CUDA wheel. MACE-Osaka26 is stored as a float64 model, so the
generator defaults to `--dtype float64`; use `--dtype float32` only for faster
screening runs where the extra conversion and lower precision are acceptable.
If an existing `.venv-mace-osaka/` contains a CPU-only torch build such as
`2.x+cpu`, the setup script now detects that after installing MACE dependencies,
reinstalls CUDA PyTorch from the configured wheel index, and stops before
generation if CUDA is still unavailable.

If PyTorch prints a `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD` warning while loading the
model, it is reporting that full model unpickling is enabled by the environment.
For this workflow the model is downloaded from the official MACE-Osaka release
and SHA256-verified before use, so that warning is expected rather than a failed
calculation.

The MACE generator is not a simple conformer rotator. It samples different
coordination numbers, ligand recipes, denticities, donor atoms, radial/angular
distortions, and mixed-ligand environments, then removes duplicates using a
rotation/translation-invariant pair-distance signature before and after MACE
relaxation. The quick default is 4 structures per element; `--training-set`
uses the production target of 150 per element. Use `--setup-only` to create the
environment/model without running generation, or
`python -m ptb_lnf.mace_generate --no-relax ...` only for fast generator tests.

After ORCA jobs finish, add `reference_cube` paths to the manifest or create a
derived manifest with those paths. PTB density projection then uses:

```bash
ptb molecule.xyz -par candidate.atompara -bas .basis_vDZP -denmat molecule.denmat
python -m ptb_lnf.cli denmat-to-cube --denmat molecule.denmat --grid reference.eldens.cube --out ptb.cube
python -m ptb_lnf.cli compare-cubes --reference reference.eldens.cube --candidate ptb.cube
```

## Notes

- MACE-Osaka is a conformer/relaxation source only. This workflow records the
  conformer source in the manifest, but reference densities are always ORCA
  density cubes.
- Generated structures are seed structures. They are meant to be relaxed or
  filtered before reference calculations.
- The optimizer harness requires SciPy. Core verification and cube comparison
  only require Python plus NumPy.
