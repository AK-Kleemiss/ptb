# Lanthanide Fit Handover

Snapshot date: 2026-07-06

This file is for continuing the MACE-Osaka lanthanide structure-generation
campaign on another Windows desktop or workstation with a similar NVIDIA GPU
setup.

## Current Snapshot

- Repository: `D:\git\ptb`
- Commit at snapshot: `c015f8d`
- Target: 150 MACE-relaxed geometries per lanthanide, La-Lu
- Target total: 2250 geometries
- Current recovered total: 409 geometries
- Current recovered manifest: `tools/lanthanide_fit/work/mace_manifest.recovered.jsonl`
- CUDA smoke status on source machine: `torch 2.11.0+cu128`, CUDA 12.8, `NVIDIA GeForce RTX 3090`

Current counts recovered from `tools/lanthanide_fit/work/mace_structures`:

| Element | Count |
| --- | ---: |
| La | 150 |
| Ce | 154 |
| Pr | 85 |
| Nd | 4 |
| Pm | 4 |
| Sm | 4 |
| Eu | 4 |
| Gd | 4 |

No structures were present yet for Tb, Dy, Ho, Er, Tm, Yb, or Lu at the time of
this snapshot.

## Files To Transfer

Copy the repository checkout plus the generated work tree:

```powershell
D:\git\ptb
```

Important paths inside the checkout:

- `tools\lanthanide_fit\setup_mace_osaka.py`
- `tools\lanthanide_fit\ptb_lnf\mace_generate.py`
- `tools\lanthanide_fit\models\mace-osaka26-small.model`
- `tools\lanthanide_fit\work\mace_structures\`
- `tools\lanthanide_fit\work\mace_manifest.recovered.jsonl`

Do not rely on copying `.venv-mace-osaka`. Recreate it on the new machine so
PyTorch matches the local Python, CUDA driver, and GPU.

## New Machine Setup

From a PowerShell terminal on the new computer:

```powershell
cd D:\git\ptb
python tools\lanthanide_fit\setup_mace_osaka.py --setup-only --device cuda
```

The setup command should end with output like:

```text
torch 2.11.0+cu128
cuda_available True
cuda_version 12.8
device_count 1
device0 <NVIDIA GPU name>
```

If `cuda_available` is `False`, check the NVIDIA driver first. As a fallback for
CPU-only continuation, use `--cpu-torch --device cpu`, but expect much slower
generation.

## Continue The Current Campaign

Use `--resume` so existing XYZ files are scanned and reused. The generator will
only create missing structures up to the 150-per-element target.

```powershell
cd D:\git\ptb
python tools\lanthanide_fit\setup_mace_osaka.py --training-set --device cuda --resume
```

If you want to rebuild the manifest from the transferred XYZ tree before
continuing, run:

```powershell
cd D:\git\ptb\tools\lanthanide_fit
python -m ptb_lnf.mace_generate --model models\mace-osaka26-small.model --out work\mace_structures --manifest work\mace_manifest.recovered.jsonl --per-element 0 --no-relax --resume
```

Then continue with the production command above. The production command writes
the default manifest at `tools\lanthanide_fit\work\mace_manifest.jsonl`.

## Expected Warnings

These warnings are non-fatal:

- `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD`: PyTorch is allowing full model unpickling.
  The workflow uses the official MACE-Osaka model and verifies the SHA256 hash.
- `cuequivariance or cuequivariance_torch is not available`: optional
  acceleration is missing; MACE still runs.

Stop and investigate if the output says `cuda_available False` while running
with `--device cuda`.

## Verification Commands

Check the venv GPU status:

```powershell
tools\lanthanide_fit\.venv-mace-osaka\Scripts\python.exe -c "import torch; print(torch.__version__); print(torch.cuda.is_available()); print(torch.version.cuda); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')"
```

Check current structure counts:

```powershell
Get-ChildItem tools\lanthanide_fit\work\mace_structures -Directory | ForEach-Object {
  $n = (Get-ChildItem $_.FullName -Recurse -Filter *.xyz | Measure-Object).Count
  [pscustomobject]@{Element=$_.Name; Count=$n}
} | Sort-Object Element | Format-Table -AutoSize
```

Run the local tests:

```powershell
python -m pytest tools\lanthanide_fit\tests
```
