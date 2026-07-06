#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import urllib.request
import venv


ROOT = Path(__file__).resolve().parent
DEFAULT_VENV = ROOT / ".venv-mace-osaka"
DEFAULT_MODEL = ROOT / "models" / "mace-osaka26-small.model"
MODEL_URL = "https://github.com/qiqb-osaka/mace-osaka26/releases/download/v0.0.1/mace-osaka26-small.model"
MODEL_SHA256 = "92865eee31cdf41b639213cd82e39e6ded86fe8690f1bb3da98506c63aff6878"
PILOT_PER_ELEMENT = 4
TRAINING_PER_ELEMENT = 150


def venv_python(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def run(cmd: list[str], cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.check_call(cmd, cwd=str(cwd) if cwd else None, env=env)


def torch_probe(py: Path) -> dict[str, object]:
    code = (
        "import json, torch; "
        "print(json.dumps({"
        "'version': torch.__version__, "
        "'cuda_available': torch.cuda.is_available(), "
        "'cuda_version': torch.version.cuda, "
        "'device_count': torch.cuda.device_count(), "
        "'device0': torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none'"
        "}))"
    )
    out = subprocess.check_output([str(py), "-c", code], text=True)
    return json.loads(out)


def print_torch_probe(probe: dict[str, object]) -> None:
    print(f"torch {probe['version']}")
    print(f"cuda_available {probe['cuda_available']}")
    print(f"cuda_version {probe['cuda_version']}")
    print(f"device_count {probe['device_count']}")
    print(f"device0 {probe['device0']}")


def install_cuda_torch(py: Path, torch_index_url: str | None, force: bool = False) -> None:
    cmd = [str(py), "-m", "pip", "install", "--upgrade"]
    if force:
        cmd.append("--force-reinstall")
    cmd.extend(["torch", "--index-url", torch_index_url or "https://download.pytorch.org/whl/cu128"])
    run(cmd)


def create_venv(venv_dir: Path, torch_index_url: str | None, cpu_torch: bool, require_cuda: bool) -> Path:
    if not venv_python(venv_dir).exists():
        print(f"Creating virtual environment: {venv_dir}")
        venv.EnvBuilder(with_pip=True, clear=False).create(venv_dir)
    py = venv_python(venv_dir)
    run([str(py), "-m", "pip", "install", "--upgrade", "pip", "wheel", "setuptools<82"])
    run([
        str(py), "-m", "pip", "install",
        "numpy>=1.24",
        "scipy>=1.10",
        "ase>=3.22",
        "mace-torch>=0.3.12",
    ])
    if cpu_torch:
        run([str(py), "-m", "pip", "install", "--upgrade", "torch"])
        return py

    probe = torch_probe(py)
    if not probe["cuda_available"] or probe["cuda_version"] is None:
        print("Installed torch is CPU-only or cannot see CUDA; reinstalling CUDA-enabled PyTorch.")
        install_cuda_torch(py, torch_index_url, force=True)
        probe = torch_probe(py)
    if require_cuda and (not probe["cuda_available"] or probe["cuda_version"] is None):
        raise SystemExit(
            "CUDA was requested, but the MACE venv still cannot use CUDA after installing PyTorch. "
            "Check the NVIDIA driver, or rerun with --cpu-torch --device cpu."
        )
    return py


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_model(model_path: Path, force: bool = False) -> None:
    model_path.parent.mkdir(parents=True, exist_ok=True)
    if model_path.exists() and not force:
        digest = sha256(model_path)
        if digest == MODEL_SHA256:
            print(f"Model already present and verified: {model_path}")
            return
        raise SystemExit(
            f"Existing model has unexpected sha256 {digest}. "
            f"Delete it or rerun with --force-model-download."
        )
    print(f"Downloading {MODEL_URL}")
    tmp = model_path.with_suffix(model_path.suffix + ".tmp")
    urllib.request.urlretrieve(MODEL_URL, tmp)
    digest = sha256(tmp)
    if digest != MODEL_SHA256:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"Downloaded model sha256 mismatch: {digest}")
    tmp.replace(model_path)
    print(f"Downloaded and verified: {model_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a MACE-Osaka venv, download the official Osaka26 model, and generate training structures."
    )
    parser.add_argument("--venv", type=Path, default=DEFAULT_VENV)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--out", type=Path, default=ROOT / "work" / "mace_structures")
    parser.add_argument("--manifest", type=Path, default=ROOT / "work" / "mace_manifest.jsonl")
    parser.add_argument("--per-element", type=int, default=None)
    parser.add_argument(
        "--training-set",
        action="store_true",
        help=f"Generate the production target of {TRAINING_PER_ELEMENT} structures per lanthanide.",
    )
    parser.add_argument("--seed", type=int, default=20260706)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--dtype", choices=("float64", "float32"), default="float64")
    parser.add_argument("--max-steps", type=int, default=120)
    parser.add_argument("--torch-index-url", default="https://download.pytorch.org/whl/cu128")
    parser.add_argument("--cpu-torch", action="store_true", help="Install CPU-only PyTorch instead of a CUDA wheel.")
    parser.add_argument("--force-model-download", action="store_true")
    parser.add_argument("--setup-only", action="store_true")
    parser.add_argument("--resume", action="store_true", help="Continue from an existing manifest or XYZ tree instead of starting over.")
    args = parser.parse_args()
    per_element = args.per_element if args.per_element is not None else (
        TRAINING_PER_ELEMENT if args.training_set else PILOT_PER_ELEMENT
    )

    require_cuda = args.device.startswith("cuda") and not args.cpu_torch
    py = create_venv(args.venv, args.torch_index_url, args.cpu_torch, require_cuda=require_cuda)
    print_torch_probe(torch_probe(py))
    download_model(args.model, force=args.force_model_download)
    if args.setup_only:
        print(f"Setup complete. Venv python: {py}")
        return

    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    run([
        str(py), "-m", "ptb_lnf.mace_generate",
        "--model", str(args.model),
        "--out", str(args.out),
        "--manifest", str(args.manifest),
        "--per-element", str(per_element),
        "--seed", str(args.seed),
        "--device", args.device,
        "--dtype", args.dtype,
        "--max-steps", str(args.max_steps),
    ] + (["--resume"] if args.resume else []), cwd=ROOT, env=env)


if __name__ == "__main__":
    main()
