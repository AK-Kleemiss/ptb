from __future__ import annotations

from pathlib import Path
import json
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.optimize import differential_evolution, minimize

from .cube import compare, read_cube
from .denmat import denmat_to_cube
from .orca import read_manifest


def _copy_candidate_atompara(template: Path, out: Path, scale: float) -> None:
    # Conservative v1 parameterization: global scale of mutable lanthanide numeric
    # rows. This is intentionally narrow until a reviewed mask is supplied.
    lines = template.read_text().splitlines()
    active = False
    result = []
    for line in lines:
        stripped = line.strip()
        if stripped.isdigit():
            active = 57 <= int(stripped) <= 71
            result.append(line)
            continue
        if active and stripped and not stripped.startswith("#"):
            vals = stripped.split()
            try:
                nums = [float(v) for v in vals]
            except ValueError:
                result.append(line)
                continue
            result.append(" ".join(f"{v * scale:15.10f}" if abs(v) > 1e-12 else f"{v:15.10f}" for v in nums))
        else:
            result.append(line)
    with out.open("w", newline="\n") as handle:
        handle.write("\n".join(result) + "\n")


def score_manifest(
    manifest: str | Path,
    ptb_exe: str | Path,
    atompara_template: str | Path,
    basis: str | Path,
    scale: float,
) -> float:
    rows = [row for row in read_manifest(manifest) if row.get("reference_cube") and row.get("split") == "train"]
    if not rows:
        raise ValueError("manifest has no train rows with reference_cube")
    with tempfile.TemporaryDirectory() as tmp_s:
        tmp = Path(tmp_s)
        candidate = tmp / "candidate.atompara"
        _copy_candidate_atompara(Path(atompara_template), candidate, scale)
        metrics = []
        failures = 0
        for row in rows:
            work = tmp / row["case_id"]
            work.mkdir()
            xyz = work / Path(row["xyz_path"]).name
            shutil.copy2(row["xyz_path"], xyz)
            denmat = work / "ptb.denmat"
            cube = work / "ptb.cube"
            cmd = [str(ptb_exe), str(xyz), "-par", str(candidate), "-bas", str(basis), "-denmat", str(denmat)]
            if row.get("charge") is not None:
                cmd.extend(["-chrg", str(row["charge"])])
            proc = subprocess.run(cmd, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if proc.returncode != 0 or not denmat.exists():
                failures += 1
                continue
            denmat_to_cube(denmat, row["reference_cube"], cube)
            metrics.append(compare(read_cube(row["reference_cube"]), read_cube(cube))["rmse"])
        if not metrics:
            return 1e6
        return float(np.mean(metrics) + failures * 10.0)


def optimize_scale(manifest: str | Path, ptb_exe: str | Path, atompara: str | Path, basis: str | Path) -> dict:
    def objective(x: np.ndarray) -> float:
        return score_manifest(manifest, ptb_exe, atompara, basis, float(x[0]))

    de = differential_evolution(objective, bounds=[(0.80, 1.20)], polish=False, workers=1)
    local = minimize(objective, x0=np.array([de.x[0]]), method="Nelder-Mead")
    return {
        "scale": float(local.x[0]),
        "score": float(local.fun),
        "global_score": float(de.fun),
    }
