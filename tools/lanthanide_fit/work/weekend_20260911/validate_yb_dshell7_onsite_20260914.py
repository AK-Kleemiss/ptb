"""Independent Windows PTB/NoSpherA2 validation of the staged Yb d-shell-7 candidate.

This is intentionally a standalone gate.  It uses the actual NoSpherA2
wavefunction reader and its density-at-points mode; PTB is used only to
export the wavefunction and to provide a density-matrix trace check.
"""
import ast
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

os.environ.update(OMP_NUM_THREADS="14", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
HERE = Path(__file__).resolve().parent
OUT = Path(os.environ.get("YB_VALIDATION_OUT", str(HERE / "yb_dshell7_onsite_validation_20260914")))
OUT.mkdir(exist_ok=False)
sys.path.insert(0, str(HERE.parents[1]))
from ptb_lnf.denmat import read_denmat  # noqa: E402

HELPER = HERE / "verify_cutoff_grid_20260912T181800Z.py"
tree = ast.parse(HELPER.read_text(encoding="utf-8"))
scope: dict[str, object] = {"FortranFile": FortranFile, "np": np}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ("grid", "metrics", "boundaries")], type_ignores=[]), str(HELPER), "exec"), scope)
grid, metrics, boundaries = scope["grid"], scope["metrics"], scope["boundaries"]
PTB = Path(os.environ.get("YB_PTB", "D:/git/ptb/build/ptb_windows_ifx_physical_refine_20260914.exe"))
NOS = Path("D:/git/NoSpherA2/build/release-windows/bin/NoSpherA2.exe")
PARAM = Path(os.environ.get("YB_CANDIDATE", str(HERE / "yb_dshell7_onsite_candidate_20260914" / "candidate.atompara")))
BASIS = Path("D:/git/ptb/build/ifx_basis_vDZP")

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def call(args: list[str], cwd: Path, log: Path) -> None:
    done = subprocess.run(args, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    log.write_text(done.stdout, encoding="utf-8")
    if done.returncode:
        raise RuntimeError(f"command failed {done.returncode}: {args[0]}")

def main() -> None:
    if not PTB.exists() or not NOS.exists() or not PARAM.exists() or not BASIS.exists():
        raise FileNotFoundError("required PTB, NoSpherA2, candidate, or basis path is missing")
    states = []
    for charge in (0, 2, 3):
        expected = 24 - charge
        work = OUT / f"q{charge}"
        work.mkdir()
        xyz = work / "atom.xyz"
        xyz.write_text("1\nIndependent Yb onsite validation\nYb 0 0 0\n", encoding="utf-8")
        wfn = work / "wfn.xtb"
        call([str(PTB), str(xyz), "-par", str(PARAM), "-bas", str(BASIS), "-chrg", str(charge), "-uhf", str(expected % 2), "-stda", "-denmat", str(work / "ptb.denmat")], work, work / "ptb.log")
        dm = read_denmat(work / "ptb.denmat")
        trace = float(np.trace(dm.density @ dm.overlap))
        if abs(trace - expected) > 1e-6:
            raise RuntimeError(f"q{charge}: trace {trace}, expected {expected}")
        edges, screen = boundaries(wfn)
        if abs(screen["export_occupation_sum"] - expected) > 1e-6:
            raise RuntimeError(f"q{charge}: exported occupation mismatch")
        exports = []
        for order in (32, 64):
            coords, weights, radii = grid(edges, order)
            points = work / f"export{order}.points"
            with points.open("w", encoding="utf-8") as handle:
                handle.write(f"{len(coords)}\n")
                np.savetxt(handle, coords, fmt="%.16e")
            call([str(NOS), "-no_gpu", "-cpus", "14", "-wfn", str(wfn), "-rho_at_points", str(points)], work, work / f"nos{order}.log")
            metric = metrics(np.loadtxt(str(points) + ".rho"), weights, radii)
            if any(abs(value - expected) > 1e-3 for value in metric["counts"].values()):
                raise RuntimeError(f"q{charge}: independent count gate")
            exports.append({"order": order, **metric})
        if max(abs(exports[0]["counts"][key] - exports[1]["counts"][key]) for key in ("20", "80")) >= 1e-5:
            raise RuntimeError(f"q{charge}: quadrature convergence gate")
        states.append({"charge": charge, "expected": expected, "trace": trace, "independent_export": exports})
        result = {"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "ptb_sha256": sha(PTB), "nosphera2_sha256": sha(NOS), "candidate_sha256": sha(PARAM), "basis_sha256": sha(BASIS), "states": states, "release": False}
        temporary = OUT / "result.new"
        temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        temporary.replace(OUT / "result.json")
        print(f"Yb q{charge} passed", flush=True)

if __name__ == "__main__":
    main()
