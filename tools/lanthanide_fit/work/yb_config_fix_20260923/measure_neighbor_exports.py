"""Independent NoSpherA2 atom-size controls for Tm and Lu, all Yb charge states."""
import ast
import hashlib
import json
import os
import subprocess
from pathlib import Path

os.environ.update(OMP_NUM_THREADS="14", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
import numpy as np
from scipy.io import FortranFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
HELPER = HERE.parent / "weekend_20260911/verify_cutoff_grid_20260912T181800Z.py"
tree = ast.parse(HELPER.read_text())
scope = {"FortranFile": FortranFile, "np": np}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef)
                              and n.name in ("grid", "metrics", "boundaries")], type_ignores=[]),
             str(HELPER), "exec"), scope)
grid, metrics, boundaries = (scope[name] for name in ("grid", "metrics", "boundaries"))
PTB = ROOT / "build/ptb_windows_ifx_yb_ion_g2_20260928.exe"
NOS = Path("D:/git/NoSpherA2/build/release-windows/bin/NoSpherA2.exe")
PAR = HERE / "ion_g2/final_candidate.atompara"
BAS = ROOT / ".basis_vDZP"
OUT = HERE / "neighbor_exports_20260928"
OUT.mkdir(exist_ok=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


records = []
for element, z in (("Tm", 69), ("Lu", 71)):
    for charge in (0, 2, 3):
        n = z - 46 - charge
        work = OUT / f"{element}_q{charge}"
        work.mkdir(exist_ok=True)
        xyz = work / "atom.xyz"
        xyz.write_text(f"1\nneighbor size control\n{element} 0 0 0\n")
        done = subprocess.run([str(PTB), str(xyz), "-par", str(PAR), "-bas", str(BAS),
                               "-chrg", str(charge), "-uhf", str(n % 2), "-stda"],
                              cwd=work, capture_output=True, text=True)
        (work / "ptb.log").write_text(done.stdout + done.stderr)
        if done.returncode:
            raise RuntimeError(f"PTB {element} q{charge}")
        wfn = work / "wfn.xtb"
        edges, screen = boundaries(wfn)
        assert abs(screen["export_occupation_sum"] - n) < 1e-6
        runs = []
        for order in (32, 64):
            coords, weights, radii = grid(edges, order)
            points = work / f"export{order}.points"
            with points.open("w") as handle:
                handle.write(f"{len(coords)}\n")
                np.savetxt(handle, coords, fmt="%.16e")
            done = subprocess.run([str(NOS), "-no_gpu", "-cpus", "14", "-wfn", str(wfn),
                                   "-rho_at_points", str(points)], cwd=work, capture_output=True, text=True)
            (work / f"nos{order}.log").write_text(done.stdout + done.stderr)
            if done.returncode:
                raise RuntimeError(f"NoSpherA2 {element} q{charge} order{order}")
            metric = metrics(np.loadtxt(str(points) + ".rho"), weights, radii)
            assert all(abs(value - n) < 1e-3 for value in metric["counts"].values())
            runs.append({"order": order, **metric})
        assert abs(runs[0]["rms_radius"] - runs[1]["rms_radius"]) < 1e-5
        records.append({"element": element, "charge": charge, "expected": n, "runs": runs})
        print(element, charge, runs[1]["rms_radius"], flush=True)

(OUT / "result.json").write_text(json.dumps({"ptb_sha256": sha(PTB), "nosphera2_sha256": sha(NOS),
                                               "parameter_sha256": sha(PAR), "basis_sha256": sha(BAS),
                                               "records": records}, indent=2) + "\n")
