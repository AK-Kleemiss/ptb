"""Screen one-coordinate g3->g2 backoffs against ions and neutral radius."""
import ast
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

os.environ.update(OMP_NUM_THREADS="14", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
import numpy as np
from scipy.io import FortranFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parents[1]))
from ptb_lnf.iongate import gate_violation, load_targets  # noqa: E402

helper = HERE.parent / "weekend_20260911/verify_cutoff_grid_20260912T181800Z.py"
tree = ast.parse(helper.read_text())
scope = {"FortranFile": FortranFile, "np": np}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef)
                              and n.name in ("grid", "metrics", "boundaries")], type_ignores=[]),
             str(helper), "exec"), scope)
grid, metrics, boundaries = (scope[name] for name in ("grid", "metrics", "boundaries"))

G2 = HERE / "ion_g2/final_candidate.atompara"
G3 = HERE / "ion_g3_size/final_candidate.atompara"
PTB = ROOT / "build/ion_g3_rc_20260929/build/ptb"
NOS = Path("D:/git/NoSpherA2/build/release-windows/bin/NoSpherA2.exe")
BAS = ROOT / ".basis_vDZP"
TARGETS = load_targets(HERE / "ion_targets_g1.json")
OUT = HERE / "ion_g3_backoffs_20260929"
OUT.mkdir(exist_ok=True)
RESULT = OUT / "result.json"
BASE_Q0 = 1.8656000588851611
LIMIT_Q0 = 1.01 * BASE_Q0


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_process(args, cwd):
    done = subprocess.run(args, cwd=cwd, capture_output=True, text=True, timeout=180)
    if done.returncode:
        raise RuntimeError(f"{args[0]} failed: {done.stdout[-300:]} {done.stderr[-300:]}")


def neutral_radius(candidate, work):
    xyz = work / "atom.xyz"
    xyz.write_text("1\nbackoff radius\nYb 0 0 0\n")
    wfn = work / "wfn.xtb"
    check_process([str(PTB), str(xyz), "-par", str(candidate), "-bas", str(BAS),
                   "-chrg", "0", "-uhf", "0", "-stda"], work)
    edges, screen = boundaries(wfn)
    assert abs(screen["export_occupation_sum"] - 24) < 1e-6
    coords, weights, radii = grid(edges, 32)
    points = work / "export32.points"
    with points.open("w") as handle:
        handle.write(f"{len(coords)}\n")
        np.savetxt(handle, coords, fmt="%.16e")
    check_process([str(NOS), "-no_gpu", "-cpus", "14", "-wfn", str(wfn),
                   "-rho_at_points", str(points)], work)
    value = metrics(np.loadtxt(str(points) + ".rho"), weights, radii)
    assert all(abs(count - 24) < 1e-3 for count in value["counts"].values())
    return value["rms_radius"]


assert sha(G2) == "6eb12794d40838ce9b96050df769be7daa1c6e311d941b5dd6b048864ea51c39"
assert sha(G3) == "9a2109cadc3e463bc1d766c0f0a6abc3b30da908b91d5b25fe49c6dfa8aaf0df"
a, b = G2.read_text().splitlines(), G3.read_text().splitlines()
start = next(i for i, line in enumerate(a) if line.strip() == "70") + 1
positions = []
for row in range(6, 15):
    for col, (g2, g3) in enumerate(zip(a[start+row].split(), b[start+row].split())):
        if g2 != g3:
            positions.append((row, col, g2, g3))

def main():
    records = []
    with tempfile.TemporaryDirectory(prefix="yb_backoff_") as temp:
        base = Path(temp)
        for index, (row, col, g2, g3) in enumerate(positions, 1):
            work = base / f"r{row}c{col}"
            work.mkdir()
            lines = b.copy()
            values = lines[start+row].split()
            values[col] = g2
            lines[start+row] = " ".join(values)
            candidate = work / "candidate.atompara"
            candidate.write_text("\n".join(lines) + "\n", encoding="ascii")
            record = {"row": row, "col": col, "g2_value": g2, "g3_value": g3,
                      "candidate_sha256": sha(candidate)}
            try:
                violation, populations = gate_violation(PTB, candidate, BAS, TARGETS)
                record["violation"] = violation
                record["populations"] = populations
                if violation == 0.0:
                    record["neutral_radius"] = neutral_radius(candidate, work)
                    record["passes_radius"] = record["neutral_radius"] <= LIMIT_Q0
                    if record["passes_radius"]:
                        (OUT / f"r{row}c{col}.atompara").write_bytes(candidate.read_bytes())
            except Exception as exc:
                record["error"] = repr(exc)
            records.append(record)
            RESULT.write_text(json.dumps({"g2_sha256": sha(G2), "g3_sha256": sha(G3),
                                          "ptb_sha256": sha(PTB), "basis_sha256": sha(BAS),
                                          "q0_radius_limit": LIMIT_Q0, "completed": len(records),
                                          "total": len(positions), "records": records}, indent=2) + "\n")
            if index % 10 == 0 or record.get("passes_radius"):
                print(index, len(positions), row, col, record.get("violation"),
                      record.get("neutral_radius"), record.get("passes_radius"), flush=True)


if __name__ == "__main__":
    main()
