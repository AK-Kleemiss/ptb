"""Compare late-Ln AO shell extents on a common atomic grid."""
import hashlib
import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

os.environ.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parents[1]))
from ptb_lnf.denmat import evaluate_density, read_denmat  # noqa: E402

PTB = Path(os.environ.get("YB_SHELL_PTB", str(ROOT / "build/ptb_windows_ifx_yb_ion_g2_20260928.exe")))
PAR = Path(os.environ.get("YB_SHELL_CANDIDATE", str(HERE / "ion_g2/final_candidate.atompara")))
BAS = ROOT / ".basis_vDZP"
OUT = Path(os.environ.get("YB_SHELL_OUT", str(HERE / "shell_sizes_20260928")))
OUT.mkdir(exist_ok=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grid(order):
    # Angular rule integrates products through f exactly; radial segments
    # separate the core, valence, and far tail without fitting to a candidate.
    t, wt = np.polynomial.legendre.leggauss(8)
    az = np.arange(16) * np.pi / 8
    dirs = np.array([[np.sqrt(1 - z*z)*np.cos(a), np.sqrt(1 - z*z)*np.sin(a), z]
                     for z in t for a in az])
    wa = np.repeat(wt, 16) * np.pi / 8
    x, wx = np.polynomial.legendre.leggauss(order)
    edges = [0, .125, .5, 2, 6, 20, 80]
    r = np.concatenate([lo + (x + 1)*(hi - lo)/2 for lo, hi in zip(edges[:-1], edges[1:])])
    wr = np.concatenate([wx*(hi - lo)/2*(lo + (x + 1)*(hi - lo)/2)**2
                         for lo, hi in zip(edges[:-1], edges[1:])])
    return (r[:, None, None]*dirs).reshape(-1, 3), (wr[:, None]*wa).ravel(), np.repeat(r, len(wa))


def measure(dm, indices, coords, weights, radii):
    if not indices:
        return None
    aos = [dm.aos[i] for i in indices]
    p = np.eye(len(indices))
    subset = replace(dm, aos=aos, density=p, overlap=dm.overlap[np.ix_(indices, indices)])
    rho = evaluate_density(subset, coords)
    n = float(np.dot(rho, weights))
    return {"integral": n, "rms_radius": float(np.sqrt(np.dot(rho*weights, radii*radii)/n)) if n > 1e-9 else None,
            "outside_20": float(np.dot(rho[radii >= 20], weights[radii >= 20]))}


records = []
for element, z in (("Tm", 69), ("Yb", 70), ("Lu", 71)):
    for charge in (0, 2, 3):
        work = OUT / f"{element}_q{charge}"
        work.mkdir(exist_ok=True)
        xyz = work / "atom.xyz"
        xyz.write_text(f"1\nlate Ln shell size\n{element} 0 0 0\n")
        dmfile = work / "ptb.denmat"
        done = subprocess.run([str(PTB), str(xyz), "-par", str(PAR), "-bas", str(BAS),
                               "-chrg", str(charge), "-uhf", str((z-46-charge)%2),
                               "-denmat", str(dmfile)], cwd=work, capture_output=True, text=True)
        (work / "ptb.log").write_text(done.stdout + done.stderr)
        if done.returncode or not dmfile.exists():
            raise RuntimeError(f"{element} q{charge}: PTB failed")
        dm = read_denmat(dmfile)
        expected = z - 46 - charge
        trace = float(np.trace(dm.density @ dm.overlap))
        assert abs(trace - expected) < 1e-6, (element, charge, trace)
        by_shell = {}
        for shell in (6, 7, 8, 9):
            indices = [i for i, ao in enumerate(dm.aos) if ao.shell_index == shell]
            if not indices:
                continue
            outputs = []
            for order in (32, 64):
                coords, weights, radii = grid(order)
                outputs.append({"order": order,
                                "ao_extent": measure(dm, indices, coords, weights, radii)})
            a, b = (outputs[i]["ao_extent"] for i in (0, 1))
            assert abs(a["rms_radius"] - b["rms_radius"]) < 1e-4, (element, charge, shell)
            by_shell[str(shell)] = outputs
        records.append({"element": element, "charge": charge, "trace": trace,
                        "shells": by_shell})
        print(element, charge, {s: round(v[1]["ao_extent"]["rms_radius"], 3) for s, v in by_shell.items()}, flush=True)

(OUT / "result.json").write_text(json.dumps({"ptb_sha256": sha(PTB), "parameter_sha256": sha(PAR),
                                               "basis_sha256": sha(BAS), "records": records}, indent=2) + "\n")
