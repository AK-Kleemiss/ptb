"""Check embedded ion_g2 defaults against the candidate and unchanged controls."""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parents[1]))
from ptb_lnf.denmat import read_denmat  # noqa: E402

PTB = Path(os.environ.get("YB_RELEASE_PTB", str(ROOT / "build/ptb_windows_ifx_yb_ion_g2_20260928.exe")))
OLD = Path(os.environ.get("YB_RELEASE_OLD", str(ROOT / "build/ptb_windows_ifx_yb_d7_joint2_20260921.exe")))
PAR = Path(os.environ.get("YB_RELEASE_PAR", str(HERE / "ion_g2/final_candidate.atompara")))
BAS = ROOT / ".basis_vDZP"
OUT = Path(os.environ.get("YB_RELEASE_OUT", str(HERE / "ion_g2/embedded_validation.json")))
WORK = Path(os.environ.get("YB_RELEASE_WORK", str(HERE / "ion_g2/embedded_smoke")))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(exe: Path, work: Path, element: str, charge: int, uhf: int, extra: list[str]):
    work.mkdir(parents=True, exist_ok=True)
    xyz = work / "atom.xyz"
    xyz.write_text(f"1\nembedded validation\n{element} 0 0 0\n", encoding="ascii")
    dm = work / "ptb.denmat"
    done = subprocess.run([str(exe), str(xyz), "-chrg", str(charge), "-uhf", str(uhf),
                           "-stda", "-denmat", str(dm), *extra], cwd=work,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (work / "ptb.log").write_text(done.stdout, encoding="utf-8")
    if done.returncode or not dm.exists():
        raise RuntimeError(f"{exe.name} {element} q{charge} failed")
    data = read_denmat(dm)
    return data.density, float(np.trace(data.density @ data.overlap))


results = []
for element, z in (("Yb", 70), ("Tm", 69), ("La", 57)):
    for charge in (0, 2, 3):
        expected = z - 46 - charge
        work = WORK / element
        pe, te = run(PTB, work / f"embedded_q{charge}", element, charge, expected % 2, [])
        pp, tp = run(PTB, work / f"parfile_q{charge}", element, charge, expected % 2,
                     ["-par", str(PAR), "-bas", str(BAS)])
        po, to = run(OLD, work / f"old_q{charge}", element, charge, expected % 2, [])
        record = {"element": element, "charge": charge, "expected": expected,
                  "trace_embedded": te, "trace_parfile": tp, "trace_old": to,
                  "maxabs_embedded_vs_parfile": float(np.abs(pe - pp).max()),
                  "maxabs_new_vs_old": float(np.abs(pe - po).max())}
        assert abs(te - expected) < 1e-6 and record["maxabs_embedded_vs_parfile"] == 0, record
        if element != "Yb":
            assert record["maxabs_new_vs_old"] == 0, record
        results.append(record)
        print(record, flush=True)

OUT.write_text(json.dumps({"ptb_sha256": sha(PTB), "old_sha256": sha(OLD),
                           "candidate_sha256": sha(PAR), "basis_sha256": sha(BAS),
                           "results": results, "status": "passed"}, indent=2) + "\n")
print("PASSED")
