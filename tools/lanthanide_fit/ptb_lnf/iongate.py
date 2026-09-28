"""Free-ion configuration gate for lanthanide fits.

The molecular density objective cannot see a wrong 4f/5d split (the reference
cubes are too coarse to tell 4f from 5d), and the count gate cannot either.
The 23 Sep 2026 Yb d2e fit landed in a d10 f1.7 "Yb3+". This gate runs the
free ion at each charge, run like the molecular scoring (no -stda), and
reports how far its l-populations are from the target windows.

On a single atom the l blocks of S decouple, so the block traces of P.S are
exact l-populations, not a partitioning choice.

Targets JSON: {"element": "Yb", "charges": {"3": {"f": [12.5, 13.5], "d": [0, 0.5]}, ...}}
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from .denmat import read_denmat
from .elements import SYMBOL_TO_Z

# Valence electrons of the neutral Ln atom in PTB's vDZP/ECP partition.
VALENCE = {z: z - 46 for z in range(57, 72)}


def ion_populations(ptb: Path, atompara: Path, basis: Path, z: int, charge: int) -> dict[str, float]:
    ptb, atompara, basis = (Path(p).resolve() for p in (ptb, atompara, basis))
    for p in (ptb, atompara, basis):
        if not p.is_file():
            # PTB silently falls back to its built-in files for a missing -par/-bas.
            raise FileNotFoundError(p)
    n = VALENCE[z] - charge
    symbol = next(s for s, zz in SYMBOL_TO_Z.items() if zz == z)
    with tempfile.TemporaryDirectory(prefix="iongate") as d:
        work = Path(d)
        (work / "atom.xyz").write_text(f"1\nion\n{symbol} 0 0 0\n")
        cmd = [str(ptb), "atom.xyz", "-chrg", str(charge), "-uhf", str(n % 2), "-denmat", "ptb.denmat",
               "-par", str(atompara), "-bas", str(basis)]
        done = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=300)
        if done.returncode or not (work / "ptb.denmat").exists():
            raise RuntimeError(f"free-ion PTB run failed ({symbol} q={charge}): {done.stdout[-400:]}")
        dm = read_denmat(work / "ptb.denmat")
    ps = np.diag(dm.density @ dm.overlap)
    pops = {l: float(sum(ps[i] for i, ao in enumerate(dm.aos) if ao.angular_l == k)) for k, l in enumerate("spdf")}
    pops["total"] = float(ps.sum())
    return pops


def gate_violation(ptb, atompara, basis, targets: dict) -> tuple[float, dict]:
    """Sum over charges and l of the distance outside each [lo, hi] window
    (electrons), plus the per-charge populations for logging."""
    z = SYMBOL_TO_Z[targets["element"]]
    violation, report = 0.0, {}
    for q, windows in targets["charges"].items():
        try:
            pops = ion_populations(ptb, atompara, basis, z, int(q))
        except RuntimeError:
            return 100.0, report  # an SCF that fails on the bare ion is as bad as a wrong configuration
        report[q] = pops
        for l, (lo, hi) in windows.items():
            violation += max(0.0, lo - pops[l]) + max(0.0, pops[l] - hi)
    return violation, report


def load_targets(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())


if __name__ == "__main__":
    # Self-check of the window arithmetic on a fake report.
    t = {"element": "Yb", "charges": {"3": {"f": [12.5, 13.5], "d": [0.0, 0.5]}}}
    fake = {"f": 11.0, "d": 1.0}
    v = sum(max(0.0, lo - fake[l]) + max(0.0, fake[l] - hi) for l, (lo, hi) in t["charges"]["3"].items())
    assert abs(v - 2.0) < 1e-12, v
    print("ok")
