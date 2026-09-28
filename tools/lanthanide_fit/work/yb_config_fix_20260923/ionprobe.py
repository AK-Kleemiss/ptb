"""Free-ion probe: l-resolved populations and valence RMS radius of an Ln atom
for charges 0/2/3, run exactly like the molecular scoring (no -stda).

On a single atom the l blocks of S are decoupled, so the block traces of P.S
are exact l-populations, not a partitioning choice."""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from ptb_lnf.denmat import evaluate_density, read_denmat  # noqa: E402

VALENCE = {57: 11, 58: 12, 59: 13, 60: 14, 61: 15, 62: 16, 63: 17, 64: 18, 65: 19,
           66: 20, 67: 21, 68: 22, 69: 23, 70: 24, 71: 25}
SYMBOL = {57: "La", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu"}


def spherical_rule(order: int = 24):
    u, wu = np.polynomial.legendre.leggauss(10)
    phi = np.arange(20) * (2 * np.pi / 20)
    dirs = np.array([[np.sqrt(1 - t * t) * np.cos(p), np.sqrt(1 - t * t) * np.sin(p), t] for t in u for p in phi])
    ang = np.repeat(wu, len(phi)) * (2 * np.pi / len(phi))
    x, w = np.polynomial.legendre.leggauss(order)
    radii, weights = [], []
    for lo, hi in zip([0, .25, 1, 3, 8, 20], [.25, 1, 3, 8, 20, 80]):
        r = lo + (x + 1) * (hi - lo) / 2
        radii.extend(r)
        weights.extend(w * (hi - lo) / 2 * r * r)
    r = np.asarray(radii)
    coords = (r[:, None, None] * dirs[None, :, :]).reshape(-1, 3)
    return coords, (np.asarray(weights)[:, None] * ang).ravel(), np.repeat(r, len(ang))


RULE = spherical_rule()


def probe(ptb: Path, atompara: Path | None, basis: Path | None, z: int = 70, charges=(0, 2, 3)) -> list[dict]:
    # PTB silently falls back to its built-in files when -par/-bas is missing.
    ptb, atompara, basis = (None if p is None else Path(p).resolve() for p in (ptb, atompara, basis))
    for p in (ptb, atompara, basis):
        if p is not None and not p.is_file():
            raise FileNotFoundError(p)
    out = []
    for q in charges:
        n = VALENCE[z] - q
        with tempfile.TemporaryDirectory(prefix="ionprobe") as d:
            work = Path(d)
            (work / "atom.xyz").write_text(f"1\nprobe\n{SYMBOL[z]} 0 0 0\n")
            cmd = [str(ptb), "atom.xyz", "-chrg", str(q), "-uhf", str(n % 2), "-denmat", "ptb.denmat"]
            if atompara is not None:
                cmd += ["-par", str(atompara)]
            if basis is not None:
                cmd += ["-bas", str(basis)]
            done = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=120)
            if done.returncode or not (work / "ptb.denmat").exists():
                raise RuntimeError(done.stdout[-400:])
            dm = read_denmat(work / "ptb.denmat")
        ps = np.diag(dm.density @ dm.overlap)
        pop = {l: float(sum(ps[i] for i, ao in enumerate(dm.aos) if ao.angular_l == k)) for k, l in enumerate("spdf")}
        coords, w, r = RULE
        rho = evaluate_density(dm, coords) * w
        tot = float(rho.sum())
        shells = {}
        for i, ao in enumerate(dm.aos):
            shells[ao.shell_index] = shells.get(ao.shell_index, 0.0) + float(ps[i])
        out.append(dict(charge=q, expected=n, trace=float(ps.sum()), integral=tot, shells=shells,
                        rms=float(np.sqrt((rho * r * r).sum() / tot)), **pop))
    return out


def fmt(states: list[dict]) -> str:
    return " | ".join(f"q{s['charge']}: d {s['d']:5.2f} f {s['f']:5.2f} rms {s['rms']:.2f}" for s in states)


if __name__ == "__main__":
    ptb = Path(sys.argv[1])
    par = Path(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != "-" else None
    bas = Path(sys.argv[3]) if len(sys.argv) > 3 else None
    zs = [int(a) for a in sys.argv[4:]] or [70]
    for z in zs:
        print(SYMBOL[z], fmt(probe(ptb, par, bas, z)))


def fmt_shells(states: list[dict]) -> str:
    return "\n".join(f"  q{s['charge']}: " + " ".join(f"{k}:{v:6.3f}" for k, v in sorted(s["shells"].items())) for s in states)
