"""Free-ion scan of single Yb .atompara values: python scan_xi.py ROW COL v1 v2 ..."""
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from ionprobe import fmt, probe  # noqa: E402

ROOT = HERE.parents[3]
PTB = ROOT / "build" / "ptb_windows_ifx_yb_basis_d2e_20260923.exe"
BAS = ROOT / ".basis_vDZP"
BASE = HERE.parent / "weekend_20260911" / "release_yb_basis_d2e_20260923" / "atompara_ifx_yb_basis_d2e_20260923.atompara"


def set_values(src: Path, dst: Path, edits: dict[tuple[int, int], float], z: int = 70) -> None:
    lines = src.read_text().splitlines()
    h = next(i for i, l in enumerate(lines) if l.strip() == str(z))
    for (r, c), v in edits.items():
        vals = lines[h + 1 + r].split()
        vals[c] = f"{v:.10f}"
        lines[h + 1 + r] = " ".join(vals)
    dst.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    args = sys.argv[1:]
    base = Path(args.pop(0)) if args and args[0].endswith(".atompara") else BASE
    # args: groups "r,c=v;r,c=v" per candidate
    with tempfile.TemporaryDirectory() as d:
        cand = Path(d) / "c.atompara"
        for spec in args:
            edits = {}
            for item in filter(None, spec.split(";")):
                pos, v = item.split("=")
                r, c = map(int, pos.split(","))
                edits[(r, c)] = float(v)
            set_values(base, cand, edits)
            print(f"{spec:28s}", fmt(probe(PTB, cand, BAS)), flush=True)
