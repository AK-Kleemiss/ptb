"""Create a training-path candidate between released ion_g2 and size-failed ion_g3."""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
G2 = HERE / "ion_g2/final_candidate.atompara"
G3 = HERE / "ion_g3_size/final_candidate.atompara"
OUT = HERE / "ion_g3_path_20260929"
OUT.mkdir(exist_ok=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


assert sha(G2) == "6eb12794d40838ce9b96050df769be7daa1c6e311d941b5dd6b048864ea51c39"
assert sha(G3) == "9a2109cadc3e463bc1d766c0f0a6abc3b30da908b91d5b25fe49c6dfa8aaf0df"
a, b = G2.read_text().splitlines(), G3.read_text().splitlines()
assert len(a) == len(b)
start = next(i for i, line in enumerate(a) if line.strip() == "70") + 1
stop = next(i for i, line in enumerate(a) if line.strip() == "71")
assert stop-start == 15
assert all(a[i] == b[i] for i in range(len(a)) if i < start or i >= stop)

parser = argparse.ArgumentParser()
parser.add_argument("fractions", nargs="+", type=float)
fractions = parser.parse_args().fractions
for fraction in fractions:
    if not 0 < fraction < 1:
        raise ValueError("fraction must be strictly between ion_g2 and ion_g3")
    rows = a.copy()
    for i in range(start, stop):
        aa, bb = [float(x) for x in a[i].split()], [float(x) for x in b[i].split()]
        assert len(aa) == len(bb)
        rows[i] = " ".join(f"{u + fraction*(v-u):.10f}" for u,v in zip(aa,bb))
    assert all(a[i].split() == rows[i].split() for i in (start+6, start+11))
    assert a[start+13].split()[10] == rows[start+13].split()[10]
    name = f"path_{fraction:.4f}"
    target = OUT / f"{name}.atompara"
    target.write_text("\n".join(rows) + "\n", encoding="ascii")
    record = {"fraction": fraction, "candidate_sha256": sha(target),
              "ion_g2_sha256": sha(G2), "ion_g3_sha256": sha(G3),
              "frozen_radial_rows": True, "path": str(target)}
    (OUT / f"{name}.json").write_text(json.dumps(record, indent=2) + "\n")
    print(name, sha(target))
