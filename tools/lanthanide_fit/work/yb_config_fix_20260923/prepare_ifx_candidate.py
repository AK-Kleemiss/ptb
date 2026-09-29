"""Embed a hashed Yb atompara candidate in an isolated IFX source copy."""
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def yb_rows(path):
    lines = path.read_text(encoding="ascii").splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == "70") + 1
    stop = next(i for i, line in enumerate(lines) if line.strip() == "71")
    rows = lines[start:stop]
    if len(rows) != 15:
        raise ValueError("expected 15 Yb rows")
    return rows


parser = argparse.ArgumentParser()
parser.add_argument("--source", type=Path, required=True)
parser.add_argument("--candidate", type=Path, required=True)
parser.add_argument("--sha256", required=True)
args = parser.parse_args()
if sha(args.candidate) != args.sha256:
    raise ValueError("candidate SHA256 mismatch")
rows = yb_rows(args.candidate)
raw = args.source.read_bytes()
newline = b"\r\n" if b"\r\n" in raw else b"\n"
lines = raw.decode("ascii").splitlines()
header = next(i for i, line in enumerate(lines) if line.strip() == '"          70", &')
assert lines[header + 16].strip() == '"          71", &'
for row in rows:
    assert re.fullmatch(r"[ 0-9.+-]+", row), row
lines[header+1:header+16] = [f'      "{row}", &' for row in rows]
args.source.write_bytes((newline.decode().join(lines) + newline.decode()).encode("ascii"))
manifest = {"candidate": str(args.candidate.resolve()), "candidate_sha256": sha(args.candidate),
            "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
            "embedded_source_sha256": sha(args.source),
            "frozen_radial_rows_sha256": hashlib.sha256("\n".join(rows[r] for r in (6, 11)).encode()).hexdigest(),
            "frozen_response_13_10": rows[13].split()[10]}
(args.source.parent.parent / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(manifest, indent=2))
