"""Embed the independently validated, size-preserving g4 Yb defaults."""
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "source/default_files.f90"
PAR = Path(__file__).resolve().parent / "ion_g4_sizepass/final_candidate.atompara"
EXPECTED_SHA = "aef910c4fb85d7d668cba2c2574661f7c68cbd8cad8df8c3763c3bca3643ddf8"
FROZEN = (6, 7, 11, 13, 14)

assert hashlib.sha256(PAR.read_bytes()).hexdigest() == EXPECTED_SHA
par_lines = PAR.read_text(encoding="ascii").splitlines()
start = next(i for i, line in enumerate(par_lines) if line.strip() == "70") + 1
stop = next(i for i, line in enumerate(par_lines) if line.strip() == "71")
rows = par_lines[start:stop]
assert len(rows) == 15

raw = SOURCE.read_bytes()
newline = b"\r\n" if b"\r\n" in raw else b"\n"
lines = raw.decode("ascii").splitlines()
header = next(i for i, line in enumerate(lines) if line.strip() == '"          70", &')
assert lines[header + 16].strip() == '"          71", &'
old_rows = [line.strip().removeprefix('"').removesuffix('", &')
            for line in lines[header + 1:header + 16]]
assert all(old_rows[i] == rows[i] for i in FROZEN), "frozen Yb row changed"
replacement = [f'      "{row}", &' for row in rows]
assert all(re.fullmatch(r'      "[ 0-9.+-]+", &', row) for row in replacement)
lines[header + 1:header + 16] = replacement
SOURCE.write_bytes((newline.decode("ascii").join(lines) + newline.decode("ascii")).encode("ascii"))
print(f"embedded Yb rows from {EXPECTED_SHA}")
