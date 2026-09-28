"""Embed the verified ion_g2 Yb parameter block in PTB defaults."""
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "source/default_files.f90"
PAR = Path(__file__).resolve().parent / "ion_g2/final_candidate.atompara"
EXPECTED_SHA = "6eb12794d40838ce9b96050df769be7daa1c6e311d941b5dd6b048864ea51c39"

assert hashlib.sha256(PAR.read_bytes()).hexdigest() == EXPECTED_SHA
par_lines = PAR.read_text(encoding="ascii").splitlines()
start = next(i for i, line in enumerate(par_lines) if line.strip() == "70") + 1
stop = next(i for i, line in enumerate(par_lines) if line.strip() == "71")
rows = par_lines[start:stop]
assert len(rows) == 15
assert float(rows[6].split()[6]) == 1.1392905882

raw = SOURCE.read_bytes()
newline = b"\r\n" if b"\r\n" in raw else b"\n"
lines = raw.decode("ascii").splitlines()
header = next(i for i, line in enumerate(lines) if line.strip() == '"          70", &')
assert lines[header + 16].strip() == '"          71", &'
replacement = [f'      "{row}", &' for row in rows]
for row in replacement:
    assert re.fullmatch(r'      "[ 0-9.+-]+", &', row), row
lines[header + 1:header + 16] = replacement
SOURCE.write_bytes((newline.decode("ascii").join(lines) + newline.decode("ascii")).encode("ascii"))
print(f"embedded Yb rows from {EXPECTED_SHA}")
