"""python gate_check.py TARGETS PAR... -> gate violation and free-ion l-populations."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from ptb_lnf.iongate import gate_violation, load_targets  # noqa: E402
from scan_xi import BAS, PTB  # noqa: E402

targets = load_targets(sys.argv[1])
for par in sys.argv[2:]:
    v, rep = gate_violation(PTB, Path(par), BAS, targets)
    pops = " | ".join(f"q{q}: " + " ".join(f"{l} {p[l]:5.2f}" for l in "spdf") for q, p in rep.items())
    print(f"viol {v:6.3f}  {pops}  {Path(par).name}", flush=True)
