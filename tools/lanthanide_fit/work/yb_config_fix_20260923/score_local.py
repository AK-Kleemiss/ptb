"""Molecular train/validation RMSE of Yb candidate .atompara files, scored locally
against the downloaded reference cubes: python score_local.py [--split S] PAR..."""
import os
import sys
import time

# 46 pool workers each running multithreaded BLAS oversubscribe the 48 threads.
for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(var, "1")
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from ptb_lnf.optimize import score_manifest  # noqa: E402
from scan_xi import BAS, PTB  # noqa: E402

MANIFEST = Path("D:/lnf_yb_data/yb_manifest_local.jsonl")

if __name__ == "__main__":
    args = sys.argv[1:]
    split = "train"
    if args[0] == "--split":
        split, args = args[1], args[2:]
    for par in args:
        t = time.time()
        s = score_manifest(MANIFEST, PTB, Path(par).resolve(), BAS, {}, workers=46, element="Yb", split=split)
        print(f"{split} {s:.6f} {time.time() - t:5.0f}s {par}", flush=True)
