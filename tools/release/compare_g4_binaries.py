"""Compare Yb free-atom density matrices from the three release builds."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/lanthanide_fit"))
from ptb_lnf.denmat import read_denmat  # noqa: E402

BUILD = ROOT / "build/ion_g4_release_20261001"
WINDOWS = ROOT / "tools/lanthanide_fit/work/yb_config_fix_20260923/ion_g4_sizepass/embedded_smoke/Yb"
SOURCES = {
    "windows_x64": lambda q: WINDOWS / f"embedded_q{q}/ptb.denmat",
    "linux_x64": lambda q: BUILD / f"smoke_linux/q{q}/ptb.denmat",
    "macos_arm64": lambda q: BUILD / f"smoke_macos/arm64/q{q}/ptb.denmat",
    "macos_x64": lambda q: BUILD / f"smoke_macos/x86_64/q{q}/ptb.denmat",
}

records = []
for charge in (0, 2, 3):
    matrices = {name: read_denmat(path(charge)) for name, path in SOURCES.items()}
    reference = matrices["windows_x64"].density
    for name, matrix in matrices.items():
        count = float(np.trace(matrix.density @ matrix.overlap))
        maxdiff = float(np.max(np.abs(matrix.density - reference)))
        record = {"charge": charge, "platform": name, "electron_count": count,
                  "maxabs_density_vs_windows": maxdiff}
        assert abs(count - (24 - charge)) < 1e-6, record
        assert maxdiff < 1e-5, record
        records.append(record)

out = BUILD / "cross_platform_yb.json"
out.write_text(json.dumps({"status": "passed", "results": records}, indent=2) + "\n")
print(out.read_text())
