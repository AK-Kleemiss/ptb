from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools" / "lanthanide_fit"))

from ptb_lnf.baseline import verify
from ptb_lnf.cube import compare, read_cube, write_cube


def test_cube_roundtrip_and_self_compare(tmp_path):
    cube_path = tmp_path / "a.cube"
    with cube_path.open("w", newline="\n") as handle:
        handle.write("\n".join([
            "comment 1",
            "comment 2",
            "    1 0.0 0.0 0.0",
            "    2 1.0 0.0 0.0",
            "    2 0.0 1.0 0.0",
            "    2 0.0 0.0 1.0",
            "    57 57.0 0.0 0.0 0.0",
            " 1.0 2.0 3.0 4.0 5.0 6.0",
            " 7.0 8.0",
            "",
        ]))
    cube = read_cube(cube_path)
    out = tmp_path / "b.cube"
    write_cube(out, cube, cube.data)
    metrics = compare(read_cube(cube_path), read_cube(out))
    assert metrics["rmse"] == 0.0
    assert np.isclose(metrics["reference_integral"], 36.0)


def test_lanthanide_baseline_verifier():
    result = verify(ROOT)
    assert result["main_reads_86_blocks"]
    assert set(result["atompara_lanthanoids"]) == set(range(57, 72))
    assert set(result["basis_lanthanoids"]) == set(range(57, 72))
    assert set(result["ecp_lanthanoids"]) == set(range(57, 72))
