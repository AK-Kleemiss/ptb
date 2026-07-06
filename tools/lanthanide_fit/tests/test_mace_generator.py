from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools" / "lanthanide_fit"))

from ptb_lnf.mace_generate import generate


def test_diverse_seed_generation_without_mace(tmp_path):
    manifest = tmp_path / "manifest.jsonl"
    records = generate(
        out=tmp_path / "xyz",
        manifest=manifest,
        model_path=tmp_path / "missing.model",
        per_element=2,
        seed=123,
        device="cpu",
        max_steps=1,
        fmax=0.5,
        no_relax=True,
    )
    assert len(records) == 30
    signatures = [r.uniqueness_signature for r in records]
    assert len(signatures) == len(set(signatures))
    ligand_classes = {r.ligand_class for r in records}
    assert len(ligand_classes) > 4
    rows = [json.loads(line) for line in manifest.read_text().splitlines()]
    assert len(rows) == len(records)
    assert all(Path(row["xyz_path"]).exists() for row in rows)
