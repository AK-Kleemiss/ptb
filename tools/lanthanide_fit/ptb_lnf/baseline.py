from __future__ import annotations

from pathlib import Path
import re

from .elements import LANTHANIDES


def _atompara_blocks(path: Path) -> set[int]:
    blocks = set()
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.isdigit():
            value = int(stripped)
            if 1 <= value <= 86:
                blocks.add(value)
    return blocks


def _element_blocks(path: Path) -> set[int]:
    blocks = set()
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.isdigit():
            value = int(stripped)
            if 1 <= value <= 118:
                blocks.add(value)
    return blocks


def verify(repo_root: str | Path) -> dict[str, object]:
    root = Path(repo_root)
    expected = {z for _, z in LANTHANIDES}
    atompara = _atompara_blocks(root / ".atompara")
    basis = _element_blocks(root / ".basis_vDZP")
    ecp = _element_blocks(root / ".ecp")
    default_files = (root / "source" / "default_files.f90").read_text()
    main = (root / "source" / "main.f90").read_text()
    reference = (root / "source" / "reference.inc").read_text()

    result = {
        "atompara_lanthanoids": sorted(atompara & expected),
        "basis_lanthanoids": sorted(basis & expected),
        "ecp_lanthanoids": sorted(ecp & expected),
        "main_reads_86_blocks": bool(re.search(r"do\s+i\s*=\s*1\s*,\s*86", main, re.IGNORECASE)),
        "default_atompara_declares_1378": "default_atompara(1378)" in default_files,
        "reference_has_la_lu": all(f"REF {symbol}" in reference for symbol, _ in LANTHANIDES if symbol != "Ce") and "REF CeH3" in reference,
    }
    result["ok"] = (
        set(result["atompara_lanthanoids"]) == expected
        and set(result["basis_lanthanoids"]) == expected
        and set(result["ecp_lanthanoids"]) == expected
        and result["main_reads_86_blocks"]
        and result["default_atompara_declares_1378"]
        and result["reference_has_la_lu"]
    )
    return result
