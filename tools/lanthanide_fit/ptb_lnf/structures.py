from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import json
import math

from .elements import LANTHANIDES, COMMON_OXIDATION_STATES, DEFAULT_SPIN_MULTIPLICITY


@dataclass
class StructureRecord:
    case_id: str
    element: str
    atomic_number: int
    oxidation_state: int
    charge: int
    spin_multiplicity: int
    ligand_class: str
    conformer: int
    conformer_source: str
    xyz_path: str
    split: str


def _octahedral_points(radius: float) -> list[tuple[float, float, float]]:
    return [(radius, 0, 0), (-radius, 0, 0), (0, radius, 0), (0, -radius, 0), (0, 0, radius), (0, 0, -radius)]


def _rotate_xy(points: list[tuple[str, float, float, float]], angle: float) -> list[tuple[str, float, float, float]]:
    ca, sa = math.cos(angle), math.sin(angle)
    out = []
    for sym, x, y, z in points:
        out.append((sym, ca * x - sa * y, sa * x + ca * y, z))
    return out


def seed_atoms(element: str, ligand_class: str, conformer: int) -> list[tuple[str, float, float, float]]:
    atoms = [(element, 0.0, 0.0, 0.0)]
    angle = conformer * 0.35
    r = 2.35
    if ligand_class == "aqua":
        for x, y, z in _octahedral_points(r):
            atoms.extend([("O", x, y, z), ("H", x * 1.08, y * 1.08 + 0.22, z), ("H", x * 1.08, y, z * 1.08 + 0.22)])
    elif ligand_class == "halide":
        for x, y, z in _octahedral_points(2.55):
            atoms.append(("Cl", x, y, z))
    elif ligand_class == "oxide_hydroxide":
        for x, y, z in _octahedral_points(2.15):
            atoms.extend([("O", x, y, z), ("H", x * 1.15, y * 1.15, z * 1.15)])
    elif ligand_class == "carboxylate":
        for x, y, z in _octahedral_points(2.25)[:4]:
            atoms.extend([("O", x, y, z), ("C", 1.2 * x, 1.2 * y, 1.2 * z), ("O", 1.35 * x + 0.25, 1.35 * y, 1.35 * z)])
    elif ligand_class == "amine":
        for x, y, z in _octahedral_points(2.45):
            atoms.extend([("N", x, y, z), ("H", 1.08 * x, 1.08 * y + 0.20, 1.08 * z), ("H", 1.08 * x, 1.08 * y, 1.08 * z + 0.20)])
    elif ligand_class == "ether":
        for x, y, z in _octahedral_points(2.40):
            atoms.extend([("O", x, y, z), ("C", 1.2 * x + 0.35, 1.2 * y, 1.2 * z), ("C", 1.2 * x - 0.35, 1.2 * y, 1.2 * z)])
    elif ligand_class == "phosphine_oxide":
        for x, y, z in _octahedral_points(2.20)[:4]:
            atoms.extend([("O", x, y, z), ("P", 1.35 * x, 1.35 * y, 1.35 * z)])
    elif ligand_class == "beta_diketonate":
        for x, y, z in _octahedral_points(2.30)[:4]:
            atoms.extend([("O", x, y, z), ("C", 1.2 * x, 1.2 * y, 1.2 * z), ("C", 1.4 * x, 1.4 * y + 0.45, 1.4 * z)])
    elif ligand_class == "cp_like":
        for zsign in (-1.0, 1.0):
            for k in range(5):
                theta = angle + 2.0 * math.pi * k / 5.0
                atoms.append(("C", 1.7 * math.cos(theta), 1.7 * math.sin(theta), zsign * 2.1))
    return _rotate_xy(atoms, angle)


def write_xyz(path: Path, atoms: list[tuple[str, float, float, float]], comment: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="\n") as handle:
        handle.write(f"{len(atoms)}\n{comment}\n")
        for sym, x, y, z in atoms:
            handle.write(f"{sym:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def generate(out: str | Path, manifest: str | Path, conformers: int = 2) -> list[StructureRecord]:
    out = Path(out)
    ligand_classes = [
        "aqua", "halide", "oxide_hydroxide", "carboxylate", "amine",
        "ether", "phosphine_oxide", "beta_diketonate", "cp_like",
    ]
    records: list[StructureRecord] = []
    for symbol, z in LANTHANIDES:
        for ox in COMMON_OXIDATION_STATES[symbol]:
            for ligand in ligand_classes:
                for conf in range(conformers):
                    case_id = f"{symbol}_{ox:+d}_{ligand}_c{conf}"
                    xyz = out / symbol / case_id / f"{case_id}.xyz"
                    write_xyz(xyz, seed_atoms(symbol, ligand, conf), case_id)
                    split = "validation" if (z + ox + conf + len(ligand)) % 5 == 0 else "train"
                    records.append(StructureRecord(
                        case_id, symbol, z, ox, ox,
                        DEFAULT_SPIN_MULTIPLICITY.get((symbol, ox), 1),
                        ligand, conf, "template_seed_for_mace_osaka", str(xyz), split,
                    ))
    with Path(manifest).open("w", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    return records
