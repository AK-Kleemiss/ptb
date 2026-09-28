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


def _unit_vec(v: tuple[float, float, float]) -> tuple[float, float, float]:
    n = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    if n < 1e-12:
        return (0.0, 0.0, 1.0)
    return (v[0] / n, v[1] / n, v[2] / n)


def _cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _methyl_cap(
    base: tuple[float, float, float], outward: tuple[float, float, float]
) -> list[tuple[str, float, float, float]]:
    """3 H in a staggered sp3 tripod capping a terminal carbon at `base`,
    where `outward` is the bond direction from the carbon's parent atom."""
    z = _unit_vec(outward)
    ref = (0.0, 0.0, 1.0) if abs(z[2]) < 0.8 else (1.0, 0.0, 0.0)
    x = _unit_vec(_cross(ref, z))
    y = _unit_vec(_cross(z, x))
    bond = 1.09
    out = []
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        rx = math.cos(phi) * x[0] + math.sin(phi) * y[0]
        ry = math.cos(phi) * x[1] + math.sin(phi) * y[1]
        rz = math.cos(phi) * x[2] + math.sin(phi) * y[2]
        dx = -0.333 * z[0] + 0.943 * rx
        dy = -0.333 * z[1] + 0.943 * ry
        dz = -0.333 * z[2] + 0.943 * rz
        out.append(("H", base[0] + bond * dx, base[1] + bond * dy, base[2] + bond * dz))
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
            atoms.extend([("O", x, y, z), ("C", 1.2 * x, 1.2 * y, 1.2 * z), ("H", 1.32 * x + 0.25, 1.32 * y, 1.32 * z)])
    elif ligand_class == "amine":
        for x, y, z in _octahedral_points(2.45):
            atoms.extend([
                ("N", x, y, z),
                ("H", 1.08 * x, 1.08 * y + 0.20, 1.08 * z),
                ("H", 1.08 * x, 1.08 * y, 1.08 * z + 0.20),
                ("H", 1.08 * x + 0.20, 1.08 * y, 1.08 * z),
            ])
    elif ligand_class == "ether":
        for x, y, z in _octahedral_points(2.40):
            c1 = (1.2 * x + 0.35, 1.2 * y, 1.2 * z)
            c2 = (1.2 * x - 0.35, 1.2 * y, 1.2 * z)
            atoms.append(("O", x, y, z))
            atoms.append(("C", *c1))
            atoms.extend(_methyl_cap(c1, (0.35, 0.0, 0.0)))
            atoms.append(("C", *c2))
            atoms.extend(_methyl_cap(c2, (-0.35, 0.0, 0.0)))
    elif ligand_class == "phosphine_oxide":
        for x, y, z in _octahedral_points(2.20)[:4]:
            p = (1.35 * x, 1.35 * y, 1.35 * z)
            atoms.append(("O", x, y, z))
            atoms.append(("P", *p))
            for ox, oy, oz in ((0.9, 0.9, 0.0), (-0.9, 0.45, 0.78), (-0.9, 0.45, -0.78)):
                c = (p[0] + ox, p[1] + oy, p[2] + oz)
                atoms.append(("C", *c))
                atoms.extend(_methyl_cap(c, (ox, oy, oz)))
    elif ligand_class == "beta_diketonate":
        for x, y, z in _octahedral_points(2.30)[:4]:
            c1 = (1.2 * x, 1.2 * y, 1.2 * z)
            c2 = (1.4 * x, 1.4 * y + 0.45, 1.4 * z)
            outward = (c2[0] - c1[0], c2[1] - c1[1], c2[2] - c1[2])
            atoms.append(("O", x, y, z))
            atoms.append(("C", *c1))
            atoms.append(("C", *c2))
            atoms.extend(_methyl_cap(c2, outward))
    elif ligand_class == "cp_like":
        ring_radius, ch_bond = 1.7, 1.08
        for zsign in (-1.0, 1.0):
            for k in range(5):
                theta = angle + 2.0 * math.pi * k / 5.0
                cx, cy = ring_radius * math.cos(theta), ring_radius * math.sin(theta)
                atoms.append(("C", cx, cy, zsign * 2.1))
                hx, hy = (ring_radius + ch_bond) * math.cos(theta), (ring_radius + ch_bond) * math.sin(theta)
                atoms.append(("H", hx, hy, zsign * 2.1))
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
