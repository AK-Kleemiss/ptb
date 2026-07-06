from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import random
from typing import Iterable

import numpy as np

from .elements import LANTHANIDES, COMMON_OXIDATION_STATES, DEFAULT_SPIN_MULTIPLICITY
from .structures import write_xyz


BOHR_PER_ANGSTROM = 1.8897259886
PILOT_PER_ELEMENT = 4
TRAINING_PER_ELEMENT = 150


@dataclass
class MaceStructureRecord:
    case_id: str
    element: str
    atomic_number: int
    oxidation_state: int
    charge: int
    spin_multiplicity: int
    ligand_class: str
    ligand_recipe: str
    conformer: int
    conformer_source: str
    xyz_path: str
    split: str
    mace_model: str
    mace_relaxed: bool
    mace_energy_ev: float | None
    uniqueness_signature: str
    min_distance_angstrom: float


def read_xyz(path: Path) -> list[tuple[str, float, float, float]]:
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"empty XYZ file: {path}")
    nat = int(lines[0].strip())
    atoms: list[tuple[str, float, float, float]] = []
    for line in lines[2:2 + nat]:
        sym, x, y, z = line.split()[:4]
        atoms.append((sym, float(x), float(y), float(z)))
    if len(atoms) != nat:
        raise ValueError(f"XYZ atom count mismatch in {path}: expected {nat}, got {len(atoms)}")
    return atoms


def parse_case_id(case_id: str) -> tuple[str, int, str, int]:
    parts = case_id.split("_")
    if len(parts) < 4 or not parts[-1].startswith("mace"):
        raise ValueError(f"cannot parse MACE case id: {case_id}")
    element = parts[0]
    oxidation_state = int(parts[1])
    conformer = int(parts[-1][4:])
    ligand_class = "_".join(parts[2:-1])
    return element, oxidation_state, ligand_class, conformer


def record_from_xyz(path: Path, model_path: Path) -> MaceStructureRecord:
    case_id = path.stem
    element, ox, ligand_class, conformer = parse_case_id(case_id)
    z = dict(LANTHANIDES)[element]
    atoms = read_xyz(path)
    lines = path.read_text().splitlines()
    comment = lines[1] if len(lines) > 1 else ""
    relaxed = "relaxed=True" in comment or "relaxed=true" in comment
    split = "validation" if (z + ox + conformer) % 5 == 0 else "train"
    return MaceStructureRecord(
        case_id=case_id,
        element=element,
        atomic_number=z,
        oxidation_state=ox,
        charge=ox,
        spin_multiplicity=DEFAULT_SPIN_MULTIPLICITY.get((element, ox), 1),
        ligand_class=ligand_class,
        ligand_recipe="recovered_from_xyz",
        conformer=conformer,
        conformer_source="mace_osaka26_diverse_relaxed" if relaxed else "diverse_seed_unrelaxed",
        xyz_path=str(path),
        split=split,
        mace_model=str(model_path),
        mace_relaxed=relaxed,
        mace_energy_ev=None,
        uniqueness_signature=uniqueness_signature(atoms),
        min_distance_angstrom=min_distance(atoms),
    )


def load_existing_records(out: Path, manifest: Path, model_path: Path) -> list[MaceStructureRecord]:
    if manifest.exists():
        records = []
        for line in manifest.read_text().splitlines():
            if line.strip():
                records.append(MaceStructureRecord(**json.loads(line)))
        return records
    if not out.exists():
        return []
    records = []
    for xyz in sorted(out.glob("*/*/*.xyz")):
        try:
            records.append(record_from_xyz(xyz, model_path))
        except (KeyError, ValueError) as exc:
            print(f"Skipping unrecognized XYZ during resume scan: {xyz} ({exc})")
    return records


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("zero vector")
    return v / n


def fibonacci_directions(n: int, rng: random.Random) -> list[np.ndarray]:
    offset = 2.0 / n
    inc = math.pi * (3.0 - math.sqrt(5.0))
    dirs = []
    phase = rng.random() * 2.0 * math.pi
    for k in range(n):
        y = ((k * offset) - 1.0) + offset / 2.0
        r = math.sqrt(max(0.0, 1.0 - y * y))
        phi = ((k + 1) % n) * inc + phase
        dirs.append(np.array([math.cos(phi) * r, y, math.sin(phi) * r]))
    rng.shuffle(dirs)
    return dirs


def orthonormal_frame(direction: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z = _unit(direction)
    ref = np.array([0.0, 0.0, 1.0]) if abs(z[2]) < 0.8 else np.array([1.0, 0.0, 0.0])
    x = _unit(np.cross(ref, z))
    y = _unit(np.cross(z, x))
    return x, y, z


def add_monodentate(
    atoms: list[tuple[str, np.ndarray]],
    direction: np.ndarray,
    donor: str,
    radius: float,
    rng: random.Random,
    tail: str | None = None,
) -> None:
    x, y, z = orthonormal_frame(direction)
    jitter = rng.uniform(-0.08, 0.08)
    donor_pos = z * (radius + jitter)
    atoms.append((donor, donor_pos))
    if donor == "O" and tail == "water":
        atoms.append(("H", donor_pos + 0.78 * z + 0.52 * x))
        atoms.append(("H", donor_pos + 0.78 * z - 0.52 * x))
    elif donor == "O" and tail == "hydroxide":
        atoms.append(("H", donor_pos + 0.98 * z))
    elif donor == "N":
        atoms.append(("H", donor_pos + 0.72 * z + 0.55 * x))
        atoms.append(("H", donor_pos + 0.72 * z - 0.55 * x))
        atoms.append(("H", donor_pos + 0.72 * z + 0.55 * y))
    elif donor == "O" and tail == "phosphine_oxide":
        p = donor_pos + 1.52 * z
        atoms.append(("P", p))
        atoms.append(("C", p + 1.70 * x))
        atoms.append(("C", p - 0.85 * x + 1.47 * y))
        atoms.append(("C", p - 0.85 * x - 1.47 * y))
    elif donor == "O" and tail == "ether":
        atoms.append(("C", donor_pos + 1.35 * x + 0.35 * z))
        atoms.append(("C", donor_pos - 1.35 * x + 0.35 * z))


def add_bidentate_carboxylate(
    atoms: list[tuple[str, np.ndarray]],
    direction: np.ndarray,
    radius: float,
    rng: random.Random,
    beta_diketonate: bool = False,
) -> None:
    x, y, z = orthonormal_frame(direction)
    spread = rng.uniform(0.72, 0.95)
    center = z * (radius + rng.uniform(-0.06, 0.06))
    o1 = center + spread * x
    o2 = center - spread * x
    atoms.append(("O", o1))
    atoms.append(("O", o2))
    c = center + 1.12 * z
    atoms.append(("C", c))
    if beta_diketonate:
        atoms.append(("C", c + 1.30 * y + 0.20 * z))
        atoms.append(("C", c - 1.30 * y + 0.20 * z))
    else:
        atoms.append(("O", c + 1.18 * z))


def add_cp_ring(atoms: list[tuple[str, np.ndarray]], direction: np.ndarray, radius: float, rng: random.Random) -> None:
    x, y, z = orthonormal_frame(direction)
    center = z * (radius + 0.35)
    phase = rng.random() * 2.0 * math.pi
    for k in range(5):
        theta = phase + 2.0 * math.pi * k / 5.0
        atoms.append(("C", center + 1.42 * math.cos(theta) * x + 1.42 * math.sin(theta) * y))


RECIPES = [
    ("aqua8", ["water"] * 8),
    ("mixed_aqua_chloride", ["water"] * 5 + ["chloride"] * 2),
    ("hydroxide_aqua", ["hydroxide"] * 3 + ["water"] * 5),
    ("amine_aqua", ["amine"] * 4 + ["water"] * 4),
    ("ether_nitrate_like", ["ether"] * 3 + ["carboxylate"] * 2),
    ("carboxylate_aqua", ["carboxylate"] * 3 + ["water"] * 2),
    ("beta_diketonate_aqua", ["beta_diketonate"] * 3 + ["water"] * 1),
    ("phosphine_oxide_halide", ["phosphine_oxide"] * 3 + ["chloride"] * 3),
    ("cp_halide", ["cp"] * 2 + ["chloride"] * 3),
]


def build_complex(element: str, recipe_name: str, ligands: list[str], rng: random.Random) -> list[tuple[str, float, float, float]]:
    atoms: list[tuple[str, np.ndarray]] = [(element, np.zeros(3))]
    dirs = fibonacci_directions(max(len(ligands), 4), rng)
    base_radius = rng.uniform(2.25, 2.65)
    for direction, ligand in zip(dirs, ligands):
        if ligand == "water":
            add_monodentate(atoms, direction, "O", base_radius, rng, "water")
        elif ligand == "hydroxide":
            add_monodentate(atoms, direction, "O", base_radius - 0.10, rng, "hydroxide")
        elif ligand == "chloride":
            add_monodentate(atoms, direction, "Cl", base_radius + 0.25, rng)
        elif ligand == "amine":
            add_monodentate(atoms, direction, "N", base_radius + 0.08, rng)
        elif ligand == "ether":
            add_monodentate(atoms, direction, "O", base_radius + 0.05, rng, "ether")
        elif ligand == "phosphine_oxide":
            add_monodentate(atoms, direction, "O", base_radius, rng, "phosphine_oxide")
        elif ligand == "carboxylate":
            add_bidentate_carboxylate(atoms, direction, base_radius, rng, False)
        elif ligand == "beta_diketonate":
            add_bidentate_carboxylate(atoms, direction, base_radius + 0.05, rng, True)
        elif ligand == "cp":
            add_cp_ring(atoms, direction, base_radius + 0.20, rng)
    # Break accidental high symmetry without generating pure rotations.
    distorted = []
    for sym, xyz in atoms:
        if sym == element:
            distorted.append((sym, xyz))
        else:
            distorted.append((sym, xyz + np.array([rng.gauss(0, 0.04), rng.gauss(0, 0.04), rng.gauss(0, 0.04)])))
    return [(sym, float(pos[0]), float(pos[1]), float(pos[2])) for sym, pos in distorted]


def min_distance(atoms: list[tuple[str, float, float, float]]) -> float:
    coords = np.array([[x, y, z] for _, x, y, z in atoms], dtype=float)
    if len(coords) < 2:
        return 999.0
    dmin = 999.0
    for i in range(len(coords)):
        for j in range(i):
            dmin = min(dmin, float(np.linalg.norm(coords[i] - coords[j])))
    return dmin


def uniqueness_signature(atoms: list[tuple[str, float, float, float]], decimals: int = 2) -> str:
    symbols = [sym for sym, *_ in atoms]
    coords = np.array([[x, y, z] for _, x, y, z in atoms], dtype=float)
    parts = []
    for sym in sorted(set(symbols)):
        idx = [i for i, s in enumerate(symbols) if s == sym]
        parts.append(f"{sym}{len(idx)}")
    distances = []
    for i in range(len(coords)):
        for j in range(i):
            pair = "-".join(sorted((symbols[i], symbols[j])))
            distances.append((pair, round(float(np.linalg.norm(coords[i] - coords[j])), decimals)))
    distances.sort()
    return "|".join(parts) + "::" + ",".join(f"{p}:{d:.{decimals}f}" for p, d in distances)


def to_ase_atoms(atoms: list[tuple[str, float, float, float]]):
    from ase import Atoms

    return Atoms([sym for sym, *_ in atoms], positions=[[x, y, z] for _, x, y, z in atoms])


def from_ase_atoms(atoms) -> list[tuple[str, float, float, float]]:
    return [
        (sym, float(pos[0]), float(pos[1]), float(pos[2]))
        for sym, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions())
    ]


def make_calculator(model_path: Path, device: str, dtype: str):
    import torch
    from mace.calculators import MACECalculator

    if device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError(
            f"device={device!r} was requested, but this Python environment has "
            f"torch {torch.__version__} with cuda_available=False. Rebuild the "
            "MACE venv with CUDA PyTorch, or rerun with --device cpu."
        )
    try:
        return MACECalculator(model_paths=str(model_path), device=device, default_dtype=dtype)
    except TypeError:
        return MACECalculator(model_path=str(model_path), device=device, default_dtype=dtype)


def relax_with_mace(
    atoms: list[tuple[str, float, float, float]],
    model_path: Path,
    device: str,
    max_steps: int,
    fmax: float,
    dtype: str = "float64",
) -> tuple[list[tuple[str, float, float, float]], float | None, bool]:
    from ase.optimize import LBFGS

    ase_atoms = to_ase_atoms(atoms)
    ase_atoms.calc = make_calculator(model_path, device, dtype)
    opt = LBFGS(ase_atoms, logfile=None)
    opt.run(fmax=fmax, steps=max_steps)
    energy = float(ase_atoms.get_potential_energy())
    return from_ase_atoms(ase_atoms), energy, True


def generate(
    out: Path,
    manifest: Path,
    model_path: Path,
    per_element: int,
    seed: int,
    device: str,
    max_steps: int,
    fmax: float,
    dtype: str = "float64",
    no_relax: bool = False,
    resume: bool = False,
) -> list[MaceStructureRecord]:
    rng = random.Random(seed)
    records = load_existing_records(out, manifest, model_path) if resume else []
    seen: set[str] = {record.uniqueness_signature for record in records}
    out.mkdir(parents=True, exist_ok=True)

    for symbol, z in LANTHANIDES:
        ox_states = COMMON_OXIDATION_STATES[symbol]
        attempts = 0
        target = per_element
        while len([r for r in records if r.element == symbol]) < target:
            attempts += 1
            if attempts > target * 80:
                raise RuntimeError(f"could not generate enough unique structures for {symbol}")
            ox = rng.choice(ox_states)
            recipe_name, ligands = rng.choice(RECIPES)
            ligands = list(ligands)
            rng.shuffle(ligands)
            atoms = build_complex(symbol, recipe_name, ligands, rng)
            if min_distance(atoms) < 0.62:
                continue
            sig0 = uniqueness_signature(atoms)
            if sig0 in seen:
                continue
            relaxed = False
            energy = None
            if not no_relax:
                atoms, energy, relaxed = relax_with_mace(atoms, model_path, device, max_steps, fmax, dtype=dtype)
                if min_distance(atoms) < 0.55:
                    continue
            sig = uniqueness_signature(atoms)
            if sig in seen:
                continue
            seen.add(sig)
            conf = len([r for r in records if r.element == symbol])
            case_id = f"{symbol}_{ox:+d}_{recipe_name}_mace{conf:03d}"
            xyz_path = out / symbol / case_id / f"{case_id}.xyz"
            write_xyz(xyz_path, atoms, f"{case_id}; model={model_path.name}; relaxed={relaxed}")
            split = "validation" if (z + ox + conf) % 5 == 0 else "train"
            records.append(MaceStructureRecord(
                case_id=case_id,
                element=symbol,
                atomic_number=z,
                oxidation_state=ox,
                charge=ox,
                spin_multiplicity=DEFAULT_SPIN_MULTIPLICITY.get((symbol, ox), 1),
                ligand_class=recipe_name,
                ligand_recipe="+".join(ligands),
                conformer=conf,
                conformer_source="mace_osaka26_diverse_relaxed" if relaxed else "diverse_seed_unrelaxed",
                xyz_path=str(xyz_path),
                split=split,
                mace_model=str(model_path),
                mace_relaxed=relaxed,
                mace_energy_ev=energy,
                uniqueness_signature=sig,
                min_distance_angstrom=min_distance(atoms),
            ))

    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(asdict(record), sort_keys=True) + "\n")
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate diverse La-Lu training structures with MACE-Osaka relaxation.")
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--per-element", type=int, default=None)
    parser.add_argument(
        "--training-set",
        action="store_true",
        help=f"Generate the production target of {TRAINING_PER_ELEMENT} structures per lanthanide.",
    )
    parser.add_argument("--seed", type=int, default=20260706)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--dtype", choices=("float64", "float32"), default="float64")
    parser.add_argument("--max-steps", type=int, default=120)
    parser.add_argument("--fmax", type=float, default=0.08)
    parser.add_argument("--no-relax", action="store_true", help="Generate diverse seeds without MACE relaxation; useful for fast tests only.")
    parser.add_argument("--resume", action="store_true", help="Reuse existing manifest/XYZ records and generate only missing per-element structures.")
    args = parser.parse_args()
    per_element = args.per_element if args.per_element is not None else (
        TRAINING_PER_ELEMENT if args.training_set else PILOT_PER_ELEMENT
    )

    records = generate(
        out=args.out,
        manifest=args.manifest,
        model_path=args.model,
        per_element=per_element,
        seed=args.seed,
        device=args.device,
        dtype=args.dtype,
        max_steps=args.max_steps,
        fmax=args.fmax,
        no_relax=args.no_relax,
        resume=args.resume,
    )
    by_split: dict[str, int] = {}
    by_element: dict[str, int] = {}
    for record in records:
        by_split[record.split] = by_split.get(record.split, 0) + 1
        by_element[record.element] = by_element.get(record.element, 0) + 1
    print(json.dumps({
        "records": len(records),
        "per_element": per_element,
        "by_split": by_split,
        "by_element": by_element,
        "manifest": str(args.manifest),
        "out": str(args.out),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
