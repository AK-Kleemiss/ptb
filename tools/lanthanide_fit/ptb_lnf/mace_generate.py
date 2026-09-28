from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import random
import shutil
from typing import Iterable

import numpy as np

from .elements import LANTHANIDES, COMMON_OXIDATION_STATES, DEFAULT_SPIN_MULTIPLICITY
from .structures import write_xyz


BOHR_PER_ANGSTROM = 1.8897259886
PILOT_PER_ELEMENT = 4
TRAINING_PER_ELEMENT = 1000


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
    # Recipe names map to a fixed ligand multiset (only the order is
    # shuffled per structure), so it can be recovered from ligand_class
    # alone even without a saved manifest, which keeps is_stale_record able
    # to validate structures scanned back from a bare XYZ tree.
    recipe_ligands = RECIPE_LIGANDS_BY_NAME.get(ligand_class)
    ligand_recipe = "+".join(recipe_ligands) if recipe_ligands else "recovered_from_xyz"
    return MaceStructureRecord(
        case_id=case_id,
        element=element,
        atomic_number=z,
        oxidation_state=ox,
        charge=complex_charge(ox, ligand_recipe) if recipe_ligands else ox,
        spin_multiplicity=DEFAULT_SPIN_MULTIPLICITY.get((element, ox), 1),
        ligand_class=ligand_class,
        ligand_recipe=ligand_recipe,
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


# Per-ligand-token atom formula (excludes the central lanthanide), used to
# detect structures written before a ligand-geometry fix added atoms that
# weren't there yet (e.g. missing hydrogens). Keep in sync with build_complex
# and add_cp_ring/add_monodentate/add_bidentate_carboxylate.
LIGAND_FORMULA: dict[str, dict[str, int]] = {
    "water": {"O": 1, "H": 2},
    "hydroxide": {"O": 1, "H": 1},
    "chloride": {"Cl": 1},
    "amine": {"N": 1, "H": 3},
    "ether": {"O": 1, "C": 2, "H": 6},
    "phosphine_oxide": {"O": 1, "P": 1, "C": 3, "H": 9},
    "carboxylate": {"O": 2, "C": 1, "H": 1},
    "beta_diketonate": {"O": 2, "C": 3, "H": 7},
    "cp": {"C": 5, "H": 5},
}

# Formal charge each ligand token contributes to the overall complex.
# water/amine/ether/phosphine_oxide are neutral donors; the rest are
# monoanionic (Cl-, OH-, RCOO-, acac-, Cp-).
LIGAND_CHARGE: dict[str, int] = {
    "water": 0,
    "hydroxide": -1,
    "chloride": -1,
    "amine": 0,
    "ether": 0,
    "phosphine_oxide": 0,
    "carboxylate": -1,
    "beta_diketonate": -1,
    "cp": -1,
}


def complex_charge(oxidation_state: int, ligand_recipe: str) -> int:
    """Overall molecular charge: metal oxidation state plus the sum of
    each ligand's formal charge. Using the bare oxidation state alone
    (ignoring anionic ligands) previously made ORCA's charge/multiplicity
    parity check fail for any recipe with an odd number of anionic
    ligands (chloride/hydroxide/carboxylate/beta_diketonate/cp)."""
    charge = oxidation_state
    for ligand in ligand_recipe.split("+"):
        charge += LIGAND_CHARGE.get(ligand, 0)
    return charge


def expected_formula(ligand_recipe: str) -> dict[str, int] | None:
    formula: dict[str, int] = {}
    for ligand in ligand_recipe.split("+"):
        per_ligand = LIGAND_FORMULA.get(ligand)
        if per_ligand is None:
            return None
        for sym, n in per_ligand.items():
            formula[sym] = formula.get(sym, 0) + n
    return formula


def actual_formula(atoms: list[tuple[str, float, float, float]], element: str) -> dict[str, int]:
    formula: dict[str, int] = {}
    for sym, *_ in atoms:
        if sym == element:
            continue
        formula[sym] = formula.get(sym, 0) + 1
    return formula


def is_stale_record(record: "MaceStructureRecord") -> bool:
    """True if the on-disk XYZ predates a ligand-geometry fix and no longer
    matches what its own recorded ligand_recipe should produce (e.g. a
    pre-fix structure missing hydrogens on a terminal ligand carbon), or if
    the relaxed geometry has a dissociated atom (MACE occasionally flings
    one off during optimization; predates the max-nearest-neighbor filter
    added to generate())."""
    try:
        atoms = read_xyz(Path(record.xyz_path))
    except (OSError, ValueError):
        return True
    if max_nearest_neighbor_distance(atoms) > 2.8:
        return True
    expected = expected_formula(record.ligand_recipe)
    if expected is None:
        return False
    return actual_formula(atoms, record.element) != expected


def _dedupe_by_case_id(records: list[MaceStructureRecord]) -> list[MaceStructureRecord]:
    """Collapse manifest rows that share a case_id (e.g. from merging
    manifests across machines) down to one row each, keeping whichever row's
    uniqueness_signature matches the structure actually on disk today (the
    other geometry was silently overwritten when both used the same path)."""
    by_id: dict[str, list[MaceStructureRecord]] = {}
    for record in records:
        by_id.setdefault(record.case_id, []).append(record)
    kept = []
    dropped = 0
    for group in by_id.values():
        if len(group) == 1:
            kept.append(group[0])
            continue
        dropped += len(group) - 1
        try:
            current_sig = uniqueness_signature(read_xyz(Path(group[0].xyz_path)))
        except (OSError, ValueError):
            current_sig = None
        match = next((r for r in group if r.uniqueness_signature == current_sig), None)
        kept.append(match if match is not None else group[-1])
    if dropped:
        print(
            f"Resume: dropped {dropped} duplicate manifest row(s) sharing a case_id "
            "with another record; kept the row matching the structure currently on disk."
        )
    return kept


def load_existing_records(out: Path, manifest: Path, model_path: Path) -> list[MaceStructureRecord]:
    if manifest.exists():
        records = []
        for line in manifest.read_text().splitlines():
            if line.strip():
                records.append(MaceStructureRecord(**json.loads(line)))
    elif out.exists():
        records = []
        for xyz in sorted(out.glob("*/*/*.xyz")):
            try:
                records.append(record_from_xyz(xyz, model_path))
            except (KeyError, ValueError) as exc:
                print(f"Skipping unrecognized XYZ during resume scan: {xyz} ({exc})")
    else:
        records = []

    records = _dedupe_by_case_id(records)

    kept = []
    discarded = 0
    for record in records:
        if is_stale_record(record):
            discarded += 1
            case_dir = Path(record.xyz_path).parent
            if case_dir.exists():
                shutil.rmtree(case_dir)
            continue
        kept.append(record)
    if discarded:
        print(
            f"Resume: discarded {discarded} pre-fix structure(s) whose atoms no "
            "longer match their recorded ligand_recipe; they will be regenerated."
        )

    repaired = 0
    for record in kept:
        correct_charge = complex_charge(record.oxidation_state, record.ligand_recipe)
        correct_mult = DEFAULT_SPIN_MULTIPLICITY.get((record.element, record.oxidation_state), record.spin_multiplicity)
        if record.charge != correct_charge or record.spin_multiplicity != correct_mult:
            record.charge = correct_charge
            record.spin_multiplicity = correct_mult
            repaired += 1
    if repaired:
        print(
            f"Resume: repaired charge/spin_multiplicity metadata on {repaired} "
            "record(s) (geometry untouched; only the ORCA charge/mult were stale)."
        )
    return kept


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


def add_methyl_cap(
    atoms: list[tuple[str, np.ndarray]],
    base: np.ndarray,
    outward: np.ndarray,
    rng: random.Random,
) -> None:
    """Cap a terminal carbon at `base` with 3 H in a staggered sp3 tripod.

    `outward` is the bond direction from the carbon's parent atom to `base`;
    the H atoms are placed at the tetrahedral angle from the reversed bond.
    """
    x, y, z = orthonormal_frame(outward)
    bond = 1.09
    phase = rng.random() * 2.0 * math.pi
    for k in range(3):
        phi = phase + 2.0 * math.pi * k / 3.0
        radial = math.cos(phi) * x + math.sin(phi) * y
        direction = -0.333 * z + 0.943 * radial
        atoms.append(("H", base + bond * direction))


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
        for offset in (1.70 * x, -0.85 * x + 1.47 * y, -0.85 * x - 1.47 * y):
            c = p + offset
            atoms.append(("C", c))
            add_methyl_cap(atoms, c, offset, rng)
    elif donor == "O" and tail == "ether":
        for offset in (1.35 * x + 0.35 * z, -1.35 * x + 0.35 * z):
            c = donor_pos + offset
            atoms.append(("C", c))
            add_methyl_cap(atoms, c, offset, rng)


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
        for offset in (1.30 * y + 0.20 * z, -1.30 * y + 0.20 * z):
            c2 = c + offset
            atoms.append(("C", c2))
            add_methyl_cap(atoms, c2, offset, rng)
        # Central-carbon vinylic H (as in real acac-'s backbone CH). Without
        # it, O2C3H6 has an even proton count, so as the intended -1 anion
        # it comes out an odd-electron (open-shell) fragment -- exactly the
        # kind of charge/multiplicity parity mismatch ORCA rejects.
        atoms.append(("H", c + 1.09 * z))
    else:
        atoms.append(("H", c + 1.09 * z))


def add_cp_ring(atoms: list[tuple[str, np.ndarray]], direction: np.ndarray, radius: float, rng: random.Random) -> None:
    x, y, z = orthonormal_frame(direction)
    center = z * (radius + 0.35)
    phase = rng.random() * 2.0 * math.pi
    ring_radius = 1.42
    ch_bond = 1.09
    for k in range(5):
        theta = phase + 2.0 * math.pi * k / 5.0
        radial = math.cos(theta) * x + math.sin(theta) * y
        atoms.append(("C", center + ring_radius * radial))
        atoms.append(("H", center + (ring_radius + ch_bond) * radial))


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

RECIPE_LIGANDS_BY_NAME: dict[str, tuple[str, ...]] = {name: tuple(ligands) for name, ligands in RECIPES}


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


def max_nearest_neighbor_distance(atoms: list[tuple[str, float, float, float]]) -> float:
    """Worst-case distance from any atom to its closest neighbor. A large
    value means some atom drifted off during relaxation and is effectively
    dissociated from the rest of the molecule -- a failure mode MACE's
    optimizer occasionally produces that min_distance (too-close check)
    can't catch."""
    coords = np.array([[x, y, z] for _, x, y, z in atoms], dtype=float)
    if len(coords) < 2:
        return 0.0
    worst = 0.0
    for i in range(len(coords)):
        nearest = min(
            float(np.linalg.norm(coords[i] - coords[j])) for j in range(len(coords)) if j != i
        )
        worst = max(worst, nearest)
    return worst


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

    # conf numbers must be tracked explicitly rather than derived from
    # len(records for element): once purge/dedup remove arbitrary entries,
    # existing conf numbers are no longer a contiguous 0..len-1 range, and a
    # length-based index can collide with an already-used higher conf.
    used_confs: dict[str, set[int]] = {}
    for record in records:
        used_confs.setdefault(record.element, set()).add(record.conformer)

    for symbol, z in LANTHANIDES:
        ox_states = COMMON_OXIDATION_STATES[symbol]
        attempts = 0
        target = per_element
        symbol_confs = used_confs.setdefault(symbol, set())
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
                if max_nearest_neighbor_distance(atoms) > 2.8:
                    continue
            sig = uniqueness_signature(atoms)
            if sig in seen:
                continue
            seen.add(sig)
            conf = 0
            while conf in symbol_confs:
                conf += 1
            symbol_confs.add(conf)
            case_id = f"{symbol}_{ox:+d}_{recipe_name}_mace{conf:03d}"
            xyz_path = out / symbol / case_id / f"{case_id}.xyz"
            write_xyz(xyz_path, atoms, f"{case_id}; model={model_path.name}; relaxed={relaxed}")
            split = "validation" if (z + ox + conf) % 5 == 0 else "train"
            ligand_recipe = "+".join(ligands)
            records.append(MaceStructureRecord(
                case_id=case_id,
                element=symbol,
                atomic_number=z,
                oxidation_state=ox,
                charge=complex_charge(ox, ligand_recipe),
                spin_multiplicity=DEFAULT_SPIN_MULTIPLICITY.get((symbol, ox), 1),
                ligand_class=recipe_name,
                ligand_recipe=ligand_recipe,
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
