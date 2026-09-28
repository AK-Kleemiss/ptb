from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np

from .cube import CubeGrid, read_cube, write_cube


@dataclass
class Term:
    powers: tuple[int, int, int]
    transform: float
    primitives: list[tuple[float, float]]


@dataclass
class AO:
    index: int
    atom_index: int
    shell_index: int
    angular_l: int
    m_index: int
    norm: float
    terms: list[Term]


@dataclass
class Denmat:
    atoms: list[tuple[int, np.ndarray]]
    aos: list[AO]
    density: np.ndarray
    overlap: np.ndarray


def _unpack(values: list[float], n: int) -> np.ndarray:
    out = np.zeros((n, n), dtype=float)
    k = 0
    for i in range(n):
        for j in range(i + 1):
            out[i, j] = values[k]
            out[j, i] = values[k]
            k += 1
    return out


def read_denmat(path: str | Path) -> Denmat:
    lines = Path(path).read_text().splitlines()
    if not lines or lines[0].strip() != "# PTB_DENMAT 1":
        raise ValueError("not a PTB_DENMAT v1 file")
    i = 1
    natoms = int(lines[i].split()[1]); i += 1
    atoms = []
    for _ in range(natoms):
        parts = lines[i].split(); i += 1
        atoms.append((int(parts[2]), np.array([float(parts[3]), float(parts[4]), float(parts[5])], dtype=float)))
    ndim = int(lines[i].split()[1]); i += 1
    aos: list[AO] = []
    while lines[i].startswith("AO "):
        parts = lines[i].split(); i += 1
        index = int(parts[1])
        atom_index = int(parts[2])
        shell_index = int(parts[3])
        angular_l = int(parts[4])
        m_index = int(parts[5])
        nterms = int(parts[6])
        norm = float(parts[7])
        terms: list[Term] = []
        for _ in range(nterms):
            parts = lines[i].split(); i += 1
            if parts[0] != "TERM":
                raise ValueError(f"expected TERM at line {i}")
            powers = (int(parts[1]), int(parts[2]), int(parts[3]))
            nprim = int(parts[4])
            transform = float(parts[5])
            primitives = []
            for _ in range(nprim):
                exp_s, coeff_s = lines[i].split(); i += 1
                primitives.append((float(exp_s), float(coeff_s)))
            if abs(transform) > 1e-14:
                terms.append(Term(powers, transform, primitives))
        aos.append(AO(index, atom_index, shell_index, angular_l, m_index, norm, terms))
    if len(aos) != ndim:
        raise ValueError(f"AO count {len(aos)} does not match NDIM {ndim}")
    if not lines[i].startswith("DENSITY_PACKED"):
        raise ValueError("missing DENSITY_PACKED section")
    packed_n = int(lines[i].split()[1]); i += 1
    density = [float(lines[i + j]) for j in range(packed_n)]
    i += packed_n
    if not lines[i].startswith("OVERLAP_PACKED"):
        raise ValueError("missing OVERLAP_PACKED section")
    packed_n = int(lines[i].split()[1]); i += 1
    overlap = [float(lines[i + j]) for j in range(packed_n)]
    return Denmat(atoms, aos, _unpack(density, ndim), _unpack(overlap, ndim))


def evaluate_density(denmat: Denmat, coords: np.ndarray, chunk_size: int = 25000) -> np.ndarray:
    flat_coords = coords.reshape((-1, 3))
    out = np.empty(flat_coords.shape[0], dtype=float)
    centers = [xyz for _, xyz in denmat.atoms]
    for start in range(0, flat_coords.shape[0], chunk_size):
        stop = min(start + chunk_size, flat_coords.shape[0])
        pts = flat_coords[start:stop]
        phi = np.empty((stop - start, len(denmat.aos)), dtype=float)
        # All AOs of an atom share rel/r2, and all Cartesian terms of a shell
        # share one contracted radial part and a handful of monomials: build
        # each once per chunk, not once per AO and term (this loop was >90 %
        # of a fitting trial; PTB itself is ~0.4 s of ~6 s per structure).
        rel_cache: dict = {}
        radial_cache: dict = {}
        poly_cache: dict = {}
        for col, ao in enumerate(denmat.aos):
            a = ao.atom_index
            if a not in rel_cache:
                rel = pts - centers[a - 1]
                rel_cache[a] = (rel, np.einsum("ij,ij->i", rel, rel))
            rel, r2 = rel_cache[a]
            val = np.zeros(stop - start, dtype=float)
            for term in ao.terms:
                rkey = (a, tuple(term.primitives))
                contracted = radial_cache.get(rkey)
                if contracted is None:
                    contracted = np.zeros(stop - start, dtype=float)
                    for exponent, coeff in term.primitives:
                        contracted += coeff * np.exp(-exponent * r2)
                    radial_cache[rkey] = contracted
                pkey = (a, term.powers)
                poly = poly_cache.get(pkey)
                if poly is None:
                    px, py, pz = term.powers
                    poly = (rel[:, 0] ** px) * (rel[:, 1] ** py) * (rel[:, 2] ** pz)
                    poly_cache[pkey] = poly
                val += term.transform * poly * contracted
            phi[:, col] = ao.norm * val
        # einsum("pi,ij,pj->p") with three operands runs as an unblocked loop;
        # the matrix product goes through BLAS.
        out[start:stop] = np.einsum("pi,pi->p", phi @ denmat.density, phi)
    return out.reshape(coords.shape[:-1])


def denmat_to_cube(denmat_path: str | Path, grid_cube: str | Path, out_cube: str | Path) -> None:
    denmat = read_denmat(denmat_path)
    grid = read_cube(grid_cube)
    density = density_on_grid(denmat, grid)
    write_cube(
        out_cube,
        grid,
        density,
        comments=("PTB density projected from -denmat", f"grid copied from {Path(grid_cube).name}"),
    )


def density_on_grid(denmat: Denmat, grid: CubeGrid) -> np.ndarray:
    """Evaluate a PTB density matrix at every point of an existing cube grid.

    Optimisation only needs the numerical density values for an RMSE, not a
    persistent candidate cube.  Keeping this operation in memory avoids a
    large format/write/read cycle for every fitted trial while retaining the
    reference cube's exact origin, axes, and resolution.
    """
    return evaluate_density(denmat, grid.coordinates())
