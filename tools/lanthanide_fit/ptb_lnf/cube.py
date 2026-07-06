from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import numpy as np


@dataclass(frozen=True)
class CubeGrid:
    comments: tuple[str, str]
    atoms: list[tuple[int, float, float, float, float]]
    origin: np.ndarray
    shape: tuple[int, int, int]
    axes: np.ndarray
    data: np.ndarray

    @property
    def voxel_volume(self) -> float:
        return float(abs(np.linalg.det(self.axes)))

    def coordinates(self) -> np.ndarray:
        nx, ny, nz = self.shape
        ix, iy, iz = np.meshgrid(
            np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij"
        )
        return (
            self.origin
            + ix[..., None] * self.axes[0]
            + iy[..., None] * self.axes[1]
            + iz[..., None] * self.axes[2]
        )


def read_cube(path: str | Path) -> CubeGrid:
    lines = Path(path).read_text().splitlines()
    comments = (lines[0], lines[1])
    nat_origin = lines[2].split()
    natoms = abs(int(nat_origin[0]))
    origin = np.array([float(x) for x in nat_origin[1:4]], dtype=float)
    shape = []
    axes = []
    for i in range(3):
        parts = lines[3 + i].split()
        shape.append(abs(int(parts[0])))
        axes.append([float(x) for x in parts[1:4]])
    atoms = []
    for i in range(natoms):
        parts = lines[6 + i].split()
        atoms.append((int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])))
    values = []
    for line in lines[6 + natoms :]:
        values.extend(float(x) for x in line.split())
    data = np.array(values, dtype=float).reshape(tuple(shape))
    return CubeGrid(comments, atoms, origin, tuple(shape), np.array(axes), data)


def write_cube(path: str | Path, grid: CubeGrid, data: np.ndarray, comments: tuple[str, str] | None = None) -> None:
    comments = comments or grid.comments
    with Path(path).open("w", newline="\n") as handle:
        handle.write(comments[0].rstrip() + "\n")
        handle.write(comments[1].rstrip() + "\n")
        handle.write(f"{len(grid.atoms):5d} {grid.origin[0]:12.6f} {grid.origin[1]:12.6f} {grid.origin[2]:12.6f}\n")
        for n, axis in zip(grid.shape, grid.axes):
            handle.write(f"{n:5d} {axis[0]:12.6f} {axis[1]:12.6f} {axis[2]:12.6f}\n")
        for z, charge, x, y, zcoord in grid.atoms:
            handle.write(f"{z:5d} {charge:12.6f} {x:12.6f} {y:12.6f} {zcoord:12.6f}\n")
        flat = np.asarray(data, dtype=float).ravel()
        for i in range(0, flat.size, 6):
            handle.write("".join(f" {value:13.5e}" for value in flat[i : i + 6]) + "\n")


def compare(reference: CubeGrid, candidate: CubeGrid) -> dict[str, float]:
    if reference.shape != candidate.shape or not np.allclose(reference.axes, candidate.axes) or not np.allclose(reference.origin, candidate.origin):
        raise ValueError("cube grids do not match")
    diff = candidate.data - reference.data
    vol = reference.voxel_volume
    return {
        "rmse": float(math.sqrt(np.mean(diff * diff))),
        "mae": float(np.mean(np.abs(diff))),
        "linf": float(np.max(np.abs(diff))),
        "reference_integral": float(np.sum(reference.data) * vol),
        "candidate_integral": float(np.sum(candidate.data) * vol),
        "integral_delta": float(np.sum(diff) * vol),
    }
