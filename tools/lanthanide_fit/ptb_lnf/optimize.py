from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy.optimize import differential_evolution, minimize, minimize_scalar

from .cube import compare, read_cube
from .denmat import density_on_grid, read_denmat
from .iongate import gate_violation
from .orca import read_manifest

# Every flag source/main.f90's argument loop recognises. A file path containing
# any of these as a substring is unusable with binaries built before the loop
# switched from index() to an exact match -- see the check in score_manifest.
PTB_CLI_FLAGS = ("-stda", "-purify", "-par", "-bas", "-chrg", "-uhf", "-denmat")

# Each element occupies exactly 15 data rows in .atompara, read in this exact
# order by source/main.f90:131-145. Grouped by physical role instead of
# treated as one undifferentiated block, since a single uniform scale factor
# across e.g. orbital exponents and energy levels alike has no physical
# justification -- these are independent axes in the original PTB fit for
# every non-lanthanide element (see parascript.f90, which perturbs one named
# parameter/group at a time).
#
# "energy" (ener_par1/2/4/5/6) is intentionally EXCLUDED from GROUP_NAMES
# (left at ROW_GROUPS below only so the row-counting logic stays correct --
# _copy_candidate_atompara's `group_scales.get(group, 1.0)` leaves any group
# absent from GROUP_NAMES untouched at scale=1.0). Verified 2026-07-26: these
# rows are ~1e-6 (negligible placeholder) for every d-block AND f-block metal
# checked across the periodic table (Ti, Fe, Zn, Ag, La, Gd, Lu, Au, Pb all
# show max|value| ~1e-6, vs real substantial values -- up to ~30 -- for
# main-group elements H/C/O/Si), and scaling this group across a huge range
# (0.1 to 60x) produced BIT-IDENTICAL RMSE against real ORCA reference
# densities for two different lanthanides (Eu, Pr) -- i.e. it has zero
# measurable effect on the density PTB outputs for metals as a class, not
# just a lanthanide-specific gap. Fitting it was pure wasted search
# dimensionality: multiple elements' DE runs pegged at whatever arbitrary
# bound was set (e.g. Eu/Gd/Ho/Lu/Tb/Tm/Yb all identically hit 6.4999 in one
# run) since there was no gradient signal to constrain the search at all.
GROUP_NAMES = ["exponent_scaling", "shell_exponents", "shell_config", "shell_response"]
ROW_GROUPS = [
    "energy",              # ener_par1
    "energy",              # ener_par2
    "exponent_scaling",    # expscal(1)
    "energy",              # ener_par6
    "energy",              # ener_par4
    "energy",              # ener_par5
    "exponent_scaling",    # expscal(2)
    "shell_exponents",     # shell_xi
    "shell_config",        # shell_cnf1
    "shell_config",        # shell_cnf2
    "shell_config",        # shell_cnf3
    "exponent_scaling",    # expscal(3)
    "shell_config",        # shell_cnf4
    "shell_response",      # shell_resp(...,1)
    "shell_response",      # shell_resp(...,2)
]


def _copy_candidate_atompara(template: Path, out: Path, group_scales: dict[str, float], target_z: int | None = None) -> None:
    # If target_z is given, only that element's block is scaled (used when
    # fitting per-element, since a structure's PTB behavior only depends on
    # the elements actually present in it -- other lanthanide blocks are
    # irrelevant to that structure's score either way, but leaving them
    # untouched keeps the candidate file meaningful to inspect standalone).
    lines = template.read_text().splitlines()
    active = False
    row_idx = -1
    result = []
    for line in lines:
        stripped = line.strip()
        if stripped.isdigit():
            z = int(stripped)
            active = (z == target_z) if target_z is not None else (57 <= z <= 71)
            row_idx = -1
            result.append(line)
            continue
        if active and stripped and not stripped.startswith("#"):
            row_idx += 1
            group = ROW_GROUPS[row_idx] if 0 <= row_idx < len(ROW_GROUPS) else None
            scale = group_scales.get(group, 1.0) if group else 1.0
            vals = stripped.split()
            try:
                nums = [float(v) for v in vals]
            except ValueError:
                result.append(line)
                continue
            result.append(" ".join(f"{v * scale:15.10f}" if abs(v) > 1e-12 else f"{v:15.10f}" for v in nums))
        else:
            if active:
                row_idx += 1
            result.append(line)
    with out.open("w", newline="\n") as handle:
        handle.write("\n".join(result) + "\n")


def free_parameter_positions(template: str | Path, target_z: int) -> list[tuple[int, int]]:
    """Returns [(row_idx, col_idx), ...] for every nonzero, non-"energy" value
    in target_z's .atompara block -- the free-fit parameter space (one
    independent multiplicative scale per real, nonzero value, instead of one
    shared scale per physical-role group). "energy" rows and structural-zero
    placeholders are excluded for the same reason they're excluded from
    GROUP_NAMES: no measurable effect on output density (energy) or not a
    real parameter at all (zero is a fixed shell-absence marker, not a fitted
    value -- see _copy_candidate_atompara's `abs(v) > 1e-12` guard)."""
    lines = Path(template).read_text().splitlines()
    active = False
    row_idx = -1
    positions: list[tuple[int, int]] = []
    for line in lines:
        stripped = line.strip()
        if stripped.isdigit():
            active = int(stripped) == target_z
            row_idx = -1
            continue
        if active and stripped and not stripped.startswith("#"):
            row_idx += 1
            group = ROW_GROUPS[row_idx] if 0 <= row_idx < len(ROW_GROUPS) else None
            if group == "energy" or group is None:
                continue
            try:
                vals = [float(v) for v in stripped.split()]
            except ValueError:
                continue
            positions.extend((row_idx, col) for col, v in enumerate(vals) if abs(v) > 1e-12)
        elif active:
            row_idx += 1
    return positions


def _copy_candidate_atompara_free(
    template: Path, out: Path, target_z: int, positions: list[tuple[int, int]], value_scales
) -> None:
    """Like _copy_candidate_atompara, but scales each (row_idx, col_idx)
    position in positions by its own independent factor from value_scales
    (same order), instead of one factor per physical-role group. Every value
    not in positions (energy rows, structural zeros) is copied unchanged."""
    scale_map = dict(zip(positions, (float(s) for s in value_scales)))
    lines = template.read_text().splitlines()
    active = False
    row_idx = -1
    result = []
    for line in lines:
        stripped = line.strip()
        if stripped.isdigit():
            active = int(stripped) == target_z
            row_idx = -1
            result.append(line)
            continue
        if active and stripped and not stripped.startswith("#"):
            row_idx += 1
            vals = stripped.split()
            try:
                nums = [float(v) for v in vals]
            except ValueError:
                result.append(line)
                continue
            scaled = [v * scale_map[(row_idx, col)] if (row_idx, col) in scale_map else v for col, v in enumerate(nums)]
            result.append(" ".join(f"{v:15.10f}" for v in scaled))
        else:
            if active:
                row_idx += 1
            result.append(line)
    with out.open("w", newline="\n") as handle:
        handle.write("\n".join(result) + "\n")


def build_merged_atompara(template: str | Path, out: str | Path, element_group_scales: dict[str, dict[str, float]]) -> None:
    """Writes one final .atompara where each lanthanide element's block has
    its own independently-fitted per-group scale factors applied (elements
    not present in element_group_scales are left unscaled)."""
    from .elements import SYMBOL_TO_Z

    z_scales = {SYMBOL_TO_Z[symbol]: scales for symbol, scales in element_group_scales.items()}
    lines = Path(template).read_text().splitlines()
    active_scales: dict[str, float] | None = None
    row_idx = -1
    result = []
    for line in lines:
        stripped = line.strip()
        if stripped.isdigit():
            z = int(stripped)
            active_scales = z_scales.get(z)
            row_idx = -1
            result.append(line)
            continue
        if active_scales is not None and stripped and not stripped.startswith("#"):
            row_idx += 1
            group = ROW_GROUPS[row_idx] if 0 <= row_idx < len(ROW_GROUPS) else None
            scale = active_scales.get(group, 1.0) if group else 1.0
            vals = stripped.split()
            try:
                nums = [float(v) for v in vals]
            except ValueError:
                result.append(line)
                continue
            result.append(" ".join(f"{v * scale:15.10f}" if abs(v) > 1e-12 else f"{v:15.10f}" for v in nums))
        else:
            if active_scales is not None:
                row_idx += 1
            result.append(line)
    with Path(out).open("w", newline="\n") as handle:
        handle.write("\n".join(result) + "\n")


def _score_one_row(args: tuple[dict, str, str, str]) -> float | None:
    """Runs PTB for a single manifest row against a candidate .atompara,
    returning the RMSE against its reference cube (None on failure). A
    top-level function so ProcessPoolExecutor can pickle it to workers."""
    row, candidate_str, basis_str, ptb_exe_str = args
    with tempfile.TemporaryDirectory() as work_s:
        work = Path(work_s)
        xyz = work / Path(row["xyz_path"]).name
        shutil.copy2(row["xyz_path"], xyz)
        denmat = work / "ptb.denmat"
        cmd = [ptb_exe_str, str(xyz), "-par", candidate_str, "-bas", basis_str, "-denmat", str(denmat)]
        if row.get("charge") is not None:
            cmd.extend(["-chrg", str(row["charge"])])
        # PTB defaults to a closed-shell (nopen=0) reference unless told
        # otherwise (source/main.f90:55); without this every lanthanide
        # here would run as the wrong electronic state, making density
        # comparison against the ORCA reference meaningless.
        spin_multiplicity = row.get("spin_multiplicity")
        if spin_multiplicity is not None:
            cmd.extend(["-uhf", str(spin_multiplicity - 1)])
        proc = subprocess.run(cmd, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if proc.returncode != 0 or not denmat.exists():
            return None
        # The reference grid defines the comparison exactly.  Do not write a
        # candidate cube only to parse it again: the optimiser needs its RMSE
        # alone, and the in-memory density is both faster and more precise.
        reference = read_cube(row["reference_cube"])
        candidate_density = density_on_grid(read_denmat(denmat), reference)
        diff = candidate_density - reference.data
        return float(np.sqrt(np.mean(diff * diff)))


def score_manifest(
    manifest: str | Path,
    ptb_exe: str | Path,
    atompara_template: str | Path,
    basis: str | Path,
    group_scales: dict[str, float],
    workers: int | None = None,
    element: str | None = None,
    split: str = "train",
) -> float:
    rows = [row for row in read_manifest(manifest) if row.get("reference_cube") and row.get("split") == split]
    if element is not None:
        rows = [row for row in rows if row.get("element") == element]
    if not rows:
        raise ValueError(f"manifest has no {split} rows with reference_cube" + (f" for element {element}" if element else ""))
    # _score_one_row runs PTB with cwd set to a per-structure temp dir, so a
    # relative --basis/--ptb resolves against a directory that does not contain
    # it. PTB does not fail on a missing -bas file: it prints one line and
    # silently falls back to the basis compiled into the binary, exit code 0.
    # That turned the 2026-09-21 Yb run into a 12-hour fit of the upstream
    # basis it was supposed to replace. Resolve here, once, for every caller.
    ptb_exe = Path(ptb_exe).resolve()
    basis = Path(basis).resolve()
    if not basis.is_file():
        raise FileNotFoundError(f"basis file not found: {basis} -- PTB would silently use its built-in basis")
    if not ptb_exe.is_file():
        raise FileNotFoundError(f"ptb executable not found: {ptb_exe}")
    # main.f90 matched its CLI flags with index(), i.e. as a substring anywhere
    # in the argument, so a path merely CONTAINING one silently redirected that
    # flag to the next argument -- /tmp/yb-basis-d2.../basis_vDZP contains
    # '-bas' and cost this campaign three 12-hour jobs before it was found.
    # Fixed in main.f90, but every already-built binary still has it, so refuse
    # such a path here rather than trust the binary in use.
    for label, path in (("basis", basis), ("ptb", ptb_exe)):
        hit = next((f for f in PTB_CLI_FLAGS if f in str(path)), None)
        if hit is not None:
            raise ValueError(
                f"{label} path contains the PTB flag {hit!r}: {path} -- older ptb binaries "
                f"match flags as substrings and would silently ignore it. Rename the directory."
            )
    workers = workers or os.cpu_count() or 1
    from .elements import SYMBOL_TO_Z

    target_z = SYMBOL_TO_Z[element] if element is not None else None
    with tempfile.TemporaryDirectory() as tmp_s:
        tmp = Path(tmp_s)
        candidate = tmp / "candidate.atompara"
        _copy_candidate_atompara(Path(atompara_template), candidate, group_scales, target_z=target_z)
        # PTB itself is fast (single-digit CPU-seconds/structure) but there
        # are ~1500 train structures per trial value and differential_evolution
        # needs 100+ trials -- sequential execution is infeasible (4+ hours
        # per trial). Parallelize across structures instead, using all
        # available cores for one trial at a time (DE itself stays
        # workers=1 to avoid oversubscribing).
        args = [(row, str(candidate), str(basis), str(ptb_exe)) for row in rows]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(_score_one_row, args))
        # A single non-finite density RMSE must not turn the complete
        # objective into NaN: differential evolution cannot rank or improve
        # NaN objectives.  Treat it exactly like a failed PTB calculation,
        # retaining the finite structures' score and applying the existing
        # failure penalty below.
        metrics = [r for r in results if r is not None and np.isfinite(r)]
        failures = len(results) - len(metrics)
        if not metrics:
            return 1e6
        return float(np.mean(metrics) + failures * 10.0)


def optimize_scale(
    manifest: str | Path,
    ptb_exe: str | Path,
    atompara: str | Path,
    basis: str | Path,
    workers: int | None = None,
    element: str | None = None,
    bounds: tuple[float, float] = (0.80, 1.20),
    popsize: int = 15,
    maxiter: int = 100,
) -> dict:
    trial = 0
    t_start = time.monotonic()
    tag = f"[{element}] " if element else ""

    def objective(x: np.ndarray) -> float:
        nonlocal trial
        trial += 1
        t0 = time.monotonic()
        group_scales = dict(zip(GROUP_NAMES, (float(v) for v in x)))
        score = score_manifest(manifest, ptb_exe, atompara, basis, group_scales, workers=workers, element=element)
        elapsed = time.monotonic() - t0
        total = time.monotonic() - t_start
        scales_str = " ".join(f"{g}={s:.4f}" for g, s in group_scales.items())
        print(
            f"{tag}trial {trial:4d}  {scales_str}  score={score:.6f}  "
            f"trial_time={elapsed:6.1f}s  total_time={total/60:7.1f}min",
            file=sys.stderr,
            flush=True,
        )
        return score

    de_bounds = [bounds] * len(GROUP_NAMES)
    # popsize/maxiter default to SciPy's own defaults, but going from 1 to 5
    # free parameters means the population (popsize * ndim) is 5x larger per
    # generation than the earlier single-scalar fits -- pass smaller values
    # explicitly here when wall-clock budget matters more than exhaustively
    # covering the 5-D space.
    de = differential_evolution(objective, bounds=de_bounds, polish=False, workers=1, popsize=popsize, maxiter=maxiter)
    # Keep the local polish inside the same bounds as the global search.
    # Nelder-Mead was unconstrained and could walk far outside the requested
    # range, which made "bounds" misleading for the grouped fit.
    local = minimize(
        objective,
        x0=de.x,
        method="Powell",
        bounds=de_bounds,
        options={"maxiter": 80, "maxfev": 160},
    )
    # Bounded Powell is not guaranteed to be monotonically improving when its
    # line searches get clipped against a bound (observed in practice: it can
    # walk from de.x to a worse boundary-pegged point). Never let the "polish"
    # step regress below the DE global best -- fall back to de.x/de.fun if so.
    if local.fun <= de.fun:
        final_x, final_fun, local_success, local_message = local.x, float(local.fun), bool(local.success), str(local.message)
    else:
        final_x, final_fun = de.x, float(de.fun)
        local_success, local_message = False, f"rejected: local polish regressed ({local.fun:.6f} > {de.fun:.6f}); kept DE global best"
        print(f"{tag}WARNING: local polish regressed ({local.fun:.6f} > {de.fun:.6f}) -- keeping DE global best instead", file=sys.stderr, flush=True)
    lo, hi = bounds
    margin = 0.02 * (hi - lo)
    final_scales = dict(zip(GROUP_NAMES, (float(v) for v in final_x)))
    for group, value in final_scales.items():
        if value <= lo + margin or value >= hi - margin:
            print(
                f"{tag}WARNING: final {group}={value:.6f} is within 2% of bounds {bounds} "
                "-- the true optimum may lie outside; consider widening bounds and rerunning.",
                file=sys.stderr,
                flush=True,
            )
    return {
        "element": element,
        "scales": final_scales,
        "score": final_fun,
        "global_score": float(de.fun),
        "local_success": local_success,
        "local_message": local_message,
        "bounds": list(bounds),
    }


def _write_free_checkpoint(
    path: str | Path, element: str, positions: list[tuple[int, int]],
    x: np.ndarray, score: float, sweep: int, param_idx: int, trial: int,
    sweep_start_score: float | None = None,
) -> None:
    data = {
        "element": element,
        "positions": [[r, c] for r, c in positions],
        "value_scales": [float(v) for v in x],
        "score": float(score),
        "sweep_start_score": None if sweep_start_score is None else float(sweep_start_score),
        "sweep": sweep,
        "param_idx": param_idx,
        "trials": trial,
    }
    tmp = Path(str(path) + ".tmp")
    with tmp.open("w", newline="\n") as handle:
        json.dump(data, handle, indent=2)
    tmp.replace(path)


def optimize_free(
    manifest: str | Path,
    ptb_exe: str | Path,
    atompara: str | Path,
    basis: str | Path,
    element: str,
    workers: int | None = None,
    bounds: tuple[float, float] = (0.5, 2.0),
    passes: int = 1,
    evals_per_param: int = 8,
    checkpoint: str | Path | None = None,
    resume_from: str | Path | None = None,
    exclude_positions: set[tuple[int, int]] | None = None,
    min_sweep_relative_improvement: float | None = None,
    patience: int = 1,
    verify_resume: bool = False,
    relative_step: float | None = None,
    ion_targets: dict | None = None,
) -> dict:
    """Coordinate-descent refinement of every individual nonzero, non-"energy"
    value in one element's .atompara block, starting from `atompara` (expected
    to already carry that element's fitted group scales from optimize_scale --
    this is a local polish on top of that result, not a from-scratch search).

    DE is infeasible here: free_parameter_positions finds one free dimension
    per real value (order 50-90 per element), and a population-based global
    search at that dimensionality would need orders of magnitude more trials
    than the already multi-day 4-group fit. Coordinate descent -- bounded 1-D
    search of one parameter at a time, holding all others fixed at their
    current best, repeated for `passes` sweeps -- costs
    O(passes * ndim * evals_per_param) trials, which is controllable, and each
    parameter's result is checkpointed durably before moving to the next, so a
    SLURM time-limit kill (observed to actually happen on this cluster's 48h
    budget for the grouped fit) loses at most one in-progress parameter, not
    the whole run. `resume_from` picks back up from a prior checkpoint's exact
    sweep/parameter position instead of restarting from scale=1.0 everywhere.
    """
    from .elements import SYMBOL_TO_Z

    target_z = SYMBOL_TO_Z[element]
    positions = free_parameter_positions(atompara, target_z)
    excluded = exclude_positions or set()
    positions = [position for position in positions if position not in excluded]
    if not positions:
        raise ValueError(f"{element}: no active free parameters remain after exclusions")
    ndim = len(positions)
    trial = 0
    t_start = time.monotonic()
    tag = f"[{element}-free] "

    x = np.ones(ndim)
    start_sweep = 0
    start_param = 0
    best_score = None
    # Score at the start of the sweep a resume lands in the middle of. Without
    # it the early-stop test measured only the post-resume part of that sweep
    # (job 526221: 0.35 % instead of the true 2.1 %, stopped after one sweep).
    resumed_sweep_start = None
    resumed_mid_sweep = False
    if resume_from is not None and Path(resume_from).exists():
        ckpt = json.loads(Path(resume_from).read_text())
        if [list(p) for p in positions] != ckpt["positions"]:
            raise ValueError(f"{tag}checkpoint positions don't match current atompara/element -- refusing to resume from a stale checkpoint")
        x = np.array(ckpt["value_scales"])
        best_score = float(ckpt["score"])
        start_sweep = ckpt["sweep"]
        start_param = ckpt["param_idx"] + 1
        resumed_mid_sweep = start_param < ndim
        resumed_sweep_start = ckpt.get("sweep_start_score")
        if start_param >= ndim:
            start_param = 0
            start_sweep += 1
        print(f"{tag}resuming from checkpoint at sweep={start_sweep} param={start_param} score={best_score:.6f}", file=sys.stderr, flush=True)

    def score_at(x_vec: np.ndarray) -> float:
        nonlocal trial
        trial += 1
        t0 = time.monotonic()
        with tempfile.TemporaryDirectory() as tmp_s:
            candidate = Path(tmp_s) / "candidate.atompara"
            _copy_candidate_atompara_free(Path(atompara), candidate, target_z, positions, x_vec)
            violation = 0.0
            if ion_targets is not None:
                # Checked first: a wrong free-ion configuration is rejected
                # without the molecular scoring (see iongate.py). 1.0 is far
                # above any density RMSE, and the violation keeps a gradient
                # back towards the feasible region. Up to hard_limit electrons
                # outside the windows are instead charged weight per electron
                # on top of the molecular score, so a search can trade the two.
                violation, _ = gate_violation(ptb_exe, candidate, basis, ion_targets)
            if ion_targets is not None and violation > ion_targets.get("hard_limit", 0.0):
                score = 1.0 + violation
            else:
                score = score_manifest(manifest, ptb_exe, candidate, basis, {}, workers=workers, element=element)
                if violation > 0.0:
                    score += ion_targets.get("weight", 0.02) * violation
        elapsed = time.monotonic() - t0
        total = time.monotonic() - t_start
        print(
            f"{tag}trial {trial:5d}  score={score:.6f}  trial_time={elapsed:6.1f}s  total_time={total/60:7.1f}min",
            file=sys.stderr, flush=True,
        )
        return score

    if best_score is None:
        best_score = score_at(x)
        print(f"{tag}starting score (free scales=1.0 on top of the group-scaled input) = {best_score:.6f}", file=sys.stderr, flush=True)
    elif verify_resume:
        resumed_score = score_at(x)
        if not np.isfinite(resumed_score) or resumed_score >= 1e5:
            raise RuntimeError(f"{tag}resumed checkpoint does not produce a finite training score; refusing to search from it")
        best_score = float(resumed_score)
        print(f"{tag}verified resumed score={best_score:.6f}", file=sys.stderr, flush=True)

    quiet_sweeps = 0
    for sweep in range(start_sweep, passes):
        score_at_sweep_start = best_score
        sweep_start_unknown = False
        if sweep == start_sweep and resumed_mid_sweep:
            if resumed_sweep_start is None:
                sweep_start_unknown = True  # pre-fix checkpoint: can't judge this sweep
            else:
                score_at_sweep_start = float(resumed_sweep_start)
        param_range = range(start_param if sweep == start_sweep else 0, ndim)
        for i in param_range:
            def objective_1d(v: float, i: int = i) -> float:
                trial_x = x.copy()
                trial_x[i] = v
                return score_at(trial_x)

            parameter_bounds = bounds
            if relative_step is not None:
                if relative_step <= 0:
                    raise ValueError("relative_step must be positive")
                parameter_bounds = (
                    max(bounds[0], x[i] * (1.0 - relative_step)),
                    min(bounds[1], x[i] * (1.0 + relative_step)),
                )
                if parameter_bounds[0] >= parameter_bounds[1]:
                    raise ValueError(f"{tag}empty local bounds for parameter {i}: {parameter_bounds}")
            res = minimize_scalar(objective_1d, bounds=parameter_bounds, method="bounded", options={"maxiter": evals_per_param})
            if res.fun <= best_score:
                x[i] = res.x
                best_score = float(res.fun)
            row, col = positions[i]
            print(
                f"{tag}sweep {sweep + 1}/{passes} param {i + 1}/{ndim} row={row} col={col} "
                f"-> scale={x[i]:.4f}  best_score={best_score:.6f}",
                file=sys.stderr, flush=True,
            )
            if checkpoint is not None:
                _write_free_checkpoint(
                    checkpoint, element, positions, x, best_score, sweep, i, trial,
                    None if sweep_start_unknown else score_at_sweep_start,
                )

        relative_improvement = (score_at_sweep_start - best_score) / max(abs(score_at_sweep_start), 1e-15)
        print(
            f"{tag}sweep {sweep + 1}/{passes} complete: relative_improvement={relative_improvement:.3e}",
            file=sys.stderr, flush=True,
        )
        if sweep_start_unknown:
            print(f"{tag}sweep {sweep + 1} was resumed mid-way from a checkpoint without sweep_start_score; not counted for early stop", file=sys.stderr, flush=True)
        elif min_sweep_relative_improvement is not None:
            if relative_improvement < min_sweep_relative_improvement:
                quiet_sweeps += 1
            else:
                quiet_sweeps = 0
            if quiet_sweeps >= patience:
                print(
                    f"{tag}early stop after {quiet_sweeps} low-improvement sweep(s) "
                    f"(< {min_sweep_relative_improvement:.3e})",
                    file=sys.stderr, flush=True,
                )
                break

    return {
        "element": element,
        "positions": [[r, c] for r, c in positions],
        "value_scales": [float(v) for v in x],
        "score": float(best_score),
        "bounds": list(bounds),
        "passes": passes,
        "evals_per_param": evals_per_param,
        "ndim": ndim,
        "trials": trial,
        "excluded_positions": [[r, c] for r, c in sorted(excluded)],
        "min_sweep_relative_improvement": min_sweep_relative_improvement,
        "verify_resume": verify_resume,
        "relative_step": relative_step,
        "ion_targets": ion_targets,
    }


def score_free_manifest(
    manifest: str | Path, ptb_exe: str | Path, atompara: str | Path, basis: str | Path,
    element: str, result: dict, workers: int | None = None, split: str = "validation",
) -> float:
    """Scores an optimize_free result against a manifest split (e.g. the held-
    out validation set) -- the free-fit analogue of score_manifest, since its
    result schema (positions/value_scales) isn't the group_scales dict
    score_manifest itself expects."""
    from .elements import SYMBOL_TO_Z

    target_z = SYMBOL_TO_Z[element]
    positions = [tuple(p) for p in result["positions"]]
    value_scales = result["value_scales"]
    with tempfile.TemporaryDirectory() as tmp_s:
        candidate = Path(tmp_s) / "candidate.atompara"
        _copy_candidate_atompara_free(Path(atompara), candidate, target_z, positions, value_scales)
        return score_manifest(manifest, ptb_exe, candidate, basis, {}, workers=workers, element=element, split=split)


def apply_free_result(template: str | Path, out: str | Path, element: str, result: dict) -> None:
    """Writes a full .atompara with one element's free-fit result applied on
    top of `template` (expected to be the group-scaled merged file so the
    free-fit's positions/scales, which are relative to that starting point,
    land correctly). Other elements' blocks are copied through unchanged."""
    from .elements import SYMBOL_TO_Z

    target_z = SYMBOL_TO_Z[element]
    positions = [tuple(p) for p in result["positions"]]
    value_scales = result["value_scales"]
    _copy_candidate_atompara_free(Path(template), Path(out), target_z, positions, value_scales)
