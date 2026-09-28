from __future__ import annotations

import re
from pathlib import Path

SUCCESS_MARKER = "ORCA TERMINATED NORMALLY"
TIME_LIMIT_MARKER = "DUE TO TIME LIMIT"


def _tail_contains(path: Path, marker: str, tail_bytes: int = 8192) -> bool:
    size = path.stat().st_size
    with path.open("rb") as handle:
        if size > tail_bytes:
            handle.seek(-tail_bytes, 2)
        data = handle.read()
    return marker.encode() in data


def classify_case(case_dir: Path) -> tuple[str, str]:
    """Returns (status, detail). status is one of: not_run, ok,
    ok_missing_cube, charge_multiplicity_error, input_error, orca_error,
    unknown_failure, unreadable."""
    case_id = case_dir.name
    out_path = case_dir / f"{case_id}.out"
    if not out_path.exists():
        return "not_run", ""
    try:
        terminated_ok = _tail_contains(out_path, SUCCESS_MARKER)
    except OSError as exc:
        return "unreadable", str(exc)

    if terminated_ok:
        # orca_plot's interactive menu can silently no-op (see the
        # module-load/menu-sequence history in orca.py); a successful SCF
        # run with no .cube file means the density was never extracted.
        if any(case_dir.glob("*.cube")):
            return "ok", ""
        return "ok_missing_cube", ""

    try:
        text = out_path.read_text(errors="replace")
    except OSError as exc:
        return "unreadable", str(exc)

    for line in text.splitlines():
        if "is even and number of electrons" in line or "is odd and number of electrons" in line:
            return "charge_multiplicity_error", line.strip()
    if "INPUT ERROR" in text:
        detail_lines = [l.strip() for l in text.splitlines() if l.strip()]
        idx = next((i for i, l in enumerate(detail_lines) if l == "INPUT ERROR"), None)
        detail = detail_lines[idx + 1] if idx is not None and idx + 1 < len(detail_lines) else "INPUT ERROR"
        return "input_error", detail
    if "aborting the run" in text.lower() or "TERMINATING THE RUN" in text:
        return "orca_error", "run aborted"
    return "unknown_failure", ""


def find_timed_out_task_indices(slurm_log_dir: str | Path, array_job_id: str) -> set[int]:
    """Scan slurm-<array_job_id>_<task>.out files (SLURM's own captured
    stdout/stderr for the driver script, written to the submission
    directory -- separate from each case's own $CASE.out) for a time-limit
    kill, and return the 1-based array task indices that were hit. A
    SLURM-imposed kill happens externally to our script, so the normal
    copy-results-back-then-cleanup path never runs; such a case just shows
    up as "not_run" in assess() with no way to tell it apart from one that
    simply hasn't started yet -- this is how we tell the difference."""
    slurm_log_dir = Path(slurm_log_dir)
    pattern = re.compile(rf"^slurm-{re.escape(array_job_id)}_(\d+)\.out$")
    indices: set[int] = set()
    for path in slurm_log_dir.glob(f"slurm-{array_job_id}_*.out"):
        m = pattern.match(path.name)
        if not m:
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        if TIME_LIMIT_MARKER in text:
            indices.add(int(m.group(1)))
    return indices


def resolve_task_indices_to_case_ids(cases_file: str | Path, indices: set[int]) -> set[str]:
    lines = Path(cases_file).read_text().splitlines()
    result: set[str] = set()
    for idx in indices:
        if 1 <= idx <= len(lines):
            case_id = lines[idx - 1].strip()
            if case_id:
                result.add(case_id)
    return result


def find_timed_out_cases(jobs: list[tuple[str, "str | Path"]], slurm_log_dir: str | Path) -> set[str]:
    """jobs: list of (array_job_id, cases_file) pairs -- one per SLURM
    submission to scan, since each submission may have used a different
    cases file (e.g. the full-batch cases.txt vs a retry's cases_retry.txt)."""
    result: set[str] = set()
    for job_id, cases_file in jobs:
        indices = find_timed_out_task_indices(slurm_log_dir, job_id)
        result |= resolve_task_indices_to_case_ids(cases_file, indices)
    return result


def manifest_case_ids(manifest: str | Path) -> set[str]:
    import json

    with Path(manifest).open() as handle:
        return {json.loads(line)["case_id"] for line in handle if line.strip()}


def attach_reference_cubes(manifest: str | Path, orca_root: str | Path, out: str | Path) -> dict:
    """Write a derived manifest with reference_cube set for every case that
    has succeeded (ORCA TERMINATED NORMALLY + a produced .cube), so
    optimize.py's score_manifest has real DFT densities to compare against.
    Cases that haven't finished yet (or failed) are written through
    unchanged, with no reference_cube -- re-running this after more of the
    SLURM batch finishes just picks up the newly-succeeded ones."""
    import json

    orca_root = Path(orca_root)
    result = assess(orca_root, manifest_case_ids(manifest))
    ok_ids = {e["case_id"] for e in result["by_status"].get("ok", [])}
    total = 0
    attached = 0
    with Path(manifest).open() as handle, Path(out).open("w", newline="\n") as out_handle:
        for line in handle:
            if not line.strip():
                continue
            total += 1
            record = json.loads(line)
            case_id = record["case_id"]
            if case_id in ok_ids:
                cubes = list((orca_root / case_id).glob("*.cube"))
                if cubes:
                    record["reference_cube"] = str(cubes[0])
                    attached += 1
            out_handle.write(json.dumps(record, sort_keys=True) + "\n")
    return {"total": total, "attached": attached}


def stage_local(
    manifest: str | Path,
    out_dir: str | Path,
    out_manifest: str | Path,
    element: str | None = None,
    split: str = "train",
) -> int:
    """Copies each matching row's xyz_path/reference_cube into out_dir and
    writes a derived manifest pointing at the local copies. Used to stage a
    single element's structures onto a compute node's local scratch
    ($TMP_JOBDIR) before optimize-scale, since optimize-scale re-reads both
    files on every trial and doing that against a network /work mount is a
    severe bottleneck across the many trials differential_evolution needs."""
    import json
    import shutil

    out_dir = Path(out_dir)
    rows = []
    with Path(manifest).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            if not record.get("reference_cube") or record.get("split") != split:
                continue
            if element is not None and record.get("element") != element:
                continue
            rows.append(record)

    with Path(out_manifest).open("w", newline="\n") as out_handle:
        for record in rows:
            case_dir = out_dir / record["case_id"]
            case_dir.mkdir(parents=True, exist_ok=True)
            xyz_dst = case_dir / Path(record["xyz_path"]).name
            cube_dst = case_dir / Path(record["reference_cube"]).name
            shutil.copy2(record["xyz_path"], xyz_dst)
            shutil.copy2(record["reference_cube"], cube_dst)
            record["xyz_path"] = str(xyz_dst)
            record["reference_cube"] = str(cube_dst)
            out_handle.write(json.dumps(record, sort_keys=True) + "\n")
    return len(rows)


def assess(orca_root: str | Path, valid_case_ids: set[str] | None = None) -> dict:
    """valid_case_ids, if given (e.g. from manifest_case_ids), restricts
    which directories under orca_root are considered -- stale directories
    left over from a since-regenerated structure (different case_id, same
    slot) are silently ignored instead of needing to be deleted."""
    orca_root = Path(orca_root)
    cases = sorted(p for p in orca_root.iterdir() if p.is_dir())
    if valid_case_ids is not None:
        cases = [p for p in cases if p.name in valid_case_ids]
    counts: dict[str, int] = {}
    by_status: dict[str, list[dict]] = {}
    for case_dir in cases:
        status, detail = classify_case(case_dir)
        counts[status] = counts.get(status, 0) + 1
        by_status.setdefault(status, []).append({"case_id": case_dir.name, "detail": detail})
    return {
        "total_cases": len(cases),
        "counts": counts,
        "by_status": by_status,
    }


def format_report(result: dict, examples_per_status: int = 5) -> str:
    lines = [f"Total cases: {result['total_cases']}"]
    order = [
        "ok",
        "ok_missing_cube",
        "charge_multiplicity_error",
        "input_error",
        "orca_error",
        "unknown_failure",
        "unreadable",
        "not_run",
    ]
    counts = result["counts"]
    for status in order:
        if status in counts:
            lines.append(f"  {status}: {counts[status]}")
    extra = set(counts) - set(order)
    for status in sorted(extra):
        lines.append(f"  {status}: {counts[status]}")

    for status in order:
        entries = result["by_status"].get(status)
        if not entries or status in ("ok",):
            continue
        lines.append(f"\n--- {status} (showing up to {examples_per_status} of {len(entries)}) ---")
        for entry in entries[:examples_per_status]:
            detail = f" :: {entry['detail']}" if entry["detail"] else ""
            lines.append(f"  {entry['case_id']}{detail}")
    return "\n".join(lines)
