from __future__ import annotations

from pathlib import Path
import json


def read_manifest(path: str | Path) -> list[dict]:
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def write_orca_inputs(
    manifest: str | Path,
    out: str | Path,
    nprocs: int = 16,
    orca_module: str = "orca",
    case_ids: set[str] | None = None,
    extra_keywords: str = "",
) -> None:
    out = Path(out)
    keyword_line = "! wB97X-3c TightSCF" + (f" {extra_keywords}" if extra_keywords else "")
    for record in read_manifest(manifest):
        if case_ids is not None and record["case_id"] not in case_ids:
            continue
        case_dir = out / record["case_id"]
        case_dir.mkdir(parents=True, exist_ok=True)
        xyz_name = Path(record["xyz_path"]).name
        inp = case_dir / f"{record['case_id']}.inp"
        with inp.open("w", newline="\n") as handle:
            handle.write("\n".join([
                keyword_line,
                "%scf",
                "  MaxIter 500",
                "end",
                f"%pal nprocs {nprocs} end",
                f"* xyzfile {record['charge']} {record['spin_multiplicity']} {xyz_name}",
                "",
            ]))
        xyz_src = Path(record["xyz_path"])
        xyz_dst = case_dir / xyz_name
        if xyz_src.resolve() != xyz_dst.resolve():
            with xyz_dst.open("w", newline="\n") as handle:
                handle.write(xyz_src.read_text())
        with (case_dir / "make_cube.sh").open("w", newline="\n") as handle:
            # orca_plot must be invoked from the same directory as the orca
            # binary that produced the .gbw (a bare `orca_plot` on PATH can
            # resolve to a mismatched build), so derive it from $ORCA set by
            # the cluster's orca module rather than relying on PATH alone.
            handle.write("\n".join([
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                f"module load {orca_module}",
                # orca_plot's interactive menu: 1=enter type of plot,
                # 2=electron density; 5=select output format, 7=Gaussian
                # cube; 11=generate the plot; 12=exit. Without the 1/2 pair
                # first, it silently defaults to plotting an MO instead of
                # the density.
                f'"$(dirname "$ORCA")/orca_plot" {record["case_id"]}.gbw -i <<\'EOF\'',
                "1",
                "2",
                "5",
                "7",
                "11",
                "12",
                "EOF",
                "",
            ]))


def write_slurm(
    orca_root: str | Path,
    out: str | Path,
    orca_module: str = "orca",
    tasks: int = 16,
    remote_root: str | None = None,
    case_ids: list[str] | None = None,
    cases_filename: str = "cases.txt",
    job_name: str = "ptb-ln-orca",
    time_limit: str = "08:00:00",
) -> None:
    orca_root = Path(orca_root)
    # case_ids lets a caller (e.g. retry-orca) target a subset of cases
    # instead of every directory under orca_root, so already-succeeded
    # cases aren't resubmitted. cases_filename keeps a retry's case list
    # separate from the full-batch cases.txt.
    cases = list(case_ids) if case_ids is not None else sorted(path.name for path in orca_root.iterdir() if path.is_dir())
    case_list = orca_root / cases_filename
    with case_list.open("w", newline="\n") as handle:
        handle.write("\n".join(cases) + "\n")
    # remote_root lets this be generated on a machine (e.g. Windows) where
    # orca_root's local path is meaningless on the cluster that will run the
    # script; pass the path as seen from the cluster (e.g. a Linux mount of
    # the same share) instead of resolving orca_root locally.
    embedded_root = remote_root if remote_root is not None else str(orca_root.resolve())
    with Path(out).open("w", newline="\n") as handle:
        handle.write("\n".join([
            "#!/usr/bin/env bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --array=1-{len(cases)}",
            # --nodes=1 required, not just a default: this task stages its
            # case directory onto node-local scratch ($TMP_JOBDIR below) and
            # runs a single ORCA process there with its own %pal nprocs (no
            # MPI hostfile spanning nodes) -- any additional node's cores
            # would sit allocated but unused. (2026-07-28: user asked to
            # allow multi-node scheduling for jobs that don't stage to local
            # scratch, to pack onto fragmented/"mixed"-state nodes more
            # easily; this job is the excluded case.)
            "#SBATCH --nodes=1",
            f"#SBATCH --ntasks={tasks}",
            "#SBATCH --mem=16G",
            f"#SBATCH --time={time_limit}",
            "#SBATCH --output=slurm-%A_%a.out",
            "set -uo pipefail",
            f"module load {orca_module}",
            f"ORCA_ROOT={embedded_root}",
            f'CASE=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" "$ORCA_ROOT/{cases_filename}")',
            'CASE_DIR="$ORCA_ROOT/$CASE"',
            # A silently-failed `cd` here (e.g. empty $CASE from a
            # cases-file/array-size mismatch) would otherwise leave every
            # subsequent command running in whatever directory the job
            # started in -- $ORCA_ROOT's parent -- scattering that case's
            # entire output next to the .slurm scripts instead of its own
            # folder. Fail loudly and immediately instead.
            'if [ -z "$CASE" ] || [ ! -d "$CASE_DIR" ]; then',
            '  echo "ERROR: could not resolve case for array task $SLURM_ARRAY_TASK_ID (CASE=\'$CASE\' CASE_DIR=\'$CASE_DIR\')" >&2',
            "  exit 1",
            "fi",
            # Diagnostic instrumentation for an intermittent, still-unexplained
            # bug where a full duplicate copy of a case's output has been
            # found sitting in $ORCA_ROOT's parent directory alongside the
            # correctly-placed copy in $CASE_DIR, with no error printed by
            # any of the checks below. Cheap to leave in permanently.
            'echo "DIAG task=$SLURM_ARRAY_TASK_ID jobid=$SLURM_JOB_ID host=$(hostname) CASE=\'$CASE\' CASE_DIR=\'$CASE_DIR\' TMP_JOBDIR_pre=\'${TMP_JOBDIR:-<unset>}\' pwd=\'$(pwd)\'" >&2',
            # ORCA does heavy scratch I/O; run on the cluster's own
            # per-job scratch ($TMP_JOBDIR) instead of the /work network
            # mount. SLURM creates and tears this down itself regardless of
            # how the job ends (success, failure, timeout, scancel), so
            # unlike a manually mktemp'd directory it can never be orphaned
            # on the node -- we still copy results back to $CASE_DIR
            # ourselves since nothing in $TMP_JOBDIR survives the job.
            'cd "$CASE_DIR" || { echo "ERROR: cd to CASE_DIR=\'$CASE_DIR\' failed" >&2; exit 1; }',
            'echo "DIAG after-cd-CASE_DIR pwd=\'$(pwd)\'" >&2',
            'if [ -z "${TMP_JOBDIR:-}" ] || [ ! -d "$TMP_JOBDIR" ]; then',
            '  echo "ERROR: TMP_JOBDIR is unset or missing (\'${TMP_JOBDIR:-}\')" >&2',
            "  exit 1",
            "fi",
            "rsync -a --exclude='*.out' --exclude='*.cube' --exclude='*.gbw' * \"$TMP_JOBDIR\"",
            'cd "$TMP_JOBDIR" || { echo "ERROR: cd to TMP_JOBDIR=\'$TMP_JOBDIR\' failed" >&2; exit 1; }',
            'echo "DIAG after-cd-TMP_JOBDIR pwd=\'$(pwd)\' TMP_JOBDIR=\'$TMP_JOBDIR\' realpath_TMP_JOBDIR=\'$(realpath "$TMP_JOBDIR" 2>/dev/null)\'" >&2',
            # Best-effort salvage of a partial .out if SLURM SIGTERMs this
            # job (time-limit/scancel) before it finishes normally.
            'trap \'cp -f "$TMP_JOBDIR"/*.out "$CASE_DIR/" 2>/dev/null; exit 143\' TERM',
            "status=0",
            # ORCA must be called via its full path (not a bare name off
            # PATH) so it can locate co-located helper binaries; $ORCA is
            # set by the module load above.
            '"$ORCA" "$CASE.inp" > "$CASE.out" 2>&1 || status=$?',
            'if [ "$status" -eq 0 ]; then',
            "  bash make_cube.sh || status=$?",
            "fi",
            'echo "DIAG before-cp-back pwd=\'$(pwd)\' CASE_DIR=\'$CASE_DIR\' TMP_JOBDIR=\'$TMP_JOBDIR\' nfiles=$(ls -1 "$TMP_JOBDIR" | wc -l)" >&2',
            'cp -f "$TMP_JOBDIR"/* "$CASE_DIR"/',
            'echo "DIAG after-cp-back status=$status" >&2',
            'if [ "$status" -ne 0 ]; then',
            '  echo "Case $CASE failed (exit $status) on $(hostname)" >&2',
            "fi",
            'exit "$status"',
            "",
        ]))
