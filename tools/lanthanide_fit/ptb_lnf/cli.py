from __future__ import annotations

import argparse
import json
from pathlib import Path

from .baseline import verify
from .cube import compare, read_cube
from .denmat import denmat_to_cube
from .iongate import load_targets
from .optimize import apply_free_result, build_merged_atompara, optimize_free, optimize_scale, score_free_manifest, score_manifest
from .orca import write_orca_inputs, write_slurm
from .status import assess, attach_reference_cubes, find_timed_out_cases, format_report, manifest_case_ids, stage_local
from .structures import generate


def main() -> None:
    parser = argparse.ArgumentParser(prog="ptb_lnf")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("verify-baseline")
    p.add_argument("--repo-root", default=".")

    p = sub.add_parser("generate-structures")
    p.add_argument("--out", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--conformers", type=int, default=2)

    p = sub.add_parser("write-orca")
    p.add_argument("--manifest", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--nprocs", type=int, default=16, help="Value for %%pal nprocs in each .inp; keep in sync with write-slurm --tasks.")
    p.add_argument("--orca-module", default="orca", help="Module name for `module load` in make_cube.sh; keep in sync with write-slurm --orca-module.")

    p = sub.add_parser("write-slurm")
    p.add_argument("--orca-root", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--orca-module", default="orca", help="Module name for `module load`; ORCA is then invoked via $ORCA.")
    p.add_argument("--tasks", type=int, default=16)
    p.add_argument("--time-limit", default="08:00:00", help="SBATCH --time value, e.g. 08:00:00.")
    p.add_argument(
        "--remote-root",
        default=None,
        help=(
            "Path to orca-root as seen from the cluster that will run the script "
            "(e.g. /work/akkleemiss/<project>/orca), for when this is generated on "
            "a different machine than the one the script runs on. Defaults to the "
            "resolved local path of --orca-root."
        ),
    )

    p = sub.add_parser("denmat-to-cube")
    p.add_argument("--denmat", required=True)
    p.add_argument("--grid", required=True)
    p.add_argument("--out", required=True)

    p = sub.add_parser("compare-cubes")
    p.add_argument("--reference", required=True)
    p.add_argument("--candidate", required=True)

    p = sub.add_parser("optimize-scale")
    p.add_argument("--manifest", required=True)
    p.add_argument("--ptb", required=True)
    p.add_argument("--atompara", required=True)
    p.add_argument("--basis", required=True)
    p.add_argument("--workers", type=int, default=None, help="Parallel PTB processes per trial value; defaults to os.cpu_count().")
    p.add_argument(
        "--element",
        default=None,
        help="Fit only this element's block (e.g. Pr), scoring only against its own train structures. "
        "The fit fully decomposes per element, so this can run independently/concurrently per element.",
    )
    p.add_argument(
        "--bounds",
        nargs=2,
        type=float,
        default=[0.1, 8.0],
        metavar=("LOW", "HIGH"),
        help="Search bounds applied to all group scale factors, e.g. --bounds 0.4 1.3. If a fitted "
        "result lands within 2%% of either bound, a warning is printed (the true optimum likely lies outside).",
    )
    p.add_argument(
        "--popsize",
        type=int,
        default=15,
        help="differential_evolution popsize (population = popsize * number of groups per generation). "
        "Lower this to control wall-clock cost.",
    )
    p.add_argument(
        "--maxiter",
        type=int,
        default=100,
        help="differential_evolution maxiter (generations).",
    )

    p = sub.add_parser("optimize-free")
    p.add_argument("--manifest", required=True)
    p.add_argument("--ptb", required=True)
    p.add_argument("--atompara", required=True, help="Starting template -- should already carry this element's fitted group scales (e.g. from merge-atompara).")
    p.add_argument("--basis", required=True)
    p.add_argument("--element", required=True, help="Fit only this element's block, scoring only against its own train structures.")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument(
        "--bounds", nargs=2, type=float, default=[0.5, 2.0], metavar=("LOW", "HIGH"),
        help="Per-value multiplicative bounds relative to --atompara's already-fitted starting point (a local refine, not a from-scratch search).",
    )
    p.add_argument("--passes", type=int, default=1, help="Number of full coordinate-descent sweeps over every free parameter.")
    p.add_argument("--evals-per-param", type=int, default=8, help="Bounded 1-D minimize_scalar budget per parameter per sweep.")
    p.add_argument("--checkpoint", default=None, help="Path to write a resumable per-parameter checkpoint after each parameter (recommended -- these runs are long and SLURM time-limit kills are common on this cluster).")
    p.add_argument("--resume-from", default=None, help="Path to a prior --checkpoint file to resume from, if it exists.")
    p.add_argument("--exclude-positions", default=None, help="JSON file with a positions array ([[row, col], ...]) to hold at template scale 1.0. Use only density-validated inactive coordinates.")
    p.add_argument("--min-sweep-relative-improvement", type=float, default=None, help="Stop after --patience completed sweeps whose relative training-score improvement is below this threshold.")
    p.add_argument("--patience", type=int, default=1, help="Number of consecutive low-improvement sweeps allowed before early stop.")
    p.add_argument("--verify-resume", action="store_true", help="Score the exact resumed vector before searching; abort if it is non-finite.")
    p.add_argument("--relative-step", type=float, default=None, help="Search each coordinate only within +/- this relative fraction of its resumed/current scale, clipped to --bounds.")
    p.add_argument("--ion-targets", default=None, help="Free-ion l-population windows JSON (see iongate.py); trials outside them score 1 + violation without molecular scoring.")

    p = sub.add_parser("evaluate-free")
    p.add_argument("--manifest", required=True)
    p.add_argument("--ptb", required=True)
    p.add_argument("--atompara", required=True, help="Same group-scaled template the optimize-free result's scales are relative to.")
    p.add_argument("--basis", required=True)
    p.add_argument("--element", required=True)
    p.add_argument("--result", required=True, help="Path to an optimize-free result JSON to evaluate.")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--split", default="validation")

    p = sub.add_parser("apply-free")
    p.add_argument("--template", required=True, help="Group-scaled merged .atompara the optimize-free result's scales are relative to.")
    p.add_argument("--out", required=True)
    p.add_argument("--element", required=True)
    p.add_argument("--result", required=True, help="Path to an optimize-free result JSON (or checkpoint file, same schema).")

    p = sub.add_parser("merge-atompara")
    p.add_argument("--template", required=True)
    p.add_argument("--out", required=True)
    p.add_argument(
        "--scales",
        required=True,
        help="Path to a JSON file mapping element symbol -> {group: fitted scale}, e.g. "
        "{\"Pr\": {\"exponent_scaling\": 1.05, \"shell_exponents\": 0.85, "
        "\"shell_config\": 1.0, \"shell_response\": 0.95}}. Matches the \"scales\" field of "
        "each optimize-scale result. Any group not present (e.g. \"energy\", which has no "
        "effect on the output density and is excluded from the fit) is left unscaled.",
    )

    p = sub.add_parser("evaluate-scale")
    p.add_argument("--manifest", required=True)
    p.add_argument("--ptb", required=True)
    p.add_argument("--atompara", required=True)
    p.add_argument("--basis", required=True)
    p.add_argument("--element", required=True)
    p.add_argument(
        "--result",
        required=True,
        help="Path to an optimize-scale result JSON (reads its \"scales\" field) to evaluate.",
    )
    p.add_argument("--workers", type=int, default=None)
    p.add_argument(
        "--split",
        default="validation",
        help="Which manifest split to score against (default: validation, the held-out ~20%% "
        "never used during fitting -- compare its score against the result's train \"score\" "
        "to check for overfitting).",
    )

    p = sub.add_parser("assess-orca")
    p.add_argument("--orca-root", required=True)
    p.add_argument("--json-out", default=None, help="Optional path to write the full per-case breakdown as JSON.")
    p.add_argument("--examples", type=int, default=5, help="Number of example case_ids to print per failure category.")
    p.add_argument(
        "--manifest",
        default=None,
        help="If given, only consider directories whose name is a case_id in this manifest (ignores stale leftover directories from regenerated structures).",
    )

    p = sub.add_parser("retry-orca")
    p.add_argument("--orca-root", required=True)
    p.add_argument("--out", required=True, help="Path to write the retry SLURM script, e.g. work/orca_array_retry.slurm.")
    p.add_argument("--orca-module", default="orca")
    p.add_argument("--tasks", type=int, default=16)
    p.add_argument("--time-limit", default="08:00:00", help="SBATCH --time value, e.g. 08:00:00.")
    p.add_argument("--remote-root", default=None)
    p.add_argument(
        "--cases-filename",
        default="cases_retry.txt",
        help="Filename (inside orca-root) to write the retry case list to; kept separate from cases.txt.",
    )
    p.add_argument(
        "--manifest",
        default=None,
        help="If given, only consider directories whose name is a case_id in this manifest (ignores stale leftover directories from regenerated structures).",
    )

    p = sub.add_parser("retry-timeouts")
    p.add_argument("--orca-root", required=True)
    p.add_argument("--manifest", required=True, help="Used both to cross-check current status and to regenerate .inp for retried cases with more cores.")
    p.add_argument("--slurm-log-dir", required=True, help="Directory containing slurm-<jobid>_<task>.out files (the sbatch submission dir).")
    p.add_argument(
        "--job",
        action="append",
        required=True,
        metavar="JOBID:CASES_FILENAME",
        help=(
            "One SLURM array job to scan for time-limit kills, as jobid:cases_filename "
            "(cases_filename is relative to --orca-root, e.g. 255802:cases.txt). Repeatable "
            "for multiple submissions (e.g. the full batch and a prior retry)."
        ),
    )
    p.add_argument("--out", required=True, help="Path to write the timeout-retry SLURM script.")
    p.add_argument("--orca-module", default="orca")
    p.add_argument("--tasks", type=int, default=8, help="Bumped core count (both %%pal nprocs and --ntasks) for these slow cases.")
    p.add_argument("--time-limit", default="16:00:00", help="Bumped SBATCH --time value for these slow cases.")
    p.add_argument("--remote-root", default=None)
    p.add_argument("--cases-filename", default="cases_timeout_retry.txt")

    p = sub.add_parser("stage-local")
    p.add_argument("--manifest", required=True)
    p.add_argument("--out-dir", required=True, help="Directory to copy each matching row's xyz/cube into (e.g. $TMP_JOBDIR/data).")
    p.add_argument("--out-manifest", required=True, help="Path to write the derived manifest pointing at the local copies.")
    p.add_argument("--element", default=None, help="Only stage this element's rows.")
    p.add_argument("--split", default="train")

    p = sub.add_parser("attach-cubes")
    p.add_argument("--manifest", required=True)
    p.add_argument("--orca-root", required=True)
    p.add_argument("--out", required=True, help="Path to write the derived manifest with reference_cube set.")

    args = parser.parse_args()
    if args.command == "verify-baseline":
        print(json.dumps(verify(args.repo_root), indent=2, sort_keys=True))
    elif args.command == "generate-structures":
        records = generate(args.out, args.manifest, args.conformers)
        print(json.dumps({"records": len(records), "manifest": args.manifest}, indent=2))
    elif args.command == "write-orca":
        write_orca_inputs(args.manifest, args.out, args.nprocs, args.orca_module)
    elif args.command == "write-slurm":
        write_slurm(args.orca_root, args.out, args.orca_module, args.tasks, args.remote_root, time_limit=args.time_limit)
    elif args.command == "denmat-to-cube":
        denmat_to_cube(args.denmat, args.grid, args.out)
    elif args.command == "compare-cubes":
        print(json.dumps(compare(read_cube(args.reference), read_cube(args.candidate)), indent=2, sort_keys=True))
    elif args.command == "optimize-scale":
        result = optimize_scale(
            args.manifest, args.ptb, args.atompara, args.basis,
            workers=args.workers, element=args.element, bounds=tuple(args.bounds),
            popsize=args.popsize, maxiter=args.maxiter,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    elif args.command == "optimize-free":
        excluded = set()
        if args.exclude_positions:
            with open(args.exclude_positions) as handle:
                exclusion_data = json.load(handle)
            raw_positions = exclusion_data.get("positions", exclusion_data.get("exact_no_response_positions", exclusion_data)) if isinstance(exclusion_data, dict) else exclusion_data
            excluded = {tuple(map(int, position)) for position in raw_positions}
        result = optimize_free(
            args.manifest, args.ptb, args.atompara, args.basis, args.element,
            workers=args.workers, bounds=tuple(args.bounds), passes=args.passes,
            evals_per_param=args.evals_per_param, checkpoint=args.checkpoint, resume_from=args.resume_from,
            exclude_positions=excluded, min_sweep_relative_improvement=args.min_sweep_relative_improvement,
            patience=args.patience, verify_resume=args.verify_resume, relative_step=args.relative_step,
            ion_targets=load_targets(args.ion_targets) if args.ion_targets else None,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    elif args.command == "evaluate-free":
        with open(args.result) as handle:
            fitted = json.load(handle)
        score = score_free_manifest(
            args.manifest, args.ptb, args.atompara, args.basis, args.element, fitted,
            workers=args.workers, split=args.split,
        )
        print(json.dumps({
            "element": args.element,
            "split": args.split,
            "score": score,
            "train_score": fitted.get("score"),
        }, indent=2, sort_keys=True))
    elif args.command == "apply-free":
        with open(args.result) as handle:
            result = json.load(handle)
        apply_free_result(args.template, args.out, args.element, result)
        print(f"Wrote {args.out} with {args.element}'s free-fit result applied ({len(result['positions'])} parameter(s), score={result['score']:.6f})")
    elif args.command == "merge-atompara":
        with open(args.scales) as handle:
            element_scales = json.load(handle)
        build_merged_atompara(args.template, args.out, element_scales)
        print(f"Wrote {args.out} with {len(element_scales)} element-specific scale(s): {element_scales}")
    elif args.command == "evaluate-scale":
        with open(args.result) as handle:
            fitted = json.load(handle)
        group_scales = fitted["scales"]
        score = score_manifest(
            args.manifest, args.ptb, args.atompara, args.basis, group_scales,
            workers=args.workers, element=args.element, split=args.split,
        )
        print(json.dumps({
            "element": args.element,
            "split": args.split,
            "score": score,
            "train_score": fitted.get("score"),
            "scales": group_scales,
        }, indent=2, sort_keys=True))
    elif args.command == "assess-orca":
        valid_ids = manifest_case_ids(args.manifest) if args.manifest else None
        result = assess(args.orca_root, valid_ids)
        print(format_report(result, args.examples))
        if args.json_out:
            with open(args.json_out, "w", newline="\n") as handle:
                json.dump(result, handle, indent=2, sort_keys=True)
    elif args.command == "retry-orca":
        valid_ids = manifest_case_ids(args.manifest) if args.manifest else None
        result = assess(args.orca_root, valid_ids)
        retry_ids = sorted(
            entry["case_id"]
            for status, entries in result["by_status"].items()
            if status != "ok"
            for entry in entries
        )
        if not retry_ids:
            print("Nothing to retry: every case is already ok.")
        else:
            write_slurm(
                args.orca_root, args.out, args.orca_module, args.tasks, args.remote_root,
                case_ids=retry_ids, cases_filename=args.cases_filename, job_name="ptb-ln-orca-retry",
                time_limit=args.time_limit,
            )
            print(f"Wrote retry script for {len(retry_ids)} case(s) (of {result['total_cases']} total) to {args.out}")
    elif args.command == "retry-timeouts":
        jobs = []
        for spec in args.job:
            job_id, cases_filename = spec.split(":", 1)
            jobs.append((job_id, Path(args.orca_root) / cases_filename))
        timed_out = find_timed_out_cases(jobs, args.slurm_log_dir)
        valid_ids = manifest_case_ids(args.manifest)
        result = assess(args.orca_root, valid_ids)
        ok_ids = {e["case_id"] for e in result["by_status"].get("ok", [])}
        retry_ids = sorted(c for c in timed_out if c in valid_ids and c not in ok_ids)
        print(f"Found {len(timed_out)} case(s) killed by time limit across {len(jobs)} job(s); {len(retry_ids)} still need retry (rest already succeeded since).")
        if not retry_ids:
            print("Nothing to retry.")
        else:
            write_orca_inputs(args.manifest, args.orca_root, nprocs=args.tasks, orca_module=args.orca_module, case_ids=set(retry_ids))
            write_slurm(
                args.orca_root, args.out, args.orca_module, args.tasks, args.remote_root,
                case_ids=retry_ids, cases_filename=args.cases_filename, job_name="ptb-ln-orca-timeout-retry",
                time_limit=args.time_limit,
            )
            print(f"Regenerated .inp for {len(retry_ids)} case(s) with nprocs={args.tasks}; wrote {args.out} (--tasks={args.tasks}, --time={args.time_limit}).")
    elif args.command == "stage-local":
        n = stage_local(args.manifest, args.out_dir, args.out_manifest, element=args.element, split=args.split)
        print(f"Staged {n} row(s) into {args.out_dir}; wrote {args.out_manifest}")
    elif args.command == "attach-cubes":
        result = attach_reference_cubes(args.manifest, args.orca_root, args.out)
        print(f"Attached reference_cube to {result['attached']} of {result['total']} record(s); wrote {args.out}")


if __name__ == "__main__":
    main()
