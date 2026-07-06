from __future__ import annotations

from pathlib import Path
import json


def read_manifest(path: str | Path) -> list[dict]:
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def write_orca_inputs(manifest: str | Path, out: str | Path) -> None:
    out = Path(out)
    for record in read_manifest(manifest):
        case_dir = out / record["case_id"]
        case_dir.mkdir(parents=True, exist_ok=True)
        xyz_name = Path(record["xyz_path"]).name
        inp = case_dir / f"{record['case_id']}.inp"
        with inp.open("w", newline="\n") as handle:
            handle.write("\n".join([
                "! wB97X-3c TightSCF Grid5 FinalGrid6 RIJCOSX",
                "%scf",
                "  MaxIter 500",
                "end",
                "%pal nprocs 16 end",
                f"* xyzfile {record['charge']} {record['spin_multiplicity']} {xyz_name}",
                "",
            ]))
        xyz_src = Path(record["xyz_path"])
        xyz_dst = case_dir / xyz_name
        if xyz_src.resolve() != xyz_dst.resolve():
            with xyz_dst.open("w", newline="\n") as handle:
                handle.write(xyz_src.read_text())
        with (case_dir / "make_cube.sh").open("w", newline="\n") as handle:
            handle.write("\n".join([
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                f"orca_plot {record['case_id']}.gbw -i <<'EOF'",
                "5",
                "7",
                "10",
                "11",
                "EOF",
                "",
            ]))


def write_slurm(orca_root: str | Path, out: str | Path, orca_bin: str = "orca", tasks: int = 16) -> None:
    orca_root = Path(orca_root)
    cases = sorted(path.name for path in orca_root.iterdir() if path.is_dir())
    case_list = orca_root / "cases.txt"
    with case_list.open("w", newline="\n") as handle:
        handle.write("\n".join(cases) + "\n")
    with Path(out).open("w", newline="\n") as handle:
        handle.write("\n".join([
            "#!/usr/bin/env bash",
            "#SBATCH --job-name=ptb-ln-orca",
            f"#SBATCH --array=1-{len(cases)}",
            f"#SBATCH --ntasks={tasks}",
            "#SBATCH --mem=16G",
            "#SBATCH --time=04:00:00",
            "#SBATCH --output=slurm-%A_%a.out",
            "set -euo pipefail",
            f"ORCA_ROOT={orca_root.resolve()}",
            'CASE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ORCA_ROOT/cases.txt")',
            'cd "$ORCA_ROOT/$CASE"',
            f"{orca_bin} \"$CASE.inp\" > \"$CASE.out\"",
            "bash make_cube.sh",
            "",
        ]))
