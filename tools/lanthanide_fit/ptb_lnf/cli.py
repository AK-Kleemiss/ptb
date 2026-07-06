from __future__ import annotations

import argparse
import json

from .baseline import verify
from .cube import compare, read_cube
from .denmat import denmat_to_cube
from .optimize import optimize_scale
from .orca import write_orca_inputs, write_slurm
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

    p = sub.add_parser("write-slurm")
    p.add_argument("--orca-root", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--orca-bin", default="orca")
    p.add_argument("--tasks", type=int, default=16)

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

    args = parser.parse_args()
    if args.command == "verify-baseline":
        print(json.dumps(verify(args.repo_root), indent=2, sort_keys=True))
    elif args.command == "generate-structures":
        records = generate(args.out, args.manifest, args.conformers)
        print(json.dumps({"records": len(records), "manifest": args.manifest}, indent=2))
    elif args.command == "write-orca":
        write_orca_inputs(args.manifest, args.out)
    elif args.command == "write-slurm":
        write_slurm(args.orca_root, args.out, args.orca_bin, args.tasks)
    elif args.command == "denmat-to-cube":
        denmat_to_cube(args.denmat, args.grid, args.out)
    elif args.command == "compare-cubes":
        print(json.dumps(compare(read_cube(args.reference), read_cube(args.candidate)), indent=2, sort_keys=True))
    elif args.command == "optimize-scale":
        print(json.dumps(optimize_scale(args.manifest, args.ptb, args.atompara, args.basis), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
