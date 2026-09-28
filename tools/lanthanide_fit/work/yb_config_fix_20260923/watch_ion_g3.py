"""Watch job 596350, reproduce its candidate, and write a gated verdict.

This does not edit PTB defaults or commit a candidate.
"""
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
JOB = "596350"
HOST = "florian@AKL007"
KEY = Path("C:/Users/florian/.ssh/id_rsa2")
REMOTE = "/work/akkleemiss/florian/lnf_orca/yb_config_fix_20260923"
DEST = HERE / "ion_g3_size"
PTB = ROOT / "build/ptb_windows_ifx_yb_ion_g2_20260928.exe"
BAS = ROOT / ".basis_vDZP"
TARGETS = HERE / "ion_targets_g1.json"
PY = [sys.executable]
SSH = ["ssh", "-i", str(KEY), "-o", "BatchMode=yes", "-o", "IdentitiesOnly=yes",
       "-o", "ForwardAgent=no", "-o", "ConnectTimeout=15", HOST]
SCP = ["scp", "-r", "-i", str(KEY), "-o", "BatchMode=yes", "-o", "IdentitiesOnly=yes",
       "-o", "ForwardAgent=no", "-o", "ConnectTimeout=15"]
FLAGS = getattr(subprocess, "CREATE_NO_WINDOW", 0)


def log(message):
    print(datetime.now(timezone.utc).isoformat(), message, flush=True)


def call(args, *, env=None, check=True):
    done = subprocess.run(args, cwd=ROOT, env=env, capture_output=True, text=True,
                          timeout=300, creationflags=FLAGS)
    if check and done.returncode:
        raise RuntimeError(f"command {args[0]} failed {done.returncode}: {done.stdout[-500:]} {done.stderr[-500:]}")
    return done.stdout + done.stderr


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def score(split, candidate):
    args = PY + [str(HERE / "score_local.py")]
    if split == "validation":
        args += ["--split", split]
    output = call(args + [str(candidate)])
    (DEST / f"windows_{split}.log").write_text(output)
    hit = re.search(rf"\b{split} ([0-9.]+)\b", output)
    if not hit:
        raise RuntimeError(f"no {split} score in Windows output: {output[-500:]}")
    return float(hit.group(1))


def validate():
    candidate = DEST / "final_candidate.atompara"
    result = json.loads((DEST / "result.json").read_text())
    validation = json.loads((DEST / "validation.json").read_text())
    ions = json.loads((DEST / "ions.json").read_text())
    base = HERE / "ion_g2/final_candidate.atompara"
    def rows(path):
        lines = path.read_text().splitlines()
        i = next(i for i, line in enumerate(lines) if line.strip() == "70")
        return lines[i+1:i+16]
    a, b = rows(base), rows(candidate)
    frozen = all(a[r] == b[r] for r in (6, 11)) and a[13].split()[10] == b[13].split()[10]
    if not frozen:
        raise RuntimeError("radial parameters changed despite exclusions")
    sys.path.insert(0, str(ROOT / "tools/lanthanide_fit"))
    from ptb_lnf.iongate import gate_violation, load_targets
    violation, pops = gate_violation(PTB, candidate, BAS, load_targets(TARGETS))
    train_local = score("train", candidate)
    val_local = score("validation", candidate)
    env = os.environ.copy()
    env.update(YB_SHELL_CANDIDATE=str(candidate),
               YB_SHELL_OUT=str(DEST / "shell_sizes"))
    call(PY + [str(HERE / "measure_shell_sizes.py")], env=env)
    env.update(YB_PTB=str(PTB), YB_CANDIDATE=str(candidate),
               YB_VALIDATION_OUT=str(DEST / "independent_windows_nosphera2_validation"))
    call(PY + [str(HERE.parent / "weekend_20260911/validate_yb_dshell7_onsite_20260914.py")], env=env)
    sizes = json.loads((DEST / "shell_sizes/result.json").read_text())
    old_sizes = json.loads((HERE / "shell_sizes_20260928/result.json").read_text())
    def shell_radii(data):
        return {(x["element"], x["charge"], s): x["shells"][s][1]["ao_extent"]["rms_radius"]
                for x in data["records"] for s in ("6", "7", "8", "9")}
    new_radii, old_radii = shell_radii(sizes), shell_radii(old_sizes)
    max_shell_delta = max(abs(new_radii[k] - old_radii[k]) for k in old_radii)
    exports = json.loads((DEST / "independent_windows_nosphera2_validation/result.json").read_text())
    old_exports = json.loads((HERE / "ion_g2/independent_windows_nosphera2_validation/result.json").read_text())
    old_radius = {x["charge"]: x["independent_export"][1]["rms_radius"] for x in old_exports["states"]}
    new_radius = {x["charge"]: x["independent_export"][1]["rms_radius"] for x in exports["states"]}
    counts_ok = all(abs(run["counts"][key] - x["expected"]) < 1e-3
                    for x in exports["states"] for run in x["independent_export"] for key in ("20", "80"))
    converged = all(abs(x["independent_export"][0]["counts"][key] - x["independent_export"][1]["counts"][key]) < 1e-5
                    for x in exports["states"] for key in ("20", "80"))
    checks = {
        "train_improved": result["score"] < 0.019426166409816967,
        "train_reproduced": abs(train_local - result["score"]) < 2e-6,
        "validation_reproduced": abs(val_local - validation["score"]) < 2e-6,
        "validation_no_regression": validation["score"] <= 0.018851480493416128 * 1.01,
        "radial_parameters_frozen": frozen,
        "ao_radii_unchanged": max_shell_delta < 1e-4,
        "free_ion_gate": violation == 0.0 and ions["violation"] == 0.0,
        "electron_counts": counts_ok and converged,
        "total_atomic_radii": all(new_radius[q] <= old_radius[q] * 1.01 for q in old_radius),
    }
    verdict = {"created_utc": datetime.now(timezone.utc).isoformat(),
               "job": int(JOB), "candidate_sha256": sha(candidate),
               "ptb_sha256": sha(PTB), "basis_sha256": sha(BAS),
               "train_cluster": result["score"], "validation_cluster": validation["score"],
               "train_windows": train_local, "validation_windows": val_local,
               "ion_violation_windows": violation, "ion_populations_windows": pops,
               "max_ao_radius_change_bohr": max_shell_delta,
               "atomic_radii_before": old_radius, "atomic_radii_after": new_radius,
               "checks": checks, "status": "passed_for_review" if all(checks.values()) else "rejected"}
    (DEST / "verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    log(verdict["status"] + " " + str(checks))


def main():
    while True:
        output = call(SSH + [f"sacct -j {JOB} -o JobID,State -n -P | head -1; tail -4 {REMOTE}/yb_ion_g3-{JOB}.out 2>/dev/null"], check=False)
        if "DONE ion_g3_size" in output:
            log("job completed, fetching results")
            call(SCP + [f"{HOST}:{REMOTE}/results/ion_g3_size", str(HERE)])
            validate()
            return
        if re.search(rf"^{JOB}\|(FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)\b", output, re.M):
            raise RuntimeError("job ended without complete result: " + output)
        log(output.strip().replace("\n", " | ")[-300:])
        time.sleep(900)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        DEST.mkdir(exist_ok=True)
        (DEST / "watcher_error.txt").write_text(repr(exc) + "\n")
        log("ERROR " + repr(exc))
        raise
