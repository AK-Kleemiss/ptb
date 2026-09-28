#!/usr/bin/env bash
# Continue ion_g2 while freezing rows 6 and 11 and shell_resp(11,1) at
# (13,10), which control the initial and responsive radial scales.
# Free-ion populations remain a scored constraint.
# The 717 training cases select; 169 validation cases are an overfit check.
#SBATCH --job-name=yb-ion-g3-size
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=92
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=2-00:00:00
#SBATCH --qos=ent20,overflow
#SBATCH --output=yb_ion_g3-%j.out
set -euo pipefail
ROOT=/work/akkleemiss/florian/lnf_orca
D7=$ROOT/yb_d7_constraint_20260916
EXP=$ROOT/yb_basis_d2_20260921
FIX=$ROOT/yb_config_fix_20260923
TAG=ion_g3_size
RUN=${SLURM_JOB_ID:-manual_$(date +%Y%m%dT%H%M%S)}
WORKERS=${SLURM_NTASKS:-92}
source "$ROOT/venv/bin/activate"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
TASK=$(mktemp -d "/tmp/yb_ion_g3.${RUN}.XXXXXX")
finish() {
  if [ -n "${PUBLISHER:-}" ]; then
    kill "$PUBLISHER" 2>/dev/null || true
    wait "$PUBLISHER" 2>/dev/null || true
  fi
  [ ! -s "$TASK/checkpoint.json" ] || cp "$TASK/checkpoint.json" "$FIX/results/${TAG}.checkpoint.publish.${RUN}.json"
  rm -rf "$TASK"
}
trap finish EXIT
cp -a "$ROOT/ptb_lnf_iongate_20260923/ptb_lnf" "$TASK/ptb_lnf"
cp "$FIX/results/ion_g2/final_candidate.atompara" "$TASK/template.atompara"
cp "$FIX/ion_targets.json" "$TASK/ion_targets.json"
cp "$D7/manifest_source.jsonl" "$TASK/manifest_source.jsonl"
cp "$ROOT/build/ptb_valel_fix_ifx" "$TASK/ptb"
cp "$EXP/basis_vDZP_yb_d2_20260921" "$TASK/basis_vDZP"
cd "$TASK"
test "$(sha256sum template.atompara | cut -d' ' -f1)" = 6eb12794d40838ce9b96050df769be7daa1c6e311d941b5dd6b048864ea51c39
test "$(sha256sum basis_vDZP | cut -d' ' -f1)" = 28fce65c2dc48cd47c57ef4984c0ac15135b47ff90ed39ec1f43fbc909eb964c
mkdir -p "$FIX/results/$TAG"
python -m ptb_lnf.cli stage-local --manifest manifest_source.jsonl --out-dir data_train --out-manifest train.jsonl --element Yb
python -m ptb_lnf.cli stage-local --manifest manifest_source.jsonl --out-dir data_validation --out-manifest validation.jsonl --element Yb --split validation
python - <<'PY'
import json
excluded=[[r,c] for r in (6,11) for c in range(13)]+[[13,10]]
positions=[[r,c] for r in range(6,15) for c in range(13) if [r,c] not in excluded]
json.dump({'positions': excluded}, open('excluded_positions.json','w'))
json.dump({'element':'Yb','positions':positions,'value_scales':[1.0]*len(positions),'score':None,'ndim':len(positions)}, open('identity.json','w'))
PY
python -m ptb_lnf.cli evaluate-free --manifest train.jsonl --ptb "$TASK/ptb" --atompara "$TASK/template.atompara" --basis "$TASK/basis_vDZP" --element Yb --result identity.json --workers "$WORKERS" --split train > baseline_train.json
cp baseline_train.json "$FIX/results/$TAG/"
python - <<'PY'
import json
s=json.load(open('baseline_train.json'))['score']
print('START train score',s,flush=True)
assert abs(s-0.019426166409816967)<2e-6
PY
(while true; do [ ! -s checkpoint.json ] || cp checkpoint.json "$FIX/results/${TAG}.checkpoint.publish.${RUN}.json"; sleep 60; done) &
PUBLISHER=$!
python -m ptb_lnf.cli optimize-free \
  --manifest train.jsonl --ptb "$TASK/ptb" --atompara "$TASK/template.atompara" --basis "$TASK/basis_vDZP" --element Yb \
  --workers "$WORKERS" --bounds 0.9 1.1 --passes 2 --evals-per-param 8 \
  --checkpoint checkpoint.json --exclude-positions excluded_positions.json \
  --ion-targets ion_targets.json --min-sweep-relative-improvement 0.002 --patience 1 > result.json
python -m ptb_lnf.cli apply-free --template template.atompara --out final_candidate.atompara --element Yb --result result.json
python - <<'PY'
from pathlib import Path
def block(path):
    lines=Path(path).read_text().splitlines()
    i=next(i for i,line in enumerate(lines) if line.strip()=='70')
    return lines[i+1:i+16]
a,b=block('template.atompara'),block('final_candidate.atompara')
assert all(a[r]==b[r] for r in (6,11)), 'radial rows changed'
assert a[13].split()[10] == b[13].split()[10], 'radial response changed'
print('SIZE ROWS UNCHANGED',flush=True)
PY
python - <<'PY'
import json
from ptb_lnf.iongate import gate_violation
v,rep=gate_violation('ptb','final_candidate.atompara','basis_vDZP',json.load(open('ion_targets.json')))
json.dump({'violation':v,'pops':rep},open('ions.json','w'),indent=1)
print('FINAL ion violation',v,flush=True)
assert v == 0.0
PY
python -m ptb_lnf.cli evaluate-free --manifest validation.jsonl --ptb "$TASK/ptb" --atompara "$TASK/template.atompara" --basis "$TASK/basis_vDZP" --element Yb --result result.json --workers "$WORKERS" --split validation > validation.json
cp result.json validation.json ions.json final_candidate.atompara excluded_positions.json "$FIX/results/$TAG/"
echo "DONE $TAG"
