# Yb constrained-d7 handover — 2026-09-16

## Completed first constrained screen

- Slurm array `409130` (`yb-d7-refit`) completed cleanly on 17 September.
- Final outputs: `/work/akkleemiss/florian/lnf_orca/yb_d7_constraint_20260916/results/`.

| Fixed `expscal(2,7)` | Final train RMSE | Validation RMSE | Atomic gate |
|---:|---:|---:|---|
| 0.75 | 0.023309414198 | 0.022205569509 | pass |
| 1.00 | 0.023287237686 | 0.022180957841 | pass |
| 1.25 | **0.023245607853** | **0.022130024067** | pass |

The 1.25 branch selects on training, but validation is 0.139633% worse than
the prior unconstrained candidate (0.022099166391). It remains staged only.

## Interior refits and onsite-only follow-up (17-18 September)

Preflight `439456` and refit array `439458` completed; both interior
values pass the atomic gate but do not beat 1.25, so the radial optimum is at
the permitted upper bound.

Job `479348` (`yb-d7-onsite`) then kept d7=1.25 fixed and relaxed only the
three non-radial `shell_xi` onsite-energy coordinates (row 7, cols 9-11),
bounds `[0.92, 1.08]`, starting from the d7=1.25 candidate. Converged after
sweep 2 (relative improvement 4.26e-4 < 5e-4). Scales 0.9302 / 0.9302 /
1.0723 -- interior, not at bounds.

| Branch | Train RMSE | Validation RMSE | Atomic gate |
|---|---:|---:|---|
| d7=1.15 | 0.023264064111 | 0.022152471579 | pass |
| d7=1.20 | 0.023254574933 | 0.022141484752 | pass |
| d7=1.25 | 0.023245607853 | 0.022130024067 | pass |
| d7=1.25 + onsite | **0.023221477108** | **0.022108013726** | pass |

`d7=1.25 + onsite` is the new constrained leader. Its validation is still
0.040% worse than the prior unconstrained candidate (0.022099166391), so it
remains **staged only**. Candidate SHA256
`6793cff547ab0b7c1a8c0f30f1a84b0f16a5d03a5179468302d9101f72aab88b`; local
copy under
`tools/lanthanide_fit/work/weekend_20260911/yb_tm_reseed_20260915/d7_constraint_20260916/cluster_results/d7_1p25_onsite/`.

## Independent Windows reproduction of `6793cff5...` (18 September)

`validate_yb_dshell7_onsite_20260914.py` (env `YB_CANDIDATE`,
`YB_VALIDATION_OUT`) with PTB `ptb_windows_ifx_physical_refine_20260914.exe`
(`17d9b089...`), `ifx_basis_vDZP` (`7f30d98e...`) and the current
NoSpherA2.exe (`0305dc29...`, rebuilt since the 16 Sep run): Yb(0)/(+2)/(+3)
give 24/22/21 electrons at 20 and 80 bohr on both 32- and 64-order grids
(|error| < 1e-8), density-matrix traces exact, RMS radii 1.3401 / 1.2854 /
1.2936 bohr. Output in
`.../d7_constraint_20260916/independent_windows_nosphera2_validation_d7_1p25_onsite/`.

## Joint refit round 1 (job 520292, finished 19 Sep 09:31)

All 13 responsive coordinates refit from the onsite candidate, d7 locked at
1.25, bounds `[0.97, 1.03]`. Ran 11 h on 92 cores; atomic gate passed.

| Branch | Train RMSE | Validation RMSE |
|---|---:|---:|
| prior unconstrained Tm-informed candidate | 0.023297734153 | 0.022099166391 |
| d7=1.25 + onsite | 0.023221477108 | 0.022108013726 |
| d7=1.25 joint round 1 | **0.023134453771** | **0.022016403169** |

**First constrained candidate to beat the unconstrained one on validation**
(-0.374%). SHA256
`bd2ca8a3de184a39018b7efafb0f3ba1e6779395169f285a3f169a63977c2ac4`, local
copy in `cluster_results/d7_1p25_joint/`. Windows IFX/NoSpherA2 reproduction
(20 Sep) passed: 24/22/21 electrons, traces exact, RMS radii
1.3541 / 1.3068 / 1.2918 bohr
(`independent_windows_nosphera2_validation_d7_1p25_joint/`). Still staged
only.

7 of 13 scales ended at the +-3% bound ([6,4], [6,7], [6,8], [7,9], [7,10],
[7,11], [11,8]), so the optimum is not yet interior.

## Joint refit round 2 (manual on AKL007, 20 Sep 10:47 - 22:28)

Slurm controller was down, so `run_yb_d7_1p25_joint2_manual_20260920.sh`
ran under `nohup setsid` on AKL007 with 92 workers. Template = round-1
candidate, bounds re-centred `[0.97, 1.03]`. Three full sweeps, early stop
on convergence (sweep-3 relative improvement 1.2e-4). Temp dir cleaned, no
processes left; AKL007 free well before the 21 Sep 08:00 deadline.

| Branch | Train RMSE | Validation RMSE |
|---|---:|---:|
| d7=1.25 joint round 1 | 0.023134453771 | 0.022016403169 |
| d7=1.25 joint round 2 | **0.023092021505** | **0.021979124020** |

Validation -0.169% vs round 1, -0.543% vs the unconstrained candidate.
Atomic gate passed. SHA256
`9b312272eb5d81614db1fd0b15a52495ac0a03806fa356a5b5e3ad53ff4bb3f9`, local
copy in `cluster_results/d7_1p25_joint2/` (with the run log). Windows
IFX/NoSpherA2 reproduction passed: 24/22/21 electrons, traces exact, RMS
radii 1.3621 / 1.3061 / 1.2875 bohr
(`independent_windows_nosphera2_validation_d7_1p25_joint2/`). Staged only.

9 of 13 scales at the +-3% bound ([6,2], [6,3], [6,4], [6,5], [6,7], [7,9],
[7,10], [7,11], [11,8]); the cumulative drift from the Tm-informed start on
the row-7 onsite columns is now ~1.10x / 1.10x / 0.90x. The optimum is still
not interior, but each round buys less (0.37% -> 0.17% on validation).

## Windows IFX build with the round-2 Yb block (21 Sep)

- Merged atompara: 14 Sep physical-refine release (`71983dec...`) with only
  the Yb (Z=70) block replaced by the round-2 candidate. The Yb fit lineage
  sat on the older La-Lu blocks, but those never enter a Yb-molecule fit, so
  the merge is consistent.
  `tools/lanthanide_fit/work/weekend_20260911/release_yb_d7_joint2_20260921/atompara_ifx_yb_d7_joint2_20260921.atompara`,
  SHA256 `ae2dbc7e0538fe241ea44f6fb310c6f2355897c709e575cf0ea49cbab1600ca0`
  (`manifest.json` alongside).
- Embedded via `tools/lanthanide_fit/work/embed_atompara_defaults.py` into
  `source/default_files.f90` (previous version kept as
  `release_yb_d7_joint2_20260921/default_files.f90.before_20260921`; 18
  lines differ, all Yb).
- Built with ifx/icx 2025.1 + MKL 2025.1 (oneAPI `setvars.bat` is broken on
  this machine - call `compiler5.1\envars.bat` and
  `mkl5.1\envars.bat` after `vcvarsall.bat x64`; 2026.1 has no `icx`
  driver; stale Linux `.o` files in `source/` must be deleted by hand since
  `make clean` does not run under cmd).
  `build/ptb_windows_ifx_yb_d7_joint2_20260921.exe`, SHA256
  `ee173d206f7a1eef4ba81c1d0f1cbefac1fa76c291ddb61c92ddbc5a59814612`.
- Embedded-defaults smoke test (`release_yb_d7_joint2_20260921/embedded_validation.json`):
  no-`-par` run vs `-par` merged file gives bit-identical density matrices
  for Yb/Tm/La at q=0/+2/+3; Yb differs from the 14 Sep exe (max |dP| 1.07),
  Tm and La are identical to it; all traces exact (24/22/21, 23/21/20,
  11/9/8). `build/ptb.exe` was NOT replaced; the dated exe is the artefact.

## Live work

Nothing running. Next round needs either the Slurm controller back
(`run_yb_d7_1p25_joint_20260918.slurm` with template changed to
`results/d7_1p25_joint2/final_candidate.atompara`) or another manual window
on a node.

## Scientific constraint

`expscal(2,7)` is Yb's second-d-shell initial radial scale (shell sequence
`s,s,s,p,p,d,d,f,f`).  The user requires it in `[0.75, 1.25]`; it is locked
absolutely at each tested value.  Do not confuse it with the optional,
default-off `-yb-d7-qexp` charge-response experiment.

The new seed SHA256 values are:

- 0.75: `f87b34134e5637fe4903a5cebb2b04e3ff182f1c7e408e6e60ca3cb282ea78fe`
- 1.00: `b870572418249a455bb4464e5938a17aa326c39c2d95fd283528b436ef96ae56`
- 1.25: `06a582fac1536b4beb0804df9b1e7424618323483b010dacac05213bd37b7047`

All three passed atomic preflight against the original anchor: expected
24/22/21 valence electrons and allowed RMS radii.

## Fit protocol

- Base candidate: Tm-informed Yb candidate SHA256
  `b0b1ea1a22ae795cb8f564a92cb889d958bd2b3291e5422444e5a99ac218b482`.
- Fixed original physical anchor SHA256
  `8f70ffd45c8a2747760e4d012d6f907428f6062213cc1fee2b54e47f1c8fa648`.
- The 13 previously finite-difference-responsive coordinates are relaxed;
  `[6,6]` is explicitly excluded. Three passes, relative bounds `[0.97,1.03]`.
- Original 717 training cases select. The 169 held-out cases only reject
  overfitting. No defaults, source, or executable have been changed.

## After `409130` completes

1. Verify every `result.json`, `final_candidate.atompara`, `physical.json`,
   and `validation.json`; reject incomplete or nonzero-exit branches.
2. Select by final 717-case training RMSE only, then require fixed-anchor
   electron-count/radius gates and 169-case validation.
3. Independently reproduce the selected candidate with Windows IFX PTB and
   NoSpherA2 before staging any release.
4. Interior refinement (1.15, 1.20) is done; retain 1.25 as the endpoint and
   never exceed 1.25.
5. Round 3 from `9b312272...` is justified by the bound hits, but check
   first that the cumulative drift on row 7 / row 6 stays physically
   sensible (RMS radius of Yb(0) has grown 1.323 -> 1.362 bohr over the
   constrained series). Consider widening bounds to `[0.94, 1.06]` to reach
   the optimum in one round instead of two.

## Interpretation

The raw 15-element parameter series is not a direct interpolation target:
14 responsive coordinates formed 13 correlation blocks. The constrained d7
approach is therefore a physically stated model restriction followed by
refitting, not a forced lanthanide trend. The early 1.25 result is better than
the prior unconstrained training value (`0.023297734`) but is not evidence of
acceptance until the final gates pass.

## 21 Sep 2026 � root cause found; d7 constraint thread superseded

The upstream vDZP basis itself is the defect, not the atompara. A Yb atom run
through ORCA 6.1.1 with `%output print[p_basis] 2` (module
`orca/avx2-6.1.1-xqm3jnz`, output `lnf_orca/yb_basis_check_20260921/yb.out`)
prints the same Yb basis as `ptb/.basis_vDZP`, exponent for exponent, including
the second d shell `0.0487 / 0.00207`. So the ORCA references and PTB share the
basis; the fit was never compensating a mismatch. The 0.00207 primitive
(~25x more diffuse than Tm/Lu) is what PTB cannot handle, and no scale inside
[0.75, 1.25] fixes a ~25x problem.

Change (commit `f2e50b9fc`, pushed): Yb second d shell replaced by the Tm/Lu
geometric mean, `0.12687724612 / 0.51052200767E-01`, contraction ratio averaged
and renormalised, in `.basis_vDZP`, `build/ifx_basis_vDZP` (sha
`28fce65c2dc48cd47c57ef4984c0ac15135b47ff90ed39ec1f43fbc909eb964c`) and the
embedded `default_basis`. Startup banner and README carry the deviation notice.
Windows binary: `build/ptb_windows_ifx_yb_basis_d2_20260921.exe` (sha
`e06ce11b0d12b5614e75b07f05cf924e1e1c038d94a2a5866ecefd8a2d09cd7d`). Yb atom
smoke test: old basis gaps 0.000 / 148.1 eV; new basis a normal near-degenerate
open-shell atom, no positive-energy block. NoSpherA2's `Src/basis_data.cpp`
still carries the upstream pair and must be updated to match.

Refit on the new basis: Slurm job **523552** (`yb-basis-d2`),
`lnf_orca/yb_basis_d2_20260921/run_yb_basis_d2_20260921.slurm`. Template =
joint2 candidate with `expscal(2,7)` reset to 1.0 (single-line diff). Step 1
writes a no-fit baseline (`results/basis_d2/baseline_{train,validation,physical}.json`)
so the basis effect is separable from the fit. Step 2: 14 coordinates (d7
free again, `[6,6]` removed from the exclusions), bounds `[0.85, 1.15]`, step
0.05, same gates. d7 therefore stays inside the mandated [0.75, 1.25].
Compare against joint2 `0.023092022 / 0.021979124`; if the baseline alone
already beats it, the d7 series is retired. On success: apply-free ->
embed via `embed_atompara_defaults.py` -> rebuild -> Windows/NoSpherA2
reproduction -> commit, as before. The embedded Yb atompara block is still
joint2 (fitted to the old basis) until then.

## 22 Sep 2026 - job 523552 was void; PTB CLI parsing bug found

**Everything the 21 Sep section says about the refit is superseded.** Jobs
523552, 524614, 524628 and 524646 all fitted the *upstream* basis, not the
replacement, and their results are void. The recorded comparison "basis change
is a wash / buys no additional accuracy" compared the old basis against itself
and must not be relied on.

Cause (fixed in `source/main.f90`): the argument loop matched every CLI flag
with `index(arg1,'-bas') .ne. 0`, a substring test. The slurm staged into
`mktemp -d /tmp/yb-basis-d2.${RUN}.XXXXXX`, whose path contains `-bas` inside
"yb-basis", so when the loop reached the path argument it re-triggered the flag
and set `bname` to the next argument (`-denmat`). PTB could not open that, fell
back to the basis compiled into the binary, and exited 0. All seven flags had
the defect (`-stda`, `-purify`, `-par`, `-bas`, `-chrg`, `-uhf`, `-denmat`);
all now match with `trim(arg1).eq.`.

An earlier diagnosis in this thread blamed an unresolved relative `--basis`
path. That was wrong. The harness guard added for it is retained (it is sound
on its own terms) but it fixed nothing.

Evidence, measured through the production `score_manifest` path on the 717-case
train split: old basis **0.023084336682517047**, new basis **0.027117179860**.
Byte-identical copies of the new basis in differently-named directories:
`/tmp/replicate.aaaaaa/` 0.0280155 (read), `/tmp/yb-basis-d2.524646.KP06F6/`
0.0221671 (not read), `/tmp/ybbasisd2.524646.KP06F6/` 0.0280155 (read),
`/tmp/yb_basis_d2.524646.KP06F6/` 0.0280155 (read). The staged sha was verified
from inside the job immediately before the call that misread it.

Guards now in place:
- `source/main.f90`: exact flag matching. Verified end-to-end on 22 Sep: the
  old binary prints "Basis file '-denmat' not found, using built-in defaults."
  for a `-bas` path containing `-bas`; the patched binary is silent and reads
  the file. **Not yet committed** -
  the cluster binary was deliberately left untouched so the running fit keeps
  one binary throughout.
- `tools/lanthanide_fit/ptb_lnf/optimize.py`: `score_manifest` refuses any
  `--basis`/`--ptb` path containing a PTB flag substring, since every
  already-built binary still has the defect. Check:
  `tools/lanthanide_fit/tests/test_score_manifest_paths.py`.
- slurm `mktemp` pattern changed to `/tmp/yb_basis_d2.${RUN}.XXXXXX`.

Full write-up:
`<vault>/Software-Notes/Debugging-Traps/Fortran-index-Matches-CLI-Flags-Inside-File-Paths-22-Sep-V1.0.md`.

### Corrected refit: job 524666

Same script, de-hyphenated task dir. Baseline confirmed at **0.027117** (new
basis) instead of 0.023084, and the coordinate descent now moves instead of
collapsing back to 1.0. Basis kept at d0511 (the Tm/Lu geometric-mean
interpolant already committed in `f2e50b9fc`); a 169-case validation scan over
outer exponents 0.0511 down to 0.0160 showed a shallow surface with no
dominant candidate, and re-picking on the held-out split would contaminate it.
d0250 is the fallback if the fit lands materially worse than joint2's
0.021979124.

Note the target: the old basis carries a 13-coordinate fit at 0.023092022 /
0.021979124, while the new basis starts unfitted at 0.027117, so the fit has
roughly 19 percent to recover before the replacement breaks even on density
RMSE. The basis change is justified on the electron-count leak and the atomic
pathology regardless (mean dN -1.983 -> +0.283 e over 169 validation cases,
RMS dipole error 55.89 -> 39.97), but if the refit does not close most of that
gap, that trade needs stating explicitly rather than presenting the change as
a pure win.

## 23 Sep 2026 - job 524666 result; continuation 524818

524666 finished 07:17 (253 trials, 3 sweeps, atomic gate passed):

| | new basis, unfitted | new basis, fitted | old basis, joint2 |
|---|---|---|---|
| train | 0.027117 | **0.023707** | 0.023092 |
| validation | 0.026403 | **0.022910** | 0.021979 |

So after 3 sweeps the replacement basis is still 2.7 % (train) / 4.2 %
(validation) behind the old basis on density RMSE. It was not converged: 13 of
14 scales ended at exactly 0.9556^3 = 0.8726 or 1.0444^3 = 1.1393, i.e. every
coordinate took the full `--relative-step` in every sweep, and the sweeps
gained 4.4 / 4.7 / 4.0 % with no decay. The run stopped on `--passes 3`.
Validation fell in step with train, so there is no overfitting signal yet.

Where the gain came from: row 6 cols 7 and 8 (expscal of the two **4f**
shells, 0.8356 -> 0.7291 and 1.4789 -> 1.6850) and row 11 col 8. Not d7.

`expscal(2,7)` is row 6 **col 6** (0-based), not col 7. Fitted value
**1.1393** (template 1.0), inside the mandated [0.75, 1.25]. It rose one full
step per sweep but was worth only ~4e-6 in score per step - a flat coordinate
drifting upward, not the old pathology where it dominated the fit.

Result files: `/work/akkleemiss/florian/lnf_orca/yb_basis_d2_20260921/results/basis_d2/`
(final_candidate.atompara sha256 6b0d5423...b591).

### Continuation: job 524818 (`run_yb_basis_d2c_20260923.slurm`)

Template = 524666's final_candidate, so scales restart at 1.0 with another
[0.85, 1.15] of room. expscal(2,7) = [6,6] added to the exclusions (frozen at
1.1393): it buys nothing and was climbing toward the 1.25 ceiling. 13 free
coordinates, up to 4 sweeps, stops early on a sweep gaining < 0.5 %. The job
first re-scores the start point and aborts unless it reproduces 0.023707 -
the check that would have caught the void jobs of 21-22 Sep.

Open question for Florian: the continuation moves the 4f exponent scales up
to ~30 % from their joint2 values. Nothing mandates a box for them, but it is
a larger move than any earlier local refine.

## 23 Sep 2026 (evening) - d2 continuation converged; gap to joint2 closed on train

- Scoring speed-up: `evaluate_density` in `ptb_lnf/denmat.py` now caches per-atom
  r/r2, contracted radial parts and monomials per chunk and contracts with BLAS.
  Density 14x faster (106 s -> 7.6 s on 12 train structures), RMSE identical to
  4e-15 relative; trial time 168 s -> 28.5 s on 92 cores. Cluster copy:
  `ptb_lnf_fastdens_20260923/`.
- Resume bug fixed in `optimize.py`: a mid-sweep resume measured that sweep's
  improvement from the resume score, not the sweep start, so 526221 early-stopped
  on 0.35 % when the real sweep gain was 2.1 %. Checkpoints now carry
  `sweep_start_score`; a pre-fix checkpoint resumed mid-sweep is not counted for
  early stop. Cluster copy: `ptb_lnf_resumefix_20260923/`.
- d2c/d2d scripts had dropped the `yb_guarded_search.py` copy, so the atomic gate
  crashed; d2e stages it again.
- Chain: 524818 (cancelled at trial 62, 0.023280) -> 526221 (sweep 1 end,
  0.023198) -> 526222 (sweep 2, 0.27 % gain, converged by the 0.5 % criterion).

| candidate | basis Yb d2 | train | validation |
|---|---|---|---|
| d2 start (524666) | 0.1269/0.0511 | 0.023707 | 0.022910 |
| **d2e (526222)** | 0.1269/0.0511 | **0.023136** | **0.022315** |
| joint2 (incumbent) | upstream | 0.023092 | 0.021979 |

- expscal(2,7) = 1.1392905882 (frozen, inside [0.75, 1.25]). Atomic gate passed
  (traces 24/22/21 to 1e-9).
- Of the 0.027117 -> 0.023092 gap the new basis had to recover (see 22 Sep),
  the refit closed 99 % on train (0.2 % behind joint2). Validation stays 1.5 %
  behind (was 4.2 %). Within the new-basis fit, train and validation fell
  together (-2.4 % / -2.6 % from d2), so there's no overfitting signal. The
  residual validation gap is the price of the basis change, which the
  electron-count leak (mean dN -1.983 -> +0.283 e) and dipoles (RMS 55.89 ->
  39.97) justify on their own.
- 4f expscal (row 6 cols 7, 8) is now 0.6780 / 1.7417 vs joint2 0.8356 /
  1.4789, i.e. -19 % / +18 %. That settles the earlier "~30 %" question
  downward, but it is still the largest local move made.
- Not embedded, not committed. Results: `yb_basis_d2_20260921/results/basis_d2e/`
  (final_candidate sha256 874be9a5...c338).
- For Florian: accept d2e as the new Yb default (-> embed, Windows IFX +
  NoSpherA2 reproduction, commit), or refit with a box on the 4f scales first.

## 23 Sep 2026 (late) - d2e release validated; free-ion configuration problem found; d3 running

- Windows build trap: the Makefile links to `build/ptb` (no `.exe`). The first
  `ptb_windows_ifx_yb_basis_d2e_20260923.exe` was a copy of the 7 Sep `ptb.exe`.
  Fixed by copying `build/ptb`; exe sha b6dd5f04...
- `release_yb_basis_d2e_20260923/embedded_validation.json`: embedded == parfile bit
  for bit (Yb/Tm/La, q0/2/3); traces 24/22/21; Tm and La unchanged vs joint2.
- NoSpherA2 rebuilt (its exe predated the basis_data.cpp edit). Independent count
  gate `release_yb_basis_d2e_20260923/independent_windows_nosphera2_validation/`
  passed (24/22/21 to 3e-8 at orders 32/64).
- **Problem:** the free-ion valence RMS radius is 2.57/2.77/2.79 bohr (q0/q2/q3),
  vs joint2 1.36/1.31/1.29 and Er/Tm/Lu neutrals 1.69/1.77/1.81; the cations
  come out larger than the neutral. Free-ion l-blocks decouple, so these
  populations are exact: d2e Yb3+ = d 10.0, f 1.7 (Tm3+: d 0, f 11.7; DFT
  Yb(III) complex Mulliken f 13.03, d 0.42). joint2 parameters on the new basis
  already give d 10, f 2, so the new d shell (not the refit) opens it. The
  molecular fit objective never sees it, and the count gate cannot see it.
- The default_files.f90 embed is NOT committed, pending Florian's decision.
- basis_d3 (job 533835) is running on rows 8-10/12-14; the start gate passed.

## 23-24 Sep 2026 - Yb free-ion configuration fix in progress

Work dir: `tools/lanthanide_fit/work/yb_config_fix_20260923/`. Local molecular
scoring on the downloaded references (`D:/lnf_yb_data`, manifest
`yb_manifest_local.jsonl`) reproduces the cluster train score exactly (d2e
0.023136); ~106 s/trial on 46 workers (`score_local.py`).

- DFT free ions (wB97X-3c, cluster `yb_config_fix_20260923/ions`, job 533901;
  Mulliken l-pops, exact for an atom): Yb0 s4 p6 d0 f14, Yb2+ s2 p6 d0 f14,
  Yb3+ s2 p6 d0 f13. Windows +-0.5 in `yb_ion_targets_draft.json`.
- d2e PTB: q0 s4 p3.3 d8 f8.8, q2 s4 p4.9 d10 f3.1, q3 s4 p5.3 d10 f1.7.
- Seed xi(7d) -5.78 -> -3 (`cand/xi7d_m3.atompara`): train 0.022885 (better than
  d2e), d0 at every charge, but q0 f18 p2 and cations keep 6s2.
- `ptb_lnf/iongate.py` + `optimize-free --ion-targets`: trials with more than
  hard_limit electrons outside the windows score 1 + violation; below it,
  RMSE + weight * violation.
- Atom-only DE (`feas.py`): with exponent rows 6/11 free the feasible points
  cost 0.049-0.055 train. Levels/responses only (rows 7,9,10,12-14):
  `cand/seed_g0_feas_lv41.atompara`, violation 0, train 0.024765.
- Local penalised free fit (`run_g1_local.ps1`, weight 0.02, hard_limit 100,
  [6,6] excluded) reached 0.022444 train, violation 0, at sweep 0 param 73/98
  (`cand/seed_g1_local_ck.atompara`): already better than d2e and physical
  (q0 s4 p5.6 d0 f14.4; q2 s2 p5.7 f14.3; q3 s2 p5.7 f13.3).
- 24 Sep: moved to AKL slurm. Job 542019 (`yb_config_fix_20260923/run_yb_ion_g1.slurm`,
  code `ptb_lnf_iongate_20260923/`) resumes that checkpoint; baseline reproduced
  (0.0247650), resume verified at 0.022444. Output -> `results/ion_g1/`,
  live checkpoint `results/ion_g1.checkpoint.publish.542019.json`.
- d3 (533835) finished: train 0.022367, validation 0.021442, but ion-gate
  violation 57.8 (Yb3+ d10 f1.8). Rejected - the d10 branch is unphysical.
- Still to do after 542019: validation split (overfit check only), atomic count
  gate 24/22/21, Windows IFX + NoSpherA2 reproduction, embed in default_files.f90.
  The uncommitted d2e embed in default_files.f90 must NOT be committed.
- Other lanthanides' free ions are also off (Ho d10, Er q3 d3 f5.9, Dy q3 f4.9,
  Lu f~11-12); the same gate would apply to them.
- ORCA module trap: job shells no longer get the spack MODULEPATH; add
  `module use /work/software/spack/share/spack/lmod/linux-rocky9-x86_64/Core`,
  and ORCA's mpirun needs `--ntasks=4`, not `--cpus-per-task=4`.

## 25 Sep 2026 - ion_g1 finished, continuation ion_g2 queued

- ion_g1 (job 542019, 14.3 h, 1761 trials): train 0.020892, validation
  0.020289 (d2e: 0.023136 / 0.022315), ion violation 0. Free ions: q0 s4 p5.73
  d0 f14.27; q2 s2 p5.71 d0 f14.29; q3 s2 p5.56 d0 f13.44 (totals 24/22/21).
  Validation improves as much as train, so no overfit.
  Local copy: `work/yb_config_fix_20260923/ion_g1/`.
- Windows IFX reproduction (`score_local.py`, `gate_check.py`): identical
  train/validation/ion numbers.
- Ended on the 3-pass limit with sweep 3 still gaining 1.8 %, so ion_g2 (job
  572261, `run_yb_ion_g2.slurm`) continues from ion_g1's final_candidate: 4
  passes, stop below 0.2 % per sweep, same gate and [6,6] exclusion.
- Still to do: take ion_g2 on its train score, falling back to ion_g1 only if g2's
  validation shows overfitting (validation is never used for selection); embed in
  default_files.f90, rebuild, NoSpherA2 run on Windows, commit.

## 28 Sep 2026 - ion_g2 validated and embedded

- Job 572261 completed in 26:02:35. Its four sweeps improved 2.306%, 2.433%,
  0.966%, and 1.495%; it reached the pass cap, not the 0.2% stop criterion.
- Selected ion_g2 on training score: 0.019426166410 versus ion_g1 0.020891689640.
  Held-out validation also improved: 0.018851480493 versus 0.020289.
  Windows IFX full train/validation scoring reproduced 0.019426/0.018851.
  Cluster candidate SHA256: 6eb12794d40838ce9b96050df769be7daa1c6e311d941b5dd6b048864ea51c39;
  the downloaded file has the same hash.
- Free ions: q0 s4 p5.932 d0 f14.068; q2 s2 p5.708 d0 f14.292;
  q3 s2 p5.627 d0 f13.373. Ion penalty 0; electron totals 24/22/21.
- `source/default_files.f90` now embeds the candidate's Yb block only.
  Windows IFX build (`build/ptb_windows_ifx_yb_ion_g2_20260928.exe`, SHA256
  308d82533af842efd7dcd2c3f3e8cf65f175816d6de6c10f5127a50e199de22f)
  gave bit-identical embedded/parameter-file density matrices for Yb q0/2/3,
  and unchanged Tm/La q0/2/3 versus the earlier joint2 executable.
- Independent NoSpherA2 density export passed for q0/2/3 at quadrature orders
  32 and 64, including trace/count convergence. RMS radii are about
  1.87/1.93/1.93 bohr, improved from d2e's 2.57/2.77/2.79.
  Results: `tools/lanthanide_fit/work/yb_config_fix_20260923/ion_g2/`.

## 28 Sep 2026 - size-preserving continuation in progress

- Direct AO extent check on the released basis/parameter file (same Windows
  IFX executable, atomic q0/2/3, Gaussian quadrature orders 32/64): second d
  shell RMS radius Tm/Yb/Lu = 3.607/4.394/5.301 bohr; f shells 8 and 9 are
  Tm/Yb/Lu = 0.710/0.698/0.686 and 2.037/1.983/1.943 bohr. The Yb shell
  extents are bracketed by its neighbors, and each AO shell's density outside
  20 bohr is negligible. `shell_sizes_20260928/result.json` gives all shells
  and charges. These are AO extents, not occupied l-population radii.
- The free-atom Yb AO overlap has p/f cross terms (maximum 0.482), so the
  gate's P.S p/f numbers are a basis partition, not exact angular-momentum
  occupations. Keep the DFT-matched population windows as a diagnostic gate;
  do not interpret f=14.07 as an exact f electron count.
- Independent NoSpherA2 atom radii for Tm q0/2/3: 1.771/1.542/1.828 bohr;
  Lu: 1.807/1.741/1.735; Yb ion_g2: 1.866/1.933/1.935. Yb cations remain
  larger, so a fit must not worsen their sizes. Neighbor control result:
  `neighbor_exports_20260928/result.json`.
- Size-preserving continuation `ion_g3_size`, job 596350, started from the
  ion_g2 candidate. It freezes all row 6 and row 11 radial scales and
  (13,10), the shell response term used by the self-consistent radial update;
  retains the ion population gate, 717 training/169 validation split, and
  exact basis hash. Two sweeps with [0.9,1.1] per-value bounds. Job 596349
  was cancelled after 2:57 because the first script omitted (13,10); its
  output must not be used. Job 596350 reproduced 0.019426166410 before fitting.
- `ion_g3_size_experiment.json` contains the acceptance criteria. The Windows
  background watcher `watch_ion_g3.py` polls the job,
  copies completed results, reproduces scores, checks ions/AO extents, runs
  independent NoSpherA2 exports, and writes `ion_g3_size/verdict.json`.
  A passing verdict is staged for review; no default update is automatic.

## 29 Sep 2026 - g3 rejected on neutral size; size-passing candidate and g4

- Job 596350 completed normally in 10:45:21. g3 train 0.019140684992,
  validation 0.018611089035, ion violation 0, AO shell radii unchanged,
  independent NoSpherA2 counts passed. Its neutral Yb RMS radius rose
  1.865600 -> 1.909010 bohr (+2.3%), beyond the preregistered +1% limit;
  `ion_g3_size/verdict.json` rejects it. Committed ion_g2 defaults remain.
- Isolated Windows IFX build of g3 (`build/ion_g3_rc_20260929/build/ptb`, SHA256
  fb7756d19669b677fc2baa5eb23c5ee01013d0da2b3d97de4cfd7c7aaaaf2027)
  exactly reproduced g3's Yb parameter-file atom densities and unchanged
  Tm/La controls. This is a diagnostic binary, not the default executable.
- Straight interpolation g2->g3 at fractions 0.25/0.4/0.5/0.75 was rejected:
  Yb2+ switched to s4/p~5.7/f~12.3, ion violation ~2.6. Single-coordinate
  backoffs from g3: 79 tested, 71 with ion violation 0, none with acceptable
  neutral radius. No molecular score was used for these rejected states.
- Coordinated row backoffs: 127 combinations of rows 7/8/9/10/12/13/14 were
  ion/radius-screened; 24 passed both. Five minimal feasible groups were scored
  on the same 717-case training set. Best was case_043, reverting rows 7,13,14
  to g2 while retaining g3's rows 8,9,10,12: train 0.019309259096,
  validation 0.018774555259, versus g2 0.019426166410/0.018851480493.
  Candidate SHA256 cac3f24eba7fa4238cfc3c9cd5809e2548203e4307f2facea8e6f5fde53011c1.
- case_043 passed Windows free-ion gate (24/22/21) and independent NoSpherA2
  counts at orders 32/64. q0/q2/q3 radii: 1.870714/1.929403/1.934039 bohr;
  all within +1% of ion_g2. The d/f AO radii are exactly unchanged.
  Isolated embedded IFX candidate `build/ion_g3_sizepass_rc_20260929/build/ptb`
  SHA256 d29020ab8d6507d005fdc2c4bc092059a451c9fd2427c329b82ea8495388559c
  reproduced train/validation and bit-identical embedded/parameter-file Yb
  atoms; Tm/La remained bit-identical to ion_g2.
- Further bounded fit job 603385 (`ion_g4_sizepass`) was submitted from case_043
  on 29 Sep after checking for queued/running jobs. It freezes rows 6,7,11,13,14;
  only rows 8,9,10,12 are free. The 717/169 split and ion gate remain. At
  submission it was PENDING (Priority). `watch_ion_g4.py` polls and independently
  verifies its final candidate. No default replacement is automatic.

## 1 Oct 2026 - ion_g4_sizepass completed and independently verified

- Job 603385 completed normally in 5:59:59 (two sweeps, 705 trials). It started from case_043 and kept Yb rows 6,7,11,13,14 fixed. Final SHA256: `aef910c4fb85d7d668cba2c2574661f7c68cbd8cad8df8c3763c3bca3643ddf8`.
- Train 0.019127715299 (case_043 0.019309259096); held-out validation 0.018606902088 (case_043 0.018774555259). Windows IFX parfile reproduction: 0.019127708791 / 0.018606901649.
- Ion violation zero; free-ion totals 24/22/21. Fixed radial rows unchanged; maximum AO shell radius delta versus ion_g2 zero. Independent NoSpherA2 export at radial orders 32/64 passed counts and convergence for Yb 0/2+/3+. Atomic RMS radii 1.864584/1.929254/1.934823 bohr versus ion_g2 1.865600/1.933280/1.934601.
- Full gate verdict: `tools/lanthanide_fit/work/yb_config_fix_20260923/ion_g4_sizepass/verdict.json` says `passed_for_review`. The first watcher stopped on transient SSH timeout; a direct rerun completed every local check. Candidate is staged only; default_files.f90 has not been changed.

## 1 Oct 2026 - g4 embedded and release binaries validated

- Promoted g4 SHA256 `aef910c4fb85d7d668cba2c2574661f7c68cbd8cad8df8c3763c3bca3643ddf8` into `source/default_files.f90`; only Yb rows 8,9,10,12 changed from released ion_g2. The frozen d/f radial and branch-control rows remain identical.
- Clean Windows IFX build SHA256 `06c7e871eb5da36ee9fedbd34695c376be0312c26e96df7e4a89b603fa2bc13a`. Embedded Yb q0/2/3 density matrices match the parameter file bitwise; Tm/La match the prior release bitwise. Full Windows train/validation reproduction: 0.019127708791 / 0.018606901649. Independent NoSpherA2 q0/2/3 export at orders 32/64 passed.
- Release files under `build/release_g4_20261001/`: `ptb.exe` (Windows x64), `ptb_linux_static` (Linux x86-64, fully static), `ptb_macos` (universal arm64/x86-64, Accelerate and libSystem only), plus SHA256SUMS and release notes. All architectures passed Yb q0/2/3 default runs; electron counts 24/22/21 and density matrices within 2.19e-6 of Windows. `RELEASE_G4_20261001.md` and `ion_g4_sizepass/release_manifest.json` preserve hashes and gates.