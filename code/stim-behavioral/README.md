# STIM Behavioral Analysis Pipeline

Behavioral + surprise-memory analysis for the flanker ("stim_v1") task.
Four scripts, run in this exact order (each depends on the previous one's
output):

```
1. stim-behavioral_extract-trials.py    (Python)  raw CSVs -> trial data + summary
2. stim-behavioral_select-participants.py       (Python)  summary  -> valid-subject list
3. stim-behavioral_prepare-memory-data.R      (R)       raw CSVs -> organized flanker/surprise CSVs   [valid subjects only]
4. stim-behavioral_compute-memory-measures.R   (R)       organized CSVs -> memory measures             [valid subjects only]
```

`ndclab_py.py` is a helper module, never run directly; it must sit in the
same folder as the two Python scripts (i.e. in `code/`).

These scripts are adapted by Ana García Morazzani from other lab scripts (NDCLab @ FIU)
written by Kianoosh Hosseini and Felix Zakirov.

The R scripts only process subjects that the Python pipeline marked as fully valid: they read
`stim_valid_behav.csv` and stop with an error if it does not exist yet.

---

## Repository layout

All code and outputs live inside the dataset (BIDS-style), with each
pipeline getting its own named `derivatives/` subfolder:

```
datasets/stim-dataset/
  sourcedata/raw/s1_r1/psychopy/sub-*/      inputs (raw PsychoPy CSVs, per-subject subfolders)
  code/                                      these scripts + ndclab_py.py + this README
  derivatives/
    stim-behavioral/s1_r1/                   steps 1-2 output (trial data, summary, valid list)
    stim-memory/s1_r1/
      csv_output/                            step 3 output (organized flanker/surprise CSVs)
      stat_output/                           step 4 output (memory measures)
```

`{dataset}` below is `datasets/stim-dataset`.

---

## Prerequisites

- **Python** (3.x) with `pandas` and `numpy` (the PyCharm venv works).
- **R** with `tidyverse` (includes `dplyr`, `stringr`) and `psycho`.
- The two R output folders must exist before running the R scripts (only
  the Python scripts create their own output dirs):
  `derivatives/stim-memory/s1_r1/csv_output/` and
  `derivatives/stim-memory/s1_r1/stat_output/`.
- The raw PsychoPy CSVs live in per-subject subfolders
  (`.../psychopy/sub-410002/sub-410002_stim_v1_...csv`). All four scripts
  read this layout (the R scripts use `list.files(..., recursive = TRUE)`),
  so no flat copy is needed.

All paths are hardcoded at the top of each script (Section 2/3, marked as
the only part you should need to touch). Session is hardcoded to `s1_r1`.

---

## Step 1 -- stim-behavioral_extract-trials.py

**What it does:** finds every `sub-*` folder, extracts the 384 main flanker
trials per subject (12 shuffled blocks x 32 trials), recodes all variables,
links each trial to its temporal neighbors (5 s adjacency rule) for the
post-error measures, and computes one summary row per subject.

**Input:**
- `{dataset}/sourcedata/raw/s1_r1/psychopy/sub-*/sub-*_stim_v1_psychopy_s1_r1_e1.csv`

**Output** (all in `{dataset}/derivatives/stim-behavioral/s1_r1/`):
- `sub-{id}_trial_data.csv` -- one per subject; trial-level data with
  `pre_*`/`next_*` neighbor columns
- `summary_s1_r1_{timestamp}.csv` -- one row per subject FOUND (skipped
  subjects included with `success = 0`)
- `full_df_{timestamp}.csv` -- all trial_data files stacked (built from
  every `sub-*_trial_data.csv` in the folder -- delete stale ones if a
  subject is removed from the raw data)
- `{timestamp}_log.txt` -- console output copy (skip reasons live here)

**Skip rules:** no matching CSV / multiple CSVs / fewer than 346 (90%)
main-task trials.

## Step 2 -- stim-behavioral_select-participants.py

**What it does:** applies the usability criteria to the newest summary and
prints a participant breakdown (total / valid / removed, with IDs per
reason).

**Criteria, in order:** processed in step 1 -> accuracy >= 0.6 ->
at least 6 valid commission errors -> not in the manual `exclude_id_list` ->
within +/-3 SD on `invalid_rt_percent` and `skipped_percent`.
(Note: the `skipped_percent` outlier criterion is one the lab manual says to
confirm with the PI.)

**Input:**
- newest `summary_s1_r1_*.csv` from step 1 (picked automatically by
  modification time)

**Output:**
- `stim_valid_behav.csv` (same folder) -- the inclusion list; row set =
  fully valid subjects. **This file gates steps 3 and 4.**

## Step 3 -- stim-behavioral_prepare-memory-data.R

**What it does:** for each PYTHON-VALID subject, re-extracts the flanker
trials and the surprise memory-test trials from the raw CSV into two tidy
per-subject files (conventions: accuracy 0 = error or no response;
legit response = responded and RT >= 150 ms; memory trial kept if
RT > 200 ms).

**Input:**
- raw CSVs in `{dataset}/sourcedata/raw/s1_r1/psychopy/sub-*/` (read
  recursively)
- `stim_valid_behav.csv` from step 2 (gating)

**Output** (in `{dataset}/derivatives/stim-memory/s1_r1/csv_output/`):
- `{id}_stim_flankerDat.csv` -- one row per flanker trial
- `{id}_stim_surpriseDat.csv` -- one row per surprise memory trial

## Step 4 -- stim-behavioral_compute-memory-measures.R

**What it does:** for each Python-valid subject, applies the memory-study
inclusion criteria, computes flanker performance measures and the memory
hit-rate contrast (faces from incongruent-error vs incongruent-correct
trials), assigns the between-subject condition from hardcoded ID lists,
and prints group summary stats.

**Additional inclusion criteria** (on top of the Python gating; counters
for each are kept in the R environment but not printed):
1. <= 20% of surprise trials removed for RT <= 200 ms (attentiveness)
2. flanker accuracy >= 0.6 (fast trials excluded, no-responses count as 0)
3. >= 8 legit incongruent errors
4. >= 8 of those error faces present in the RT-filtered surprise data

**Input:**
- raw CSVs (for the participant id only)
- `{id}_stim_flankerDat.csv` / `{id}_stim_surpriseDat.csv` from step 3
- `stim_valid_behav.csv` from step 2 (gating)

**Output** (in `{dataset}/derivatives/stim-memory/s1_r1/stat_output/`):
- `processed_data_stim_Proj_v1.csv` -- one row per included participant:
  accuracies, RT means (raw + log(1+rt)), flanker effects, error_hitRate,
  correct_hitRate, hitRate_error_minus_correct, condition.
  NOTE: overwritten on every run (no timestamp in the filename).

**Maintenance:** new participants must be added to the hardcoded condition
lists (Section 5) or they get condition NA in the group summary.

---

## Conventions (shared across the pipeline, per the lab manual)

- accuracy: 1 = correct, 0 = error; missing responses count as 0 for
  accuracy but are excluded from RT/post-error analyses
- congruency: 1 = congruent, 0 = incongruent
- RT < 150 ms = invalid; excluded from accuracy, RT, and post-error
  analyses (RT = exactly 150 ms is valid)
- keys: "1" = left, "8" = right; only the first response counts toward
  accuracy/RT; trials with multiple responses are excluded from post-error
  analyses (as the N trial; allowed as the N+1 trial)
- rt_corr / rt_err in the Python summary are INCONGRUENT-only by design
  (congruency held fixed when studying accuracy effects), so
  rt_corr == rt_incon
- error/correct analyses require >= 6 errors (Python); the memory study
  tightens this to >= 8 legit incongruent errors (R, study-specific)

## Known differences that remain between the Python and R stages

- log-RT: Python uses mean(log(rt)) x 1000; R uses mean(log(1 + rt)).
  Not comparable; changing R to log(rt) is pending a decision with the
  original author.
- The R stage applies stricter inclusion (criteria 1-4 above), so its
  final N is a subset of the Python-valid N. This is intentional.
.
