# ==============================================================================
# stim-behavioral_extract-trials.py  --  STAGE 1 of the STIM Behavioral Analysis Pipeline
# ==============================================================================
#   This script reads each subject's raw PsychoPy CSV, extracts and recodes the
#   12 blocks x 32 trials = 384 main-task trials, links each trial to its temporal
#   neighbors (for post-error measures), and writes, into {output_dataset_path}{output_path}:
#     - sub-{id}_trial_data.csv        one per subject, trial-level data
#     - summary_{session}_{time}.csv   one row per subject, incl. skipped subs
#     - full_df_{time}.csv             all trial_data files stacked together
#     - {time}_log.txt                 copy of everything printed to console
#
#   ADAPTED FROM: thrive-theta-ddm behavior_analysis.py.
#   All deviations from the original are marked with "# Edit by Ana:" (except for section banners
#   and reader's notes.)
#
#   KEY CONVENTIONS:
#   accuracy:        1 = correct, 0 = error
#   congruent:       1 = congruent, 0 = incongruent
#   *_R suffix:      0 = left, 1 = right  (keys: "1" = left, "8" = right)
#   rt:              seconds; NaN = no response
#   valid_rt:        1 = RT >= 150 ms  (note: no-response trials count as 1!)
#   block_num:       TRUE block identity (1-12) from the stimulus file name
#   block_order:     position in which the block was presented (1-12);
#                    differs from block_num because blocks are shuffled
#                    within each group of 3
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: IMPORTS & PANDAS SETTINGS
# ------------------------------------------------------------------------------
import sys
import os
import glob
import pandas as pd
import numpy as np
import re
import time
import datetime

pd.options.mode.chained_assignment = None

# ------------------------------------------------------------------------------
# SECTION 2: HELPER FUNCTIONS
# ------------------------------------------------------------------------------
# - PsychoPy stores RTs and key presses as STRINGS that look like lists,
#   e.g. "[0.6625208999030292]" or "['1', '8']". These two parsers convert
#   them into usable values.
# - convert_to_list_rt: returns ONE float per trial. If the participant
#   pressed multiple keys, only the FIRST RT is kept. NaN = no response.
# - convert_to_list_resp: returns the full list of key codes per trial
#   (e.g. [1, 8]) so that multiple responses can be detected later
#   (extra_resp column). NaN = no response.
# - Tee: lets print() output go to BOTH the console and the log file
#   simultaneously
# ------------------------------------------------------------------------------
def convert_to_list_rt(series):
    float_list = []
    for value in series:
        if isinstance(value, str):
            if "," in value.strip("[]"):
                float_list.append([float(v) for v in value.strip("[]").split(",")][0]) # Check if the value is a string
            else:
                float_list.append(float(value.strip("[]"))) # Convert string to float and remove brackets
        elif isinstance(value, list): # Check if the value is a list
            float_list.extend([float(v) for v in value]) # Convert each element of the list to float
        else: # Handle NaN values
            float_list.append(np.nan) # Append NaN if value is NaN
    return float_list

def convert_to_list_resp(series):
    resp_list = []
    for value in series:
        if isinstance(value, str):
            converted_row = list(map(int, re.findall(r'\d+', value)))
            resp_list.append(converted_row)
        else: # Handle NaN values
            resp_list.append(np.nan) # Append NaN if value is NaN
    return resp_list

# Edit by Ana: Tee stdout to both the console and a log file (original redirected
# stdout entirely into the log, which hides all output when running locally.)
class Tee:
    def __init__(self, *streams):
        self.streams = streams
    def write(self, msg):
        for s in self.streams:
            s.write(msg)
    def flush(self):
        for s in self.streams:
            s.flush()

# ------------------------------------------------------------------------------
# SECTION 3: CONFIGURATION -- paths, session, task constants
# ------------------------------------------------------------------------------
# - THIS IS THE ONLY SECTION YOU SHOULD NEED TO TOUCH for a normal re-run.
# - input:  {input_dataset_path}sourcedata/raw/{session}/psychopy/sub-*/
# - output: {output_dataset_path}derivatives/behavior/{session}/
#           (created automatically if it does not exist)
# - valid_trial_count: a subject is skipped entirely if fewer than 90% of
#   the expected 384 main-task trials are present (= 346).
# - adjacency_thresh (5 s): two consecutive trials are treated as neighbors
#   (for pre_*/next_* linking) only if their response-window offsets are
#   within 5 s. Within-block gaps here are ~3.6-4.0 s (jittered ISI) and
#   block breaks are >= 22 s, so 5 s cleanly separates "same block" from
#   "across a break". DO NOT lower this below ~4 s or all trials will be
#   treated as non-adjacent and PES/PEA will silently become NaN.
# ------------------------------------------------------------------------------
start = time.time()
session = "s1_r1"  # Edit by Ana: hardcoded; was sys.argv[1]

# Edit by Ana: local macOS paths; layout is raw/{session}/psychopy/sub-*/
input_dataset_path = "/Users/anagarci/Documents/NDClab/datasets/stim-dataset/"
output_dataset_path = "/Users/anagarci/Documents/NDClab/datasets/stim-dataset/"  # Edit by Ana: derivatives now live inside the dataset (BIDS-style), not a separate analyses/ tree
data_path = f"sourcedata/raw/{session}/psychopy/"
output_path = f"derivatives/stim-behavioral/{session}/"  # Edit by Ana: pipeline-named folder (was derivatives/behavioral/)

os.makedirs(f"{output_dataset_path}{output_path}", exist_ok=True)  # Edit by Ana: create output dir locally

# ------------------------------------------------------------------------------
# SECTION 4: LOG FILE SETUP
# ------------------------------------------------------------------------------
# - Everything printed from here on is written to the console AND to
#   {time}_log.txt in the output folder. The timestamp in the log filename
#   matches the timestamp in the summary/full_df filenames of the same run.
# ------------------------------------------------------------------------------
date_time = datetime.datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
log_file = open(f"{output_dataset_path}{output_path}{date_time}_log.txt", "wt")
sys.stdout = Tee(sys.__stdout__, log_file)  # Edit by Ana: was sys.stdout = open(...)

n_blocks = 12   # Edit by Ana: was 20
n_trials = 32   # Edit by Ana: was 40
valid_rt_thresh = 0.150
valid_trial_count = int(np.ceil(0.9 * n_blocks * n_trials))  # Edit by Ana: 90% of 384 = 346 (was 360/400 per condition)
adjacency_thresh = 5  # Edit by Ana: was 3; ISI here is jittered 3.5-4.0 s so consecutive trials
                      # are ~3.6-4.0 s apart, while block breaks are >=22 s. 5 s separates cleanly.

# ------------------------------------------------------------------------------
# SECTION 5: SUBJECT DISCOVERY & SUMMARY SCAFFOLDING
# ------------------------------------------------------------------------------
# - Subjects = every sub-* folder found under the input path; the numeric ID
#   is extracted from the folder name.
# - processing_log is a dict of lists, one entry appended per subject, that
#   becomes the summary CSV at the end. EVERY subject gets a row, including
#   skipped ones (success = 0, measures = NaN) -- stim-behavioral_select-participants.py
#   relies on this to count total participants.
# - log_skip(): appends the skip row consistently so all lists stay the same
#   length (the original crashed/misaligned when a subject was skipped).
# ------------------------------------------------------------------------------
# Edit by Ana: subject discovery from flat sub-* folders under raw/{session}/psychopy/
sub_folders = glob.glob(f"{input_dataset_path}{data_path}sub-*")
subjects = sorted(set([re.findall(r'\d+', item.split(os.sep)[-1])[0] for item in sub_folders]))
print(subjects)
processing_log = dict()
# Edit by Ana: single condition -> no _soc/_nonsoc prefixes anywhere
summary_columns = [
            "n_trials_valid", "invalid_rt_percent", "skipped_percent",
            "acc", "acc_con", "acc_incon", "rt_con", "rt_incon", "rt_corr", "rt_err",
            "rt_con_log", "rt_incon_log", "rt_corr_log", "rt_err_log",
            "pes", "pea", "peri_acc", "peri_rt", "6_or_more_err",
            ]
processing_log["sub"] = []
processing_log["success"] = []
processing_log["n_trials"] = []
for colname in summary_columns:
    processing_log[colname] = []

# Edit by Ana: helper to keep the processing_log dict aligned when a subject is skipped
# (fixes original bug where a `continue` without appending produced unequal list lengths)
def log_skip(reason):
    print(reason)
    processing_log["success"].append(0)
    processing_log["n_trials"].append(0)
    for c in summary_columns:
        processing_log[c].append(np.nan)

# ==============================================================================
# SECTION 6: PER-SUBJECT PROCESSING LOOP
# ==============================================================================
for sub in subjects:
    processing_log["sub"].append(sub)
    subject_folder = f"{input_dataset_path}{data_path}sub-{sub}{os.sep}"

    # --------------------------------------------------------------------------
    # SECTION 6a: FILE DISCOVERY / SKIP RULES
    # --------------------------------------------------------------------------
    # - A subject is skipped (success = 0 in the summary) if:
    #     * no CSV matches the expected filename, or
    #     * more than one CSV matches (protocol deviation -> handle manually)
    # - Skips are printed to console + log with the reason.
    # --------------------------------------------------------------------------
    # Edit by Ana: simplified deviation handling -- expect exactly one matching CSV, otherwise skip + log
    pattern = f"{subject_folder}sub-{sub}_stim_v1_psychopy_{session}_e1.csv"  # Edit by Ana: filename pattern
    filename = glob.glob(pattern)
    if len(filename) == 0:
        log_skip(f"sub-{sub} has no matching psychopy CSV, skipping...")
        continue
    elif len(filename) > 1:
        log_skip(f"sub-{sub} has {len(filename)} matching CSVs (deviation?), skipping...")
        continue

    # --------------------------------------------------------------------------
    # SECTION 6b: LOAD RAW CSV & TRIM TO MAIN-TASK TRIALS
    # --------------------------------------------------------------------------
    # - The raw CSV also contains practice trials, resting segments, and the
    #   surprise memory phase. Keeping only rows AFTER the first
    #   "task_blockText.started" value and dropping rows without a middleStim
    #   isolates exactly the 384 main flanker trials (verified on pilot data).
    # - Third skip rule lives here: fewer than 346 main-task trials
    #   (e.g. task aborted early) -> subject skipped.
    # - A non-fatal WARNING is printed if the count is >= 346 but != 384.
    # --------------------------------------------------------------------------
    print("Processing sub-{}...".format(sub))
    data = pd.read_csv(filename[0])
    start_index = data["task_blockText.started"].first_valid_index()
    data = data.iloc[start_index:, :].dropna(subset = "middleStim")
    data = data.reset_index(drop=True)  # Edit by Ana: no conditionText filter (single condition)
    # Edit by Ana: was a hard assert on exactly n_blocks*n_trials; now skip if below 90% threshold
    if len(data) < valid_trial_count:
        # Edit by Ana: undo the sub/success appends pattern by logging skip consistently
        log_skip(f"sub-{sub} has only {len(data)} main-task trials (<{valid_trial_count}), skipping...")
        continue
    if len(data) != n_blocks * n_trials:
        print(f"WARNING: sub-{sub} has {len(data)} trials (expected {n_blocks * n_trials}), processing anyway...")
    processing_log["success"].append(1)
    processing_log["n_trials"].append(len(data))

    # --------------------------------------------------------------------------
    # SECTION 6c: TRUE BLOCK IDENTITY (blocks are shuffled!)
    # --------------------------------------------------------------------------
    # - Blocks run in 4 groups of 3 (1-3, 4-6, 7-9, 10-12) and are SHUFFLED
    #   within each group. Each trial row carries the name of its source
    #   stimulus file (e.g. "flanker_block7.csv") in one of the four
    #   whichBlock* columns; bfill(axis=1) merges them into one column.
    # - Later, block_num = the number extracted from that filename (true
    #   identity) and block_order = presentation position. Use block_num to
    #   match stimulus lists; use block_order for time-on-task effects.
    # --------------------------------------------------------------------------
    # Edit by Ana: coalesce the four whichBlock* columns (each filled only for its group of 3 blocks)
    # into one column giving the source block CSV for every trial
    data["whichBlock"] = data[["whichBlock123", "whichBlock456", "whichBlock789", "whichBlock101112"]].bfill(axis=1).iloc[:, 0]

    # --------------------------------------------------------------------------
    # SECTION 6d: BUILD TRIAL-LEVEL DATAFRAME & RECODE VARIABLES
    # --------------------------------------------------------------------------
    # - Keeps only the columns needed downstream, parses the PsychoPy
    #   list-strings (Section 2 helpers), then derives binary codes:
    #     target_R, fl_direction_R (flanker direction, inferred from
    #     target x congruency), valid_rt, no_resp, trial/block indices.
    # - main_stim_image (the trial's background face) is carried through as
    #   a passenger column for the later memory analysis; it is not used here.
    # - CAUTION if reusing valid_rt: no-response trials have rt = NaN, and
    #   NaN < 0.150 is False, so they get valid_rt = 1. They are excluded
    #   from summaries via no_resp / accuracy filters instead (original
    #   behavior, kept for comparability).
    # --------------------------------------------------------------------------
    trial_data = data[[
        "target",
        "congruency",                   # Edit by Ana: was "congruent"
        "accuracy",
        "task1_stim_keyResp.rt",        # Edit by Ana: was task_stim_keyResp.*
        "task1_stim_keyResp.stopped",
        "task1_stim_keyResp.keys",
        "whichBlock",                   # Edit by Ana: was "conditionText"
        "main_stim_image",              # Edit by Ana: carry the trial's face image for the memory analysis later
    ]]
    trial_data["rt"] = convert_to_list_rt(trial_data["task1_stim_keyResp.rt"])
    trial_data.drop("task1_stim_keyResp.rt", axis = 1, inplace = True)
    assert (np.sum([type(i) != float for i in trial_data["rt"]]) == 0), "Check your RT!"

    trial_data["resp_direction_R"] = convert_to_list_resp(trial_data["task1_stim_keyResp.keys"])
    trial_data.drop("task1_stim_keyResp.keys", axis = 1, inplace = True)

    trial_data.columns = [
        "target",
        "congruent",
        "accuracy",
        "task1_stim_keyResp.stopped",
        "whichBlock",
        "main_stim_image",
        "rt",
        "resp_direction_R",
    ]
    # Edit by Ana: true block identity from the (shuffled) block CSV filename, plus sequential order.
    # Original assigned block_num positionally, which would mislabel shuffled blocks.
    trial_data["block_num"] = trial_data["whichBlock"].str.extract(r'block(\d+)').astype(int)
    trial_data["block_order"] = (trial_data["block_num"] != trial_data["block_num"].shift()).cumsum()
    trial_data.drop("whichBlock", axis = 1, inplace = True)

    trial_data["target_R"] = [0 if i == "left" else 1 for i in trial_data["target"]]
    trial_data.drop("target", axis = 1, inplace = True)

    trial_data["fl_direction_R"] = [
                                    0 if
                                    (
                                        (trial_data.loc[i, 'target_R'] == 0 and trial_data.loc[i, 'congruent'] == 1) or
                                        (trial_data.loc[i, 'target_R'] == 1 and trial_data.loc[i, 'congruent'] == 0)
                                    )
                                    else 1 if
                                    (
                                        (trial_data.loc[i, 'target_R'] == 0 and trial_data.loc[i, 'congruent'] == 0) or
                                        (trial_data.loc[i, 'target_R'] == 1 and trial_data.loc[i, 'congruent'] == 1)
                                    )
                                    else None
                                    for i in range(len(trial_data))
                                ]
    trial_data["valid_rt"] = [0 if i < valid_rt_thresh else 1 for i in trial_data["rt"]]
    trial_data["no_resp"] = [1 if np.isnan(i) else 0 for i in trial_data["rt"]]
    trial_data["trial_num"] = [i for i in range(1, len(trial_data)+1)]
    trial_data["first_trial"] = [1 if i == 0 else 0 for i in range(len(trial_data))]
    trial_data["last_trial"] = [1 if i == (len(trial_data)-1) else 0 for i in range(len(trial_data))]

    # --------------------------------------------------------------------------
    # SECTION 6e: RESPONSE DIRECTION & MULTIPLE-RESPONSE DETECTION
    # --------------------------------------------------------------------------
    # - Key "1" = left (coded 0), key "8" = right (coded 1). Only the FIRST
    #   key press determines resp_direction_R.
    # - extra_resp = 1 if the participant pressed more than one key in the
    #   response window; used to exclude messy n-1 trials from PES/PEA.
    # --------------------------------------------------------------------------
    extra_resp = []
    resp_direction = []
    for i in range(len(trial_data)):
        row = trial_data.loc[i, "resp_direction_R"]
        if type(row) == list:
            if row[0] == 1:
                resp_direction.append(0)
            elif row[0] == 8:
                resp_direction.append(1)
            if len(row) > 1:
                extra_resp.append(1)
            else:
                extra_resp.append(0)
        elif np.isnan(row):
            extra_resp.append(np.nan)
            resp_direction.append(np.nan)

    trial_data["resp_direction_R"] = resp_direction
    trial_data["extra_resp"] = extra_resp

    # --------------------------------------------------------------------------
    # SECTION 6f: PRE/NEXT TRIAL LINKING (adjacency)
    # --------------------------------------------------------------------------
    # - Every trial variable is duplicated as pre_* (previous trial's value)
    #   and next_* (next trial's value). These feed PES/PEA/peri measures.
    # - A neighbor is linked ONLY if:
    #     * it is within adjacency_thresh (5 s) -> same block, AND
    #     * the CURRENT trial has valid_rt == 1 and no_resp == 0.
    #   Otherwise pre_*/next_* are NaN (block-first trials, block-last
    #   trials, invalid/no-response trials).
    # - Cells are initialized to the sentinel "None" and the assert below
    #   verifies every cell was overwritten (either with a value or NaN).
    # --------------------------------------------------------------------------
    current_cols = trial_data.columns
    for col_name in current_cols:
        # Edit by Ana: explicit object dtype -- on pandas >=2.x, `= "None"` creates a string-dtype
        # column that then rejects the numeric values assigned in the loops below
        trial_data["pre_" + col_name] = pd.Series(["None"] * len(trial_data), dtype="object")
        trial_data["next_" + col_name] = pd.Series(["None"] * len(trial_data), dtype="object")

    # Iterate through each row of the dataframe
    for i in range(len(trial_data)):
        # Check for previous trial (n-1) if it exists and is in the same block
        if i > 0 and (trial_data.loc[i, 'task1_stim_keyResp.stopped'] - trial_data.loc[i-1, 'task1_stim_keyResp.stopped']) <= adjacency_thresh\
        and trial_data.loc[i, 'valid_rt'] == 1 and trial_data.loc[i, 'no_resp'] == 0:  # Edit by Ana: threshold 3 -> 5 (see above)
            for col_name in current_cols:
                trial_data.loc[i, 'pre_' + col_name] = trial_data.loc[i-1, col_name]
        else:
            for col_name in current_cols:
                trial_data.loc[i, 'pre_' + col_name] = np.nan
    for i in range(len(trial_data)):
        # Check for next trial (n+1) if it exists and is in the same block
        if i < len(trial_data)-1 and (trial_data.loc[i+1, 'task1_stim_keyResp.stopped'] - trial_data.loc[i, 'task1_stim_keyResp.stopped']) <= adjacency_thresh\
        and trial_data.loc[i, 'valid_rt'] == 1 and trial_data.loc[i, 'no_resp'] == 0:  # Edit by Ana: threshold 3 -> 5
            for col_name in current_cols:
                trial_data.loc[i, 'next_' + col_name] = trial_data.loc[i+1, col_name]
        else:
            for col_name in current_cols:
                trial_data.loc[i, 'next_' + col_name] = np.nan

    # Check if the string "None" exists anywhere in the DataFrame to make sure all cells were properly populated in the above step
    assert not ((trial_data == "None").any().any()), "Check your data!"

    trial_data.drop(['pre_task1_stim_keyResp.stopped', 'next_task1_stim_keyResp.stopped'], axis = 1, inplace = True)

    # --------------------------------------------------------------------------
    # SECTION 6g: SAVE PER-SUBJECT TRIAL-LEVEL CSV
    # --------------------------------------------------------------------------
    # - "sub" is moved to the first column, then the file is written as
    #   sub-{id}_trial_data.csv. NOTE: full_df at the end of the script is
    #   built from every sub-*_trial_data.csv PRESENT IN THE FOLDER, so if a
    #   subject is removed from the raw data, delete their stale trial_data
    #   CSV before re-running.
    # --------------------------------------------------------------------------
    trial_data["sub"] = sub
    all_cols = list(trial_data.columns)[:-1]
    all_cols.insert(0, "sub")
    trial_data = trial_data[all_cols]

    trial_data.to_csv(f"{output_dataset_path}{output_path}sub-{sub}_trial_data.csv", index=False)

    # --------------------------------------------------------------------------
    # SECTION 6h: SUBJECT-LEVEL SUMMARY MEASURES
    # --------------------------------------------------------------------------
    # - Computed in this order, with progressively stricter filtering:
    #     1. skipped_percent / invalid_rt_percent on ALL trials
    #     2. then keep only valid_rt == 1 trials ->
    #        n_trials_valid, 6_or_more_err (>= 6 commission errors),
    #        acc / acc_con / acc_incon,
    #        RT means (ms): rt_con & rt_incon (correct trials by congruency),
    #        rt_corr (INCONGRUENT correct -- per lab manual, congruency is
    #        held fixed, so rt_corr == rt_incon by design),
    #        rt_err (INCONGRUENT commission errors),
    #        each also as *_log = mean of log(RT in s) x 1000, hence negative
    #     3. then keep only trials whose PREVIOUS trial was clean
    #        (pre_valid_rt == 1, pre_extra_resp == 0, pre_no_resp == 0) ->
    #        the error-monitoring measures, ALL restricted to trials that
    #        FOLLOW AN INCONGRUENT trial (pre_congruent == 0):
    #        pes      post-error slowing, log-RT after errors minus after
    #                 corrects (positive = slowing), correct trials only
    #        pea      accuracy after errors minus after corrects
    #                 (negative = worse after errors)
    #        peri_acc / peri_rt   post-error change in the congruency
    #                 effect (interaction), for accuracy and log-RT
    # - Sample-size caution: with 384 trials these cells can be small; pes
    #   requires post-incongruent-error correct trials, so subjects near the
    #   6-error floor contribute noisy estimates.
    # --------------------------------------------------------------------------
    # Edit by Ana: single condition -> condition loop removed, summary computed once per subject
    condition_data = trial_data
    processing_log["skipped_percent"].append(np.round(condition_data["no_resp"].sum() / len(condition_data) * 100, 3))
    processing_log["invalid_rt_percent"].append(np.round((1 - (sum(condition_data["valid_rt"]) / len(condition_data))) * 100, 3))
    condition_data = condition_data[(condition_data["valid_rt"] == 1)]
    processing_log["n_trials_valid"].append(len(condition_data))  # Edit by Ana: was n_trials inside condition loop
    processing_log["6_or_more_err"].append(1 if len(condition_data[(condition_data["no_resp"] == 0) & (condition_data["accuracy"] == 0)]) >= 6 else 0)
    processing_log["acc"].append(np.round(condition_data.accuracy.mean(), 3))
    processing_log["acc_con"].append(np.round(condition_data[condition_data["congruent"] == 1].accuracy.mean(), 3))
    processing_log["acc_incon"].append(np.round(condition_data[condition_data["congruent"] == 0].accuracy.mean(), 3))
    processing_log["rt_con"].append(np.round(condition_data[(condition_data["congruent"] == 1) & (condition_data["accuracy"] == 1)]["rt"].mean() * 1000, 3))
    processing_log["rt_con_log"].append(np.round(np.log(condition_data[(condition_data["congruent"] == 1) & (condition_data["accuracy"] == 1)]["rt"]).mean() * 1000, 3))
    processing_log["rt_incon"].append(np.round(condition_data[(condition_data["congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"].mean() * 1000, 3))
    processing_log["rt_incon_log"].append(np.round(np.log(condition_data[(condition_data["congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]).mean() * 1000, 3))
    # Edit by Ana: rt_corr/rt_err are INCONGRUENT-ONLY per the lab manual ("keep the
    # congruency fixed" when studying accuracy effects), so rt_corr intentionally
    # equals rt_incon.
    processing_log["rt_corr"].append(np.round(condition_data[(condition_data["congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"].mean() * 1000, 3))
    processing_log["rt_corr_log"].append(np.round(np.log(condition_data[(condition_data["congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]).mean() * 1000, 3))
    processing_log["rt_err"].append(np.round(condition_data[(condition_data["congruent"] == 0) & (condition_data["no_resp"] == 0) & (condition_data["accuracy"] == 0)]["rt"].mean() * 1000, 3))
    processing_log["rt_err_log"].append(np.round(np.log(condition_data[(condition_data["congruent"] == 0) & (condition_data["no_resp"] == 0) & (condition_data["accuracy"] == 0)]["rt"]).mean() * 1000, 3))
    condition_data = condition_data[(condition_data["pre_valid_rt"] == 1) & (condition_data["pre_extra_resp"] == 0) & (condition_data["pre_no_resp"] == 0)]
    processing_log["pes"].append(np.round(
        np.log(
            condition_data[(condition_data["accuracy"] == 1) & (condition_data["pre_accuracy"] == 0) &\
            (condition_data["pre_congruent"] == 0)]["rt"]
        ).mean()\
        - np.log(
            condition_data[(condition_data["accuracy"] == 1) & (condition_data["pre_accuracy"] == 1) &\
            (condition_data["pre_congruent"] == 0)]["rt"]
        ).mean(), 5
    ))
    processing_log["pea"].append(np.round(
         condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["pre_congruent"] == 0)]["accuracy"].mean()\
         - condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["pre_congruent"] == 0)]["accuracy"].mean(), 5
     ))

    processing_log["peri_acc"].append(np.round(
         (
             condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["congruent"] == 0) &\
             (condition_data["pre_congruent"] == 0)]["accuracy"].mean()\
          - condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["congruent"] == 1) &\
             (condition_data["pre_congruent"] == 0)]["accuracy"].mean()
         )\
         - (
             condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["congruent"] == 0) &\
             (condition_data["pre_congruent"] == 0)]["accuracy"].mean()\
          - condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["congruent"] == 1) &\
             (condition_data["pre_congruent"] == 0)]["accuracy"].mean()
         ), 5
     ))

    processing_log["peri_rt"].append(np.round(
         (
             np.log(
             condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["congruent"] == 0) &\
                 (condition_data["pre_congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]
             ).mean()\
          - np.log(
              condition_data[(condition_data["pre_accuracy"] == 0) & (condition_data["congruent"] == 1) &\
              (condition_data["pre_congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]
          ).mean()
         )\
         - (
             np.log(
             condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["congruent"] == 0) &\
             (condition_data["pre_congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]
             ).mean()\
          - np.log(
              condition_data[(condition_data["pre_accuracy"] == 1) & (condition_data["congruent"] == 1) &\
              (condition_data["pre_congruent"] == 0) & (condition_data["accuracy"] == 1)]["rt"]
          ).mean()
           ), 5
     ))

    print(f"sub-{sub} has been processed")

# ------------------------------------------------------------------------------
# SECTION 7: SESSION-LEVEL OUTPUTS
# ------------------------------------------------------------------------------
# - Sanity assert: all processing_log lists must be the same length (one
#   entry per subject) or the summary would be misaligned.
# - summary_{session}_{time}.csv : one row per subject (incl. skipped).
# - full_df_{time}.csv : concatenation of ALL sub-*_trial_data.csv files
#   found in the output folder (see caution in Section 6g).
# ------------------------------------------------------------------------------
# Edit by Ana: sanity check that all processing_log lists are aligned before building the DataFrame
assert len(set(len(v) for v in processing_log.values())) == 1, "processing_log lists misaligned!"
pd.DataFrame(processing_log).to_csv(f"{output_dataset_path}{output_path}summary_{session}_{date_time}.csv", index=False)

list_of_ind_csv = []
for df in sorted([i for i in os.listdir(f"{output_dataset_path}{output_path}") if "sub-" in i]):
    list_of_ind_csv.append(pd.read_csv(f"{output_dataset_path}{output_path}{df}"))
full_df = pd.concat(list_of_ind_csv)
full_df.to_csv(f"{output_dataset_path}{output_path}full_df_{date_time}.csv", index = False)

end = time.time()
print(f"Executed time {np.round(end - start, 2)} s")
