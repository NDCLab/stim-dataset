# ==============================================================================
# stim-behavioral_select-participants.py  --  STAGE 2 of the STIM Behavioral Analysis Pipeline
# ==============================================================================
#   This script loads the summary CSV produced by stim-behavioral_extract-trials.py, applies
#   the usability criteria, and writes stim_valid_behav.csv -- the list of
#   fully usable participants (with their summary measures). Ends with a
#   participant breakdown report: total found / valid / removed, with the
#   removed subject IDs listed per reason.
#
#   USABILITY CRITERIA (applied in this order):
#   1. success == 1        processed in stage 1 (complete data on disk)
#   2. acc >= 0.6          overall accuracy at least 60%
#   3. 6_or_more_err == 1  at least 6 valid commission errors
#   4. not in exclude_id_list (manual exclusions)
#   5. within +/- 3 SD on invalid_rt_percent and skipped_percent
#      (SDs computed over the subjects that survived steps 1-4)
#
#   ADAPTED FROM: thrive-theta-ddm create_valid_behav.py.
#   All deviations from the original are marked with "# Edit by Ana:".
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: IMPORTS
# ------------------------------------------------------------------------------
# - ndclab_py.py must sit in the SAME FOLDER as this script; it provides
#   find_newest_file() and replace_outliers_with_nan_cols().
# ------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import sys
from glob import glob
import os

# Edit by Ana: ndclab_py.py lives next to this script locally, no sys.path hack needed
from ndclab_py import *

# ------------------------------------------------------------------------------
# SECTION 2: CONFIGURATION
# ------------------------------------------------------------------------------
# - path must point to the same output folder stim-behavioral_extract-trials.py writes to.
# - exclude_id_list: manual exclusions as INTEGER IDs (the summary CSV's
#   "sub" column is read as int). IDs listed here that were already dropped
#   earlier (e.g. skipped in stage 1) or that don't exist in the data are
#   simply ignored -- the breakdown reports only the ones actually removed
#   at this step.
# ------------------------------------------------------------------------------
session = "s1_r1"  # Edit by Ana: hardcoded; was sys.argv[1]

# Edit by Ana: derivatives live inside the dataset (BIDS-style), pipeline-named folder
path = f"/Users/anagarci/Documents/NDClab/datasets/stim-dataset/derivatives/stim-behavioral/{session}/"

# Edit by Ana: no unusable/exclusion lists yet; single condition so only one list
exclude_id_list = [410024, 410031, 410077, 410076, 410067]

# ------------------------------------------------------------------------------
# SECTION 3: LOAD NEWEST SUMMARY
# ------------------------------------------------------------------------------
# - Picks the most recently MODIFIED summary_{session}_*.csv in the folder,
#   so after re-running stage 1, this script automatically uses the new run.
# - success == 0 rows (subjects skipped in stage 1) are counted for the
#   breakdown, then removed before any criteria are applied.
# ------------------------------------------------------------------------------
# Edit by Ana: removed stray undefined `matching_files` line and fixed find_newset_file -> find_newest_file
behavior_df = pd.read_csv(find_newest_file(f"{path}/summary_{session}_*.csv"))

# Edit by Ana: track each removal step for the final breakdown report
n_total = behavior_df.shape[0]
skipped_stage1 = behavior_df[behavior_df["success"] == 0]["sub"].tolist()
behavior_df = behavior_df[behavior_df["success"] == 1]

# ------------------------------------------------------------------------------
# SECTION 4: PERFORMANCE CRITERIA & MANUAL EXCLUSIONS
# ------------------------------------------------------------------------------
# - Each filter records the IDs it removes (fail_acc, fail_6err,
#   excluded_manual) before dropping them, purely for the breakdown report.
# - The "Full/New DF length" prints are kept from the original for
#   continuity with older logs; the breakdown at the end is the readable
#   summary.
# ------------------------------------------------------------------------------
# Edit by Ana: single condition -> no _soc/_nonsoc subsetting; one pass over the whole summary
fail_acc = behavior_df[behavior_df["acc"] < 0.6]["sub"].tolist()
behavior_df = behavior_df[behavior_df["acc"] >= 0.6]
fail_6err = behavior_df[behavior_df["6_or_more_err"] != 1]["sub"].tolist()
behavior_df = behavior_df[behavior_df["6_or_more_err"] == 1]

print(f"Full DF length: {behavior_df.shape[0]}")
print(f"Removing subjects {exclude_id_list}")
excluded_manual = behavior_df[behavior_df["sub"].isin(exclude_id_list)]["sub"].tolist()
behavior_df = behavior_df[~behavior_df["sub"].isin(exclude_id_list)].reset_index(drop=True)
print(f"New DF length: {behavior_df.shape[0]} \n")

# ------------------------------------------------------------------------------
# SECTION 5: OUTLIER REMOVAL & SAVE
# ------------------------------------------------------------------------------
# - replace_outliers_with_nan_cols NaNs any value beyond +/- 3 SD of the
#   column mean (computed over the CURRENT, already-filtered sample), then
#   the dropna calls remove those subjects.
# - A subject can exceed the bound on BOTH columns; the [Log] lines count
#   replaced VALUES per column, while removed_outliers counts SUBJECTS --
#   so the [Log] totals can exceed the number of subjects removed.
# - NOTE: removing participants for excessive MISSED responses
#   (the skipped_percent column here) is a criterion the manual says to
#   confirm with the PI. Currently applied; drop "skipped_percent" from the
#   list below if the PI decides against it.
# - Output: stim_valid_behav.csv (same columns as the summary; the row set
#   IS the inclusion list for downstream analyses).
# ------------------------------------------------------------------------------
# criteria-based removals
pre_outlier_subs = behavior_df["sub"].tolist()  # Edit by Ana
behavior_df = replace_outliers_with_nan_cols(behavior_df, ["invalid_rt_percent", "skipped_percent"])
behavior_df = behavior_df.dropna(subset="invalid_rt_percent")
behavior_df = behavior_df.dropna(subset="skipped_percent")
removed_outliers = [s for s in pre_outlier_subs if s not in behavior_df["sub"].tolist()]  # Edit by Ana
print(f"Final DF length: {behavior_df.shape[0]} \n")
behavior_df.to_csv(f"{path}/stim_valid_behav.csv", index=False)  # Edit by Ana: single output CSV

# ------------------------------------------------------------------------------
# SECTION 6: PARTICIPANT BREAKDOWN REPORT
# ------------------------------------------------------------------------------
# READER'S NOTES:
# - Reasons are MUTUALLY EXCLUSIVE and applied in order: a subject removed
#   for low accuracy is not re-counted under later criteria, so the five
#   reason counts always sum to "Removed".
# ------------------------------------------------------------------------------
# Edit by Ana: participant accounting breakdown
n_valid = behavior_df.shape[0]
n_removed = n_total - n_valid
print("=" * 60)
print("PARTICIPANT BREAKDOWN")
print("=" * 60)
print(f"Total participants found:            {n_total}")
print(f"Fully usable (valid):                {n_valid}")
print(f"Removed:                             {n_removed}")
print("-" * 60)
print(f"  Incomplete/missing data (stage 1): {len(skipped_stage1):>3}  {skipped_stage1}")
print(f"  Accuracy < 0.6:                    {len(fail_acc):>3}  {fail_acc}")
print(f"  Fewer than 6 valid errors:         {len(fail_6err):>3}  {fail_6err}")
print(f"  Manual exclusion list:             {len(excluded_manual):>3}  {excluded_manual}")
print(f"  Outlier (+/-3 SD on invalid_rt/skipped %): {len(removed_outliers):>3}  {removed_outliers}")
print("=" * 60)

print("Processing has been completed; 1 CSV was created!")