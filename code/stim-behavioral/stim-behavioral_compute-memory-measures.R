# ==============================================================================
# stim-behavioral_compute-memory-measures.R  --  STAGE 4 of the STIM Behavioral Analysis Pipeline
# ==============================================================================
#   For each Python-valid participant (gated on stim_valid_behav.csv, same
#   as stage 1), this script loads the organized flankerDat/surpriseDat CSVs,
#   applies four inclusion criteria (Section 4b), computes flanker
#   performance measures and the memory hit-rate measures (old faces from
#   incongruent ERROR trials vs incongruent CORRECT trials), assigns the
#   between-subject condition, and writes three CSVs to
#   {proje_wd}/derivatives/stim-memory/s1_r1/stat_output/ :
#     - stim-behavioral_memory-measures.csv              one row per included participant
#     - stim-behavioral_memory-summary-by-condition.csv  per-condition mean/SD/n
#     - stim-behavioral_memory-exclusions.csv            one row per excluded participant + reason
#   It also prints the per-condition summary and an exclusion breakdown.
#
#   ADAPTED FROM: stim_measure_computations.R (by Kianoosh Hosseini)
#   All deviations from the original are marked with "# Edit by Ana:" (except for section banners
#   and reader's notes.)
#
#   INCLUSION CRITERIA (nested, in this order -- a participant excluded at one
#   step is never evaluated on the later ones, and the four counters are
#   therefore mutually exclusive):
#   1. <= 20% of surprise memory trials removed for rt <= 0.2 s
#   2. flanker accuracy >= 0.6 (fast rt < 0.15 s trials excluded; missing
#      responses count as accuracy 0 -- per lab manual, EDIT 2026-08-17)
#   3. >= 8 legit incongruent errors in the flanker task
#      (note: one inline comment below says "at least 6" -- the CODE says 8)
#   4. >= 8 of those incongruent-error faces present in the (rt-filtered)
#      surprise data
#
#   KEY MEASURE DEFINITIONS:
#   error_hitRate    among OLD faces from legit incongruent-ERROR flanker
#                    trials that appear in the surprise task: proportion
#                    judged "Old" (hits / (hits + misses))
#   correct_hitRate  same, for faces from legit incongruent-CORRECT trials
#   hitRate_error_minus_correct  the study's main memory contrast
#   flankEff_*       incongruent minus congruent (acc / RT / log-RT)
#
# !! REMAINING CAUTIONS / OPEN DECISIONS:
#   - log-RT here is mean(log(1 + rt)); the Python pipeline and the manual's
#     PES formulas use log(rt). Left unchanged (the manual does not govern
#     the flankEff measures) -- decision pending with the script's author.
#   - Inclusion criteria 3 and 4 (>= 8) are a study-specific tightening of
#     the manual's >= 6; intentional for the memory analysis, but the two
#     pipelines can therefore produce DIFFERENT valid-participant lists.
#   - The four exclusion counters are now printed AND saved (Section 7),
#     with the dropped participant IDs per criterion. All four share the
#     "num_of_participants_removed_based_on_*" naming.
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: LIBRARIES & WORKING DIRECTORY
# ------------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(stringr)
library(psycho)


#Working directory should be the Psychopy experiment directory.
proje_wd <- "/Users/anagarci/Documents/NDClab/datasets/stim-dataset"
setwd(proje_wd)

# ------------------------------------------------------------------------------
# SECTION 2: INPUT / OUTPUT PATHS & FILENAME SUFFIXES
# ------------------------------------------------------------------------------
# - Reads raw CSVs only to recover each participant's id; the actual data
#   come from the organized csv_output files written by stage 3 (prepare-memory-data).
# - output: derivatives/stim-memory/s1_r1/stat_output/ (must ALREADY EXIST).
# - The raw-file listing uses recursive = TRUE (same as stage 3): raw CSVs
#   live in per-subject subfolders under input_raw_path, the layout the
#   Python pipeline reads. Returned paths are relative and include the
#   subfolder; the ID-extraction and per-file path lines handle that.
# ------------------------------------------------------------------------------
# Defining the input and output folders.
input_raw_path <- paste(proje_wd, "sourcedata", "raw", "s1_r1", "psychopy", sep ="/", collapse = NULL) # input data directory
input_organized_path <- paste(proje_wd, "derivatives", "stim-memory", "s1_r1", "csv_output", sep ="/", collapse = NULL) # Edit by Ana: matches prepare-memory output (was derivatives/psychopy/csv_output)
output_path <- paste(proje_wd, "derivatives", "stim-memory", "s1_r1", "stat_output", sep ="/", collapse = NULL) # Edit by Ana: was derivatives/psychopy/stat_output
flanker_csv_fileName <- "_stim_flankerDat.csv" # each output csv file will have this on its filename
surprise_csv_fileName <- "_stim_surpriseDat.csv" # each output csv file will have this on its filename


# ------------------------------------------------------------------------------
# SECTION 3: FILE LIST, EMPTY SUMMARY FRAME & EXCLUSION COUNTERS
# ------------------------------------------------------------------------------
# - main_df: one row per INCLUDED participant, 21 summary columns.
# - The four counters (and matching excluded_ids_* vectors) tally exclusions
#   per criterion, mutually exclusive; printed and saved in Section 7.
# ------------------------------------------------------------------------------
## creating a list of all raw data csv files in the input folder.
raw_datafiles_list <- c() # an empty list that will be filled in the next "for" loop!
csvSelect <- list.files(input_raw_path, pattern = ".csv", recursive = TRUE) # Edit by Ana: recursive = TRUE -- raw CSVs live in per-subject subfolders (sub-410002/...), same layout the Python pipeline reads. Paths are returned relative to input_raw_path; the ID-extraction and per-file path lines below handle the subfolder prefix correctly.
for (i in 1:length(csvSelect)){
  temp_for_file <- ifelse (str_detect(csvSelect[i], "stim_v1_", negate = FALSE), 1, 0)
  if (temp_for_file == 1){
    temp_list <- csvSelect[i]
    raw_datafiles_list <- c(raw_datafiles_list, temp_list)
  }
}

# Edit by Ana: restrict processing to participants deemed FULLY VALID by the Python
# behavioral pipeline (same gating as in stim-behavioral_prepare-memory-data.R; see note there).
valid_behav_file <- "/Users/anagarci/Documents/NDClab/datasets/stim-dataset/derivatives/stim-behavioral/s1_r1/stim_valid_behav.csv"  # Edit by Ana:
if (!file.exists(valid_behav_file)) {  # Edit by Ana:
  stop(paste("Valid-subject list not found:", valid_behav_file, "-- run the Python pipeline (stim-behavioral_extract-trials.py, then stim-behavioral_select-participants.py) first."))  # Edit by Ana:
}  # Edit by Ana:
valid_subs <- read.csv(valid_behav_file)$sub  # Edit by Ana:
file_sub_ids <- as.numeric(gsub("^sub-(\\d+).*", "\\1", raw_datafiles_list))  # Edit by Ana:
n_before_valid_filter <- length(raw_datafiles_list)  # Edit by Ana:
raw_datafiles_list <- raw_datafiles_list[file_sub_ids %in% valid_subs]  # Edit by Ana:
n_gated_in <- length(raw_datafiles_list)  # Edit by Ana: number of Python-valid subjects fed into the loop (for the exclusion breakdown)
print(paste("Restricting to Python-valid subjects:", length(raw_datafiles_list), "of", n_before_valid_filter, "data files kept."))  # Edit by Ana:
# Creating the main empty dataframe that will be filled with the data from the loop below:
main_df <- setNames(data.frame(matrix(ncol = 21, nrow = 0)), c("participant_id", "congAcc", "incongAcc",
                                                               "incongruent_dat_meanRT", "errorDat_meanRT", "congruent_dat_meanRT", "corrDat_meanRT",
                                                               "congCorr_meanRT", "incongCorr_meanRT", "congCorr_logMeanRT",
                                                               "congErr_meanRT", "incongErr_meanRT", "congErr_logMeanRT", "incongErr_logMeanRT",
                                                               "incongCorr_logMeanRT", "flankEff_meanACC", "flankEff_meanRT", "flankEff_logMeanRT",
                                                               "error_hitRate", "correct_hitRate", "hitRate_error_minus_correct"))


# Counters for the number of excluded people based on each criterion
num_of_participants_removed_based_on_memory_surp_trial_removal <- 0 # Will be updated in the loop below
num_of_participants_removed_based_on_accuracy <- 0 # Will be updated in the loop below
num_of_participants_removed_based_on_incong_error_num <- 0 # Will be updated in the loop below
num_of_participants_removed_based_on_num_incong_error_faces_in_memory_surp <- 0 # Edit by Ana: renamed (added "of_") so all four counters share one pattern; was num_participants_...

# Edit by Ana: also record WHICH participant IDs each criterion drops, for the breakdown report at the end
excluded_ids_memory_surp_trial_removal <- c()
excluded_ids_accuracy <- c()
excluded_ids_incong_error_num <- c()
excluded_ids_num_incong_error_faces_in_memory_surp <- c()

# ==============================================================================
# SECTION 4: PER-PARTICIPANT LOOP
# ==============================================================================
# Looping over all participants
for (subject in 1:length(raw_datafiles_list)){
  
  # ----------------------------------------------------------------------------
  # SECTION 4a: LOAD THIS PARTICIPANT'S ORGANIZED DATA
  # ----------------------------------------------------------------------------
  #for this participant, find the raw csv file
  psychopy_file <- paste(input_raw_path,raw_datafiles_list[subject], sep = "/", collapse = NULL)
  
  #read in the data for this participant, establish id, and remove extraneous variables
  psychopyDat <- read.csv(file = psychopy_file, stringsAsFactors = FALSE, na.strings=c("", "NA"))
  participant_id <- psychopyDat$id[1]
  
  # Load this participant's flanker and surprise data frames
  flanker_name <- paste0(participant_id, flanker_csv_fileName, sep = "", collapse = NULL)
  surprise_name <- paste0(participant_id, surprise_csv_fileName, sep = "", collapse = NULL)
  flanker_df <- read.csv(file = paste(input_organized_path, flanker_name, sep = "/", collapse = NULL), stringsAsFactors = FALSE, na.strings=c("", "NA"))
  surprise_df <- read.csv(file = paste(input_organized_path, surprise_name, sep = "/", collapse = NULL), stringsAsFactors = FALSE, na.strings=c("", "NA"))
  
  # Edit by Ana: manual-compliant trial subsets, used below.
  #   flanker_df_acc -- for ACCURACY: excludes fast trials (rt < 0.15 s) but KEEPS
  #                     no-response trials (rt NA, accuracy 0), per the manual.
  #   flanker_df_rt  -- for RT means: excludes fast trials AND no-response trials
  #                     (subset() drops rows where the condition is NA).
  flanker_df_acc <- subset(flanker_df, is.na(current_trial_rt) | current_trial_rt >= 0.15)  # Edit by Ana:
  flanker_df_rt <- subset(flanker_df, current_trial_rt >= 0.15)  # Edit by Ana:
  
  # ----------------------------------------------------------------------------
  # SECTION 4b: INCLUSION CRITERIA (nested ifs; see header for the list)
  # ----------------------------------------------------------------------------
  # - Criterion 1: surprise trials with rt <= 0.2 s count as "removed"; if
  #   more than 20% (rounded) of surprise trials are removed, the
  #   participant is excluded as inattentive.
  # - Criterion 3: "legit incongruent errors" = incongruent trials with
  #   accuracy == 0 AND legitResponse == 1. The inline comment says "at
  #   least 8"-worthy check but an earlier comment mentions 6 -- the CODE
  #   requires >= 8.
  # - Criterion 4: for each legit incongruent-error face, the (rt-filtered)
  #   surprise data must contain EXACTLY ONE row with that face for it to
  #   count (duplicate faces would be counted as absent).
  # - NOTE: surprise_df is overwritten by its rt-filtered version inside
  #   criterion 4's block, and all later face lookups use the filtered data.
  # ----------------------------------------------------------------------------
  # removing participants based on whether they just pressed the keys without actually paying attention to the task.
  # We check this by "keep_surp_memory_trial_based_on_rt" and "keep_surp_friendly_trial_based_on_rt" variables
  # in the surprise_df data frame. O means that they have responded faster than 200 ms during the given trial
  # and therefore, we need to remove that surprise trial.
  # In this study, we will remove people based on the surprise memory task not surprise friendly task.
  # If more than 20% of surprise trials are removed, we exclude that participant.
  number_of_removed_trials_in_the_memory_surp_based_on_rt <- nrow(surprise_df) - (sum(surprise_df$keep_surp_memory_trial_based_on_rt))
  number_of_faces_in_surp_memory_task <- nrow(surprise_df)
  twenty_percent_threshold <- round(0.2 * number_of_faces_in_surp_memory_task)
  
  if (number_of_removed_trials_in_the_memory_surp_based_on_rt <= twenty_percent_threshold){ # Participants who have less than twenty_percent_threshold surprise
    # trials removed, will be included.
    
    # Check to see if this participant has the flanker accuracy above 60%
    if (mean(as.numeric(flanker_df_acc$current_trial_accuracy), na.rm = TRUE ) >= 0.6){ # Edit by Ana: was flanker_df; now excludes fast trials, includes no-resp as 0 (manual)
      
      # check to see if the participant has at least 8 legit incongruent errors or not.
      incong_flankerDat <- filter(flanker_df, current_trial_congruency == 0)
      cong_flankerDat <- filter(flanker_df, current_trial_congruency == 1)
      error_incong_flankerDat <- filter(incong_flankerDat, current_trial_accuracy == 0)
      legit_error_incong_flankerDat <- filter(error_incong_flankerDat, current_trial_legitResponse == 1)
      if (nrow(legit_error_incong_flankerDat) >= 8 ){
        
        # Checking to see if there are at least 6 incongruent error faces in the surprise_df of this participant.
        # First we need to remove the trials that are marked based on rt!
        surprise_df <- filter(surprise_df, keep_surp_memory_trial_based_on_rt == 1)
        num_incong_error_faces_in_surp <- 0
        # Counting the number of legit incong error faces available in surprise_df!
        for (rr in 1:nrow(legit_error_incong_flankerDat)){
          temp_face <- legit_error_incong_flankerDat$current_trial_face[rr]
          temp_for_surp <- filter(surprise_df, face == temp_face)
          errorFace_exist_in_surpDat <- ifelse(nrow(temp_for_surp) == 1, 1,0)
          if (errorFace_exist_in_surpDat == 1){
            num_incong_error_faces_in_surp <- num_incong_error_faces_in_surp + 1
          }
        } # Closing the loop that counts the number of incong error faces available in surprise_df!
        if (num_incong_error_faces_in_surp >= 8){
          
          # --------------------------------------------------------------------
          # SECTION 4c: FLANKER PERFORMANCE MEASURES (included participants)
          # --------------------------------------------------------------------
          # - Accuracy by congruency, then mean RTs for every crossing of
          #   correct/error x congruent/incongruent, raw and log(1 + rt).
          # - Per lab manual (EDIT 2026-08-17): accuracy means use
          #   flanker_df_acc (fast trials excluded, no-resp counted as 0);
          #   RT means use flanker_df_rt (fast and no-resp excluded).
          # - flankEff_* = incongruent minus congruent, so positive =
          #   worse/slower on incongruent (the usual flanker effect for RT;
          #   for accuracy a NEGATIVE value is the typical direction).
          # --------------------------------------------------------------------
          # Edit by Ana: accuracy computed from flanker_df_acc (fast trials excluded,
          # no-resp counted as 0), per manual. incong_/cong_flankerDat (unfiltered)
          # remain in use above for the legit-error count, where the
          # legitResponse == 1 filter already handles validity.
          incong_flankerDat_acc <- filter(flanker_df_acc, current_trial_congruency == 0)  # Edit by Ana:
          cong_flankerDat_acc <- filter(flanker_df_acc, current_trial_congruency == 1)  # Edit by Ana:
          incongAcc <- mean(as.numeric(incong_flankerDat_acc$current_trial_accuracy), na.rm = TRUE)  # Edit by Ana: was incong_flankerDat
          congAcc <- mean(as.numeric(cong_flankerDat_acc$current_trial_accuracy), na.rm = TRUE)  # Edit by Ana: was cong_flankerDat
          
          
          # subset the data for correct and error trials, separately for congruent and incongruent trials, creating new data frames for each
          # Edit by Ana: all RT subsets below built from flanker_df_rt (was flanker_df),
          # so fast (< 0.15 s) and no-response trials are excluded from RT means (manual).
          corrDat <- flanker_df_rt[flanker_df_rt$current_trial_accuracy == 1,]  # Edit by Ana:
          corrDat_meanRT <- mean(corrDat$current_trial_rt, na.rm = TRUE)
          
          congruent_dat <- flanker_df_rt[flanker_df_rt$current_trial_congruency == 1,]  # Edit by Ana:
          congruent_dat_meanRT <- mean(congruent_dat$current_trial_rt, na.rm = TRUE)
          
          cong_corrDat <- corrDat[corrDat$current_trial_congruency == 1,]
          incong_corrDat <- corrDat[corrDat$current_trial_congruency == 0,]
          
          errorDat <- flanker_df_rt[flanker_df_rt$current_trial_accuracy == 0,]  # Edit by Ana:
          errorDat_meanRT <- mean(errorDat$current_trial_rt, na.rm = TRUE)
          
          incongruent_dat <- flanker_df_rt[flanker_df_rt$current_trial_congruency == 0,]  # Edit by Ana:
          incongruent_dat_meanRT <- mean(incongruent_dat$current_trial_rt, na.rm = TRUE)
          
          
          cong_errorDat <- errorDat[errorDat$current_trial_congruency == 1,]
          incong_errorDat <- errorDat[errorDat$current_trial_congruency == 0,]
          #for correct trials, compute mean RT (raw and log-corrected)
          congCorr_meanRT <- mean(cong_corrDat$current_trial_rt, na.rm = TRUE)
          incongCorr_meanRT <- mean(incong_corrDat$current_trial_rt, na.rm = TRUE)
          
          congErr_meanRT <- mean(cong_errorDat$current_trial_rt, na.rm = TRUE)
          incongErr_meanRT <- mean(incong_errorDat$current_trial_rt, na.rm = TRUE)
          
          congCorr_logMeanRT <- mean(log((1+cong_corrDat$current_trial_rt)), na.rm = TRUE)
          incongCorr_logMeanRT <- mean(log((1+incong_corrDat$current_trial_rt)), na.rm = TRUE)
          
          congErr_logMeanRT <- mean(log((1+cong_errorDat$current_trial_rt)), na.rm = TRUE)
          incongErr_logMeanRT <- mean(log((1+incong_errorDat$current_trial_rt)), na.rm = TRUE)
          
          # compute flanker-effect scores for accuracy, RT, log-RT
          flankEff_meanACC <- incongAcc - congAcc
          flankEff_meanRT <- incongCorr_meanRT - congCorr_meanRT
          flankEff_logMeanRT <- incongCorr_logMeanRT - congCorr_logMeanRT
          
          # --------------------------------------------------------------------
          # SECTION 4d: MEMORY HIT RATES -- INCONG. ERROR vs CORRECT FACES
          # --------------------------------------------------------------------
          # - From here on, only LEGIT trials are used: incong_corrDat is
          #   narrowed to legitResponse == 1, and the error set reuses
          #   legit_error_incong_flankerDat from criterion 3.
          # - Four near-identical loops count, among faces that appear
          #   exactly once in the rt-filtered surprise data:
          #     error faces judged Old  -> hits      (identified_as_new == 0)
          #     error faces judged New  -> misses    (identified_as_new == 1)
          #     correct faces judged Old -> hits
          #     correct faces judged New -> misses
          #   ("Old" is the correct answer for all these faces, since they
          #   were shown during the flanker task.)
          # - hitRate = hits / (hits + misses). If a participant has zero
          #   countable faces in a category the denominator is 0 -> NaN.
          # - hitRate_error_minus_correct > 0 means better memory for faces
          #   seen on error trials than on correct trials.
          # --------------------------------------------------------------------
          legit_corr_incong_flankerDat <- filter(incong_corrDat, current_trial_legitResponse == 1)
          incong_corrDat <- legit_corr_incong_flankerDat
          
          error_incong_flankerDat <- legit_error_incong_flankerDat
          ### number of incong error faces identified as old
          num_incong_errorFaces_reported_old <- 0 # this is the number of incongruent error faces that they report as OLD and will be updated in the loop below:
          for (iii in 1:nrow(error_incong_flankerDat)){
            if (error_incong_flankerDat$current_trial_legitResponse[iii] == 1){ # if this trial is legit
              temp_face_from_flanker <- error_incong_flankerDat$current_trial_face[iii]
              temp_face_row_in_surp <- filter(surprise_df, face == temp_face_from_flanker)
              if (nrow(temp_face_row_in_surp) == 1){ # if old incong error face exist in the surprise_df
                if (temp_face_row_in_surp$identified_as_new == 0){ # the face is identified as old
                  num_incong_errorFaces_reported_old <- num_incong_errorFaces_reported_old + 1 # The number of incongruent error faces that they report as OLD
                }
              }
            }
          } # closing the loop over error incong flankerDat
          
          ### number of incong error faces identified as new
          num_incong_errorFaces_reported_new <- 0 # this is the number of incongruent error faces that they report as new and will be updated in the loop below:
          for (iii in 1:nrow(error_incong_flankerDat)){
            if (error_incong_flankerDat$current_trial_legitResponse[iii] == 1){ # if this trial is legit
              temp_face_from_flanker <- error_incong_flankerDat$current_trial_face[iii]
              temp_face_row_in_surp <- filter(surprise_df, face == temp_face_from_flanker)
              if (nrow(temp_face_row_in_surp) == 1){ # if old incong error face exist in the surprise_df
                if (temp_face_row_in_surp$identified_as_new == 1){ # the face is identified as New
                  num_incong_errorFaces_reported_new <- num_incong_errorFaces_reported_new + 1 # The number of incongruent error faces that they report as New
                }
              }
            }
          } # closing the loop over error incong flankerDat
          
          
          ### number of incong correct faces identified as old
          num_incong_correctFaces_reported_old <- 0 # this is the number of incongruent correct faces that they report as Old and will be updated in the loop below:
          for (iii in 1:nrow(incong_corrDat)){
            if (incong_corrDat$current_trial_legitResponse[iii] == 1){ # if this trial is legit
              temp_face_from_flanker <- incong_corrDat$current_trial_face[iii]
              temp_face_row_in_surp <- filter(surprise_df, face == temp_face_from_flanker)
              if (nrow(temp_face_row_in_surp) == 1){ # if old incong correct face exist in the surprise_df
                if (temp_face_row_in_surp$identified_as_new == 0){ # the face is identified as Old
                  num_incong_correctFaces_reported_old <- num_incong_correctFaces_reported_old + 1 # The number of incongruent correct faces that they report as Old
                }
              }
            }
          } # closing the loop over correct incong flankerDat
          
          ### number of incong correct faces identified as new
          num_incong_correctFaces_reported_new <- 0 # this is the number of incongruent correct faces that they report as new and will be updated in the loop below:
          for (iii in 1:nrow(incong_corrDat)){
            if (incong_corrDat$current_trial_legitResponse[iii] == 1){ # if this trial is legit
              temp_face_from_flanker <- incong_corrDat$current_trial_face[iii]
              temp_face_row_in_surp <- filter(surprise_df, face == temp_face_from_flanker)
              if (nrow(temp_face_row_in_surp) == 1){ # if old incong correct face exist in the surprise_df
                if (temp_face_row_in_surp$identified_as_new == 1){ # the face is identified as New
                  num_incong_correctFaces_reported_new <- num_incong_correctFaces_reported_new + 1 # The number of incongruent correct faces that they report as New
                }
              }
            }
          } # closing the loop over correct incong flankerDat
          
          
          incong_error_hit_num <- num_incong_errorFaces_reported_old
          incong_error_miss_num <- num_incong_errorFaces_reported_new
          error_hitRate <- (incong_error_hit_num) / ((incong_error_hit_num) + incong_error_miss_num) # hit rate
          
          
          incong_correct_hit_num <- num_incong_correctFaces_reported_old
          incong_correct_miss_num <- num_incong_correctFaces_reported_new
          correct_hitRate <- (incong_correct_hit_num) / ((incong_correct_hit_num) + incong_correct_miss_num) # hit rate
          
          hitRate_error_minus_correct <- error_hitRate - correct_hitRate
          #### filling the main data frame
          main_df[nrow(main_df) + 1,] <-c(participant_id, congAcc, incongAcc,
                                          incongruent_dat_meanRT, errorDat_meanRT, congruent_dat_meanRT, corrDat_meanRT,
                                          congCorr_meanRT, incongCorr_meanRT, congCorr_logMeanRT,
                                          congErr_meanRT, incongErr_meanRT, congErr_logMeanRT, incongErr_logMeanRT,
                                          incongCorr_logMeanRT, flankEff_meanACC, flankEff_meanRT, flankEff_logMeanRT,
                                          error_hitRate,
                                          correct_hitRate, hitRate_error_minus_correct)
          
          
        } else { # If a participant has been excluded because they had less than 8 legit incong error faces in the surprise memory task, we add 1 to the counter below.
          num_of_participants_removed_based_on_num_incong_error_faces_in_memory_surp <- num_of_participants_removed_based_on_num_incong_error_faces_in_memory_surp + 1 # Edit by Ana: renamed counter (added "of_")
          excluded_ids_num_incong_error_faces_in_memory_surp <- c(excluded_ids_num_incong_error_faces_in_memory_surp, participant_id) # Edit by Ana
        }
      } else { # If a participant has been excluded because they had less than 8 legit incong errors, we add 1 to the counter below.
        num_of_participants_removed_based_on_incong_error_num <- num_of_participants_removed_based_on_incong_error_num + 1
        excluded_ids_incong_error_num <- c(excluded_ids_incong_error_num, participant_id) # Edit by Ana
      }
    } else { # If a participant has been excluded because they had less than 60% flanker accuracy, we add 1 to the counter below.
      num_of_participants_removed_based_on_accuracy <- num_of_participants_removed_based_on_accuracy + 1
      excluded_ids_accuracy <- c(excluded_ids_accuracy, participant_id) # Edit by Ana
    }
  } else { # If a participant has been excluded because they had more than 20% surp trial removed, we add 1 to the counter below.
    num_of_participants_removed_based_on_memory_surp_trial_removal <- num_of_participants_removed_based_on_memory_surp_trial_removal + 1
    excluded_ids_memory_surp_trial_removal <- c(excluded_ids_memory_surp_trial_removal, participant_id) # Edit by Ana
  }
} 

# Closing the loop for each participant

# ------------------------------------------------------------------------------
# SECTION 5: CONDITION ASSIGNMENT (hardcoded ID lists)
# ------------------------------------------------------------------------------
# - Between-subject condition (1, 2, or 3) is assigned from hardcoded
#   participant lists. ANY NEW PARTICIPANT MUST BE ADDED to the correct
#   list, otherwise their condition is NA and they silently drop out of
#   the group summary (though they remain in the output CSV).
# ------------------------------------------------------------------------------

main_df <- main_df %>%
  mutate(
    condition = case_when(
      participant_id %in% c(410004, 410013, 410009, 410012, 410018, 410010, 410024, 410027, 410039, 410047, 410048, 410050, 410040, 410044, 410056, 410061, 410057, 410063, 410071, 410077, 410067) ~ 1,
      participant_id %in% c(410002, 410007, 410015, 410008, 410019, 410022, 410025, 410028, 410037, 410032, 410042, 410049, 410033, 410045, 410059, 410058, 410062, 410069, 410072, 410073) ~ 2,
      participant_id %in% c(410005, 410014, 410011, 410020, 410016, 410030, 410026, 410036, 410034, 410038, 410051, 410052, 410043, 410054, 410055, 410060, 410064, 410066, 410074) ~ 3,
      TRUE ~ NA_real_
    )
  )

# ------------------------------------------------------------------------------
# SECTION 6: GROUP SUMMARY & SAVE
# ------------------------------------------------------------------------------
# - Prints AND saves mean, SD, and n (subject count) of the three hit-rate
#   measures per condition. Any condition-NA participants form their own
#   row -- an n against a blank condition flags an ID missing from the
#   Section 5 lists.
# - Writes two CSVs to stat_output/: main_df (one row per included
#   participant) and summary_stats (one row per condition). NOTE: both are
#   OVERWRITTEN on every run; there is no timestamp in the filenames.
# ------------------------------------------------------------------------------
summary_stats <- main_df %>%
  group_by(condition) %>%
  summarise(
    n = n(),  # Edit by Ana: number of subjects contributing to each condition's stats
    across(
      c(error_hitRate, correct_hitRate, hitRate_error_minus_correct),
      list(mean = ~mean(.x, na.rm = TRUE),
           sd   = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

print(summary_stats)

####################
# Save the dataset
proc_fileName <- "stim-behavioral_memory-measures.csv" # output filename
# write the extracted and computed summary scores to disk
write.csv(main_df, paste(output_path, proc_fileName, sep = "/", collapse = NULL), row.names=FALSE)
# Edit by Ana: also save the per-condition group summary (mean/SD/n) to its own CSV
summary_fileName <- "stim-behavioral_memory-summary-by-condition.csv"
write.csv(summary_stats, paste(output_path, summary_fileName, sep = "/", collapse = NULL), row.names=FALSE)
##################

# ------------------------------------------------------------------------------
# SECTION 7: EXCLUSION BREAKDOWN (memory-stage criteria)  -- Edit by Ana
# ------------------------------------------------------------------------------
# - Accounts for the participants dropped by THIS script's criteria, i.e. those
#   who were Python-valid (gated in) but did not make it into the memory output.
# - Reasons are MUTUALLY EXCLUSIVE and applied in the loop's nested order, so
#   the four reason counts always sum to (gated in - included), mirroring the
#   participant breakdown in stim-behavioral_select-participants.py.
# - Prints the report to the console AND saves it as a CSV in stat_output/.
# ------------------------------------------------------------------------------
n_included <- nrow(main_df)
n_excluded <- n_gated_in - n_included

cat(paste0(
  strrep("=", 60), "\n",
  "MEMORY-STAGE EXCLUSION BREAKDOWN\n",
  strrep("=", 60), "\n",
  sprintf("Python-valid subjects (gated in):     %d\n", n_gated_in),
  sprintf("Included in memory output:            %d\n", n_included),
  sprintf("Removed by memory-stage criteria:     %d\n", n_excluded),
  strrep("-", 60), "\n",
  sprintf("  >20%% surprise trials removed (rt):   %2d  %s\n",
          num_of_participants_removed_based_on_memory_surp_trial_removal,
          paste(excluded_ids_memory_surp_trial_removal, collapse = ", ")),
  sprintf("  Flanker accuracy < 0.6:              %2d  %s\n",
          num_of_participants_removed_based_on_accuracy,
          paste(excluded_ids_accuracy, collapse = ", ")),
  sprintf("  < 8 legit incongruent errors:        %2d  %s\n",
          num_of_participants_removed_based_on_incong_error_num,
          paste(excluded_ids_incong_error_num, collapse = ", ")),
  sprintf("  < 8 error faces in surprise data:    %2d  %s\n",
          num_of_participants_removed_based_on_num_incong_error_faces_in_memory_surp,
          paste(excluded_ids_num_incong_error_faces_in_memory_surp, collapse = ", ")),
  strrep("=", 60), "\n"
))

# Build a tidy exclusions table (one row per excluded participant) and save it.
exclusion_df <- rbind(
  if (length(excluded_ids_memory_surp_trial_removal) > 0)
    data.frame(participant_id = excluded_ids_memory_surp_trial_removal, reason = "surprise_trials_removed_gt20pct"),
  if (length(excluded_ids_accuracy) > 0)
    data.frame(participant_id = excluded_ids_accuracy, reason = "flanker_accuracy_lt_0.6"),
  if (length(excluded_ids_incong_error_num) > 0)
    data.frame(participant_id = excluded_ids_incong_error_num, reason = "lt_8_legit_incong_errors"),
  if (length(excluded_ids_num_incong_error_faces_in_memory_surp) > 0)
    data.frame(participant_id = excluded_ids_num_incong_error_faces_in_memory_surp, reason = "lt_8_error_faces_in_surprise")
)
if (is.null(exclusion_df)) {  # no one excluded -> still write an empty table with headers
  exclusion_df <- data.frame(participant_id = integer(0), reason = character(0))
}
exclusion_fileName <- "stim-behavioral_memory-exclusions.csv"
write.csv(exclusion_df, paste(output_path, exclusion_fileName, sep = "/", collapse = NULL), row.names=FALSE)
