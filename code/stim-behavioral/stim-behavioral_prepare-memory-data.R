# ==============================================================================
# stim-behavioral_prepare-memory-data.R  --  STAGE 3 of the STIM Behavioral Analysis Pipeline
# ==============================================================================
#   For every raw PsychoPy CSV ("stim_v1_" files) belonging to a subject
#   deemed FULLY VALID by the Python pipeline, this script extracts and recodes
#   (a) the 384 main flanker trials and (b) the surprise memory-test trials,
#   and writes two organized CSVs per participant into
#   {proje_wd}/derivatives/stim-memory/s1_r1/csv_output/:
#     - {id}_stim_flankerDat.csv    one row per flanker trial
#     - {id}_stim_surpriseDat.csv   one row per surprise memory trial
#
#   ADAPTED FROM: stim_data_organizer.R (by Kianoosh Hosseini)
#   All deviations from the original are marked with "# Edit by Ana:" (except for section banners
#   and reader's notes.)
#
#   OUTPUT COLUMN CONVENTIONS:
#   flankerDat:  current_trial_accuracy    1 = correct, 0 = error or no
#                                          response (per lab manual)
#                current_trial_congruency  1 = congruent, 0 = incongruent
#                current_trial_rt          seconds, first response only
#                current_trial_responded   1 = responded, 0 = no response
#                current_trial_legitResponse 1 = responded AND rt >= 0.15 s
#                current_trial_resp_nums   number of key presses in the trial
#   surpriseDat: is_new             0 = Old face (shown during flanker),
#                                   1 = New face (task's ground truth)
#                identified_as_new  1 = participant judged face as New
#                keep_surp_memory_trial_based_on_rt  1 = rt > 0.2 s
# ==============================================================================


# ------------------------------------------------------------------------------
# SECTION 1: LIBRARIES & WORKING DIRECTORY
# ------------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(stringr)

proje_wd <- "/Users/anagarci/Documents/NDClab/datasets/stim-dataset"
setwd(proje_wd)

# ------------------------------------------------------------------------------
# SECTION 2: INPUT / OUTPUT PATHS & FILENAME SUFFIXES
# ------------------------------------------------------------------------------
# - input:  sourcedata/raw/s1_r1/psychopy/     (raw PsychoPy CSVs)
# - output: derivatives/stim-memory/s1_r1/csv_output/  (must ALREADY EXIST -- the
#   script does not create it; write.csv fails otherwise)
# ------------------------------------------------------------------------------
# Defining the input and output folders.
input_path <- paste(proje_wd, "sourcedata", "raw", "s1_r1", "psychopy", sep ="/", collapse = NULL) # input data directory
output_path <- paste(proje_wd, "derivatives", "stim-memory", "s1_r1", "csv_output", sep ="/", collapse = NULL) # Edit by Ana: was derivatives/psychopy/csv_output; now pipeline-named + session-nested
flanker_csv_fileName <- "_stim_flankerDat.csv" # each output csv file will have this on its filename
surprise_csv_fileName <- "_stim_surpriseDat.csv" # each output csv file will have this on its filename


# ------------------------------------------------------------------------------
# SECTION 3: BUILD THE LIST OF RAW DATA FILES
# ------------------------------------------------------------------------------
# - Keeps every .csv in input_path whose name contains "stim_v1_".
# - list.files() uses recursive = TRUE: the raw CSVs live in per-subject
#   subfolders (sub-410002/sub-410002_...csv), the same layout the Python
#   pipeline reads, so both scripts point at the identical input_path. The
#   returned paths are RELATIVE to input_path (they include the subfolder),
#   which the ID-extraction and per-file path lines below handle correctly.
# ------------------------------------------------------------------------------
## creating a list of all data csv files in the input folder.
datafiles_list <- c() # an empty list that will be filled in the next "for" loop!
csvSelect <- list.files(input_path, pattern = ".csv", recursive = TRUE) # Edit by Ana: recursive = TRUE -- raw CSVs live in per-subject subfolders (sub-410002/...), same layout the Python pipeline reads. Paths are returned relative to input_path (e.g. "sub-410002/sub-410002_..._e1.csv").
for (i in 1:length(csvSelect)){
  temp_for_file <- ifelse (str_detect(csvSelect[i], "stim_v1_", negate = FALSE), 1, 0)
  if (temp_for_file == 1){
    temp_list <- csvSelect[i]
    datafiles_list <- c(datafiles_list, temp_list)
  }
}

# Edit by Ana: restrict processing to participants deemed FULLY VALID by the Python
# behavioral pipeline (stim-behavioral_extract-trials.py -> stim-behavioral_select-participants.py). Reads the
# inclusion list from stim_valid_behav.csv, so the Python pipeline MUST be run
# first. Stops with an informative error if the file is missing.
valid_behav_file <- "/Users/anagarci/Documents/NDClab/datasets/stim-dataset/derivatives/stim-behavioral/s1_r1/stim_valid_behav.csv"  # Edit by Ana:
if (!file.exists(valid_behav_file)) {  # Edit by Ana:
  stop(paste("Valid-subject list not found:", valid_behav_file, "-- run the Python pipeline (stim-behavioral_extract-trials.py, then stim-behavioral_select-participants.py) first."))  # Edit by Ana:
}  # Edit by Ana:
valid_subs <- read.csv(valid_behav_file)$sub  # Edit by Ana:
file_sub_ids <- as.numeric(gsub("^sub-(\\d+).*", "\\1", datafiles_list))  # Edit by Ana: extract the subject ID from each filename
n_before_valid_filter <- length(datafiles_list)  # Edit by Ana:
datafiles_list <- datafiles_list[file_sub_ids %in% valid_subs]  # Edit by Ana:
print(paste("Restricting to Python-valid subjects:", length(datafiles_list), "of", n_before_valid_filter, "data files kept."))  # Edit by Ana:



# ==============================================================================
# SECTION 4: PER-PARTICIPANT LOOP
# ==============================================================================
# will loop over all participant datafiles.
for(subject in 1:length(datafiles_list)){

  # ----------------------------------------------------------------------------
  # SECTION 4a: SET UP EMPTY OUTPUT FRAMES, LOAD & TRIM THE RAW CSV
  # ----------------------------------------------------------------------------
  # - Rows are appended one at a time with c(...), which coerces EVERYTHING
  #   to character; read the output CSVs back with type conversion (the
  #   stage-2 script's read.csv handles this).
  # - Practice trials are dropped by removing rows where
  #   prac_background_face.started is filled.
  # - newKey is recoded from "right"/"left" to 8/1 so it can be compared
  #   directly with the pressed key codes later (note: the two inline
  #   comments on those lines describe the direction backwards; the code
  #   itself is correct: 'right' -> 8, 'left' -> 1).
  # ----------------------------------------------------------------------------
  # creating an empty data frame that will store all the information drawn from the flanker task of a participant.
  # This data frame will be saved as a csv file for this participant.
  flanker_df <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("participant_id", "current_trial_face", "current_trial_accuracy",
                                                                   "current_trial_congruency", "current_trial_rt", "current_trial_responded",
                                                                   "current_trial_legitResponse", "current_trial_resp_nums"))

  # creating an empty data frame that will store all the information drawn from the surprise task of a participant.
  # This data frame will be saved as a csv file for this participant.
  surprise_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("participant_id", "face", "is_new", "identified_as_new", "memory_surp_rt"))

  #for this participant, find the csv file
  psychopy_file <- paste(input_path,datafiles_list[subject], sep = "/", collapse = NULL)

  #read in the data for this participant, establish id, and remove extraneous variables
  psychopyDat <- read.csv(file = psychopy_file, stringsAsFactors = FALSE, na.strings=c("", "NA"))
  participant_id <- psychopyDat$id[1]

  psychopyDatTrim <- psychopyDat[c("id",
                                   "status", # The displayed face is new? This column stores the correct value of the task not the response from the subject
                                   "newKey", # Which key should be pressed when the face is new!
                                   "congruency",
                                   "accuracy",
                                   "task1_stim_keyResp.keys",
                                   "surprise_key_resp.keys",
                                   "surprise_key_resp.rt",
                                   "prac_background_face.started",
                                   "image_path", # image path for surprise memory task faces
                                   "main_stim_image", # image path for flanker task faces
                                   "task1_stim_keyResp.rt")] #  this stores reaction time for each trial

  #remove practice trials and any rows that do not reflect experiment data
  remove_prac_trials <- subset(psychopyDatTrim, !complete.cases(psychopyDatTrim$prac_background_face.started)) # removes practice trials
  remove_prac_trials$newKey <- replace(remove_prac_trials$newKey, remove_prac_trials$newKey =='right', 8) # replace 8 values with right for the next loop.
  remove_prac_trials$newKey <- replace(remove_prac_trials$newKey, remove_prac_trials$newKey =='left', 1)

  # ----------------------------------------------------------------------------
  # SECTION 4b: FLANKER TRIALS -- PARSE RTs & COUNT RESPONSES
  # ----------------------------------------------------------------------------
  # - Flanker rows = rows with a congruency value (surprise rows have none).
  #   Should be exactly 384 rows per participant.
  # - PsychoPy RT strings like "[0.66, 1.02]" are stripped of brackets and
  #   everything after the first comma -> FIRST response's RT only.
  # - number_of_responses counts all key presses in the trial;
  #   task1_stim_keyResp.keys is then reduced to the FIRST key (1 or 8).
  # ----------------------------------------------------------------------------
  flankerDat <- subset(remove_prac_trials, complete.cases(remove_prac_trials$congruency)) # only keeps flanker trials and removes trials from surprise
  # memory task. # For this study, flankerDat should have only 384 rows!
  flankerDat$task1_stim_keyResp.rt <- str_replace_all(flankerDat$task1_stim_keyResp.rt,"\\[", "") # removes the bracket
  flankerDat$task1_stim_keyResp.rt <- str_replace_all(flankerDat$task1_stim_keyResp.rt,"\\]", "") # removes the bracket
  flankerDat$task1_stim_keyResp.rt <- gsub(",.*","",flankerDat$task1_stim_keyResp.rt) # removing the RT for the second response within the same trial.
  flankerDat$task1_stim_keyResp.rt <- as.numeric(flankerDat$task1_stim_keyResp.rt)


  # Lets add a column that tells how many responses have been made in a given flanker trial
  for (flanker_sub in 1:nrow(flankerDat)){
    response_keys <- str_extract_all(flankerDat$task1_stim_keyResp.keys[flanker_sub], "\\d+\\.?\\d*") # to extract all sequences of digits from the input string
    response_keys <- parse_number(response_keys[[1]]) # to convert the extracted numbers from character strings to numeric values.
    flankerDat$number_of_responses[flanker_sub] <- length(response_keys)
  }
  flankerDat$task1_stim_keyResp.keys <- as.numeric( str_extract(flankerDat$task1_stim_keyResp.keys, '[[:digit:]]')) #extracts the first number (first response) and converts them to numeric

  # ----------------------------------------------------------------------------
  # SECTION 4c: FLANKER TRIALS -- RECODE PER TRIAL & WRITE flankerDat CSV
  # ----------------------------------------------------------------------------
  # - KEY CONVENTION (per lab manual, matches the Python pipeline): when no
  #   response was made, accuracy is 0 -- missing responses count toward
  #   accuracy as errors, but such trials are excluded from RT and other
  #   analyses (their rt is NA).
  # - legitResponse = responded AND rt >= 0.15 s. This is the filter the
  #   stage-2 script uses for the memory hit-rate measures.
  # - num_flanker_trial_removed (printed per participant) counts responded
  #   trials with rt < 0.15 s.
  # ----------------------------------------------------------------------------
  num_flanker_trial_removed <- 0
  # loop over all flanker trials.
  for (trial in 1:nrow(flankerDat)){
    current_trial_face <- flankerDat$main_stim_image[trial]
    current_trial_congruency <- flankerDat$congruency[trial]
    current_trial_rt <- flankerDat$task1_stim_keyResp.rt[trial]
    current_trial_resp_nums <- flankerDat$number_of_responses[trial] # number of responses for the current trial

    if (is.na(flankerDat$task1_stim_keyResp.keys[trial])){ # When no response made in a flanker task trial
      current_trial_responded <- 0 # 0 = not responded; 1 = responded
      # Because of an error in the Python code for the Psychopy task, the accuracy values reported by Psychopy are not correct in trials with no response.
      # Thus, I am putting NAs for accuracy in trials in which no response has been made!
      # Edit by Ana: changed NA -> 0 per the lab manual ("If a response is missing, count
      # them towards accuracy (0)"). The PsychoPy accuracy bug mentioned above was
      # checked empirically on all 59 s1_r1 participants (22,656 trials): every
      # no-response trial already has accuracy = 0 in the raw data, so the bug does
      # not manifest in this dataset and the manual's convention applies.
      current_trial_accuracy <- 0 # Edit by Ana: was NA
    } else if (!is.na(flankerDat$task1_stim_keyResp.keys[trial])){ # When a response made in a flanker task trial
      current_trial_responded <- 1
      current_trial_accuracy <- flankerDat$accuracy[trial]
      if (current_trial_rt < 0.15){
        num_flanker_trial_removed <- num_flanker_trial_removed + 1
      }
    }


    if (current_trial_responded == 1 && current_trial_rt >= 0.15 ){ # Edit by Ana: was > 0.15; manual marks RT < 150 ms invalid, so exactly 150 ms is valid
      current_trial_legitResponse <- 1
    } else {
      current_trial_legitResponse <- 0
    }

    flanker_df[nrow(flanker_df) + 1,] <-c(participant_id, current_trial_face, current_trial_accuracy,
                                          current_trial_congruency, current_trial_rt, current_trial_responded,
                                          current_trial_legitResponse, current_trial_resp_nums)
  } # Closing the loop for each trial

  ### Printing output
  print(paste("Participant ", participant_id, " had ", num_flanker_trial_removed, " of flanker trial removed due to RT."))
  #### end of printing output
  flanker_name <- paste0(participant_id, flanker_csv_fileName, sep = "", collapse = NULL)
  write.csv(flanker_df, paste(output_path, flanker_name, sep = "/", collapse = NULL), row.names=FALSE) # Writing the flanker CSV file to disk

  # ----------------------------------------------------------------------------
  # SECTION 4d: SURPRISE MEMORY TRIALS -- RECODE & WRITE surpriseDat CSV
  # ----------------------------------------------------------------------------
  # - Surprise rows = rows with a newKey value AND a status value; the two
  #   status-NA rows removed correspond to the two block markers.
  # - identified_as_new compares the participant's pressed key against
  #   newKey (the key that MEANS "new" for this participant, counter-
  #   balanced). Equal -> judged New (1); different -> judged Old (0).
  # - !! CAUTION: if a surprise trial has NO response
  #   (surprise_key_resp.keys is NA), the if() comparison evaluates to NA
  #   and R stops with "missing value where TRUE/FALSE needed". The task
  #   presumably forces a response; if the script ever crashes here, that
  #   assumption failed for some participant.
  # - keep_surp_memory_trial_based_on_rt = 1 when rt > 0.2 s; used in
  #   stage 2 both to drop fast-guess trials and to exclude participants
  #   with > 20% such trials. An NA rt would propagate NA into this flag.
  # ----------------------------------------------------------------------------
  ## Creating the Surprise csv file for each participant

  surprise_memory_dat <- subset(remove_prac_trials, complete.cases(remove_prac_trials$newKey)) # keeps rows from the surprise memory task
  surprise_memory_dat <- subset(surprise_memory_dat, complete.cases(surprise_memory_dat$status)) # removes the NA rows  (there should be two NAs as we have two blocks)

  # Looping through surprise memory trials.
  for (surpTrial in 1:nrow(surprise_memory_dat)){
    face <- surprise_memory_dat$image_path[surpTrial]
    is_new <- surprise_memory_dat$status[surpTrial] # 0 = This face is Old (shown during the flanker task); 1 = This face is New (not shown during the flanker task)
    memory_surp_rt <- surprise_memory_dat$surprise_key_resp.rt[surpTrial]

    if (surprise_memory_dat$newKey[surpTrial] == surprise_memory_dat$surprise_key_resp.keys[surpTrial]){
      identified_as_new <- 1 # 1 = Participant has identified this face as New
    } else if(surprise_memory_dat$newKey[surpTrial] != surprise_memory_dat$surprise_key_resp.keys[surpTrial]){
      identified_as_new <- 0 # 0 = Participant has identified this face as Old;
    }



    surprise_df[nrow(surprise_df) + 1,] <-c(participant_id, face, is_new, identified_as_new, memory_surp_rt)

  } # closing the surprise memory trial loop
  # Adding two additional columns to surprise_df to identify whether we should keep that trial in the surprise memory task or not.

  for (kk in 1:nrow(surprise_df)){
    surprise_df$keep_surp_memory_trial_based_on_rt[kk] <- ifelse (surprise_df$memory_surp_rt[kk] > 0.2, 1, 0)
  }

  surprise_name <- paste0(participant_id, surprise_csv_fileName, sep = "", collapse = NULL)
  write.csv(surprise_df, paste(output_path, surprise_name, sep = "/", collapse = NULL), row.names=FALSE) # Writing the surprise CSV file to disk


} # closing the loop for each participant
