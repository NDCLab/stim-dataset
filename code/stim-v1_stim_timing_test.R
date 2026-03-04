

#######  www.ndclab.com   ###########
# By: Kianoosh Hosseini at NDCLab@FIU
# This script reads text files that include triggers sent to the recorder laptop.
# This script needed to be changed for different timing tests.
# These files have ".vmrk" extension.

# Last update: Mar. 4th, 2026
# Adapted by Ana García Morazzani for the STIM task.

# install.packages("pracma")
library(stringr)
library(pracma)
library(dplyr)

############ THIS SECTION MUST BE UPDATED BEFORE RUNNING THE SCRIPT ############
################################################################################

# Change the address below to where you have stored the .vmrk file.
path <- ("/Users/anagarci/Desktop")
# Change the address below to name of the .vmrk file.
mrkTxt <- readLines(paste(path, "/sys1_2026-03-04_stim_v1_timing-test_flanker-stim-timing.vmrk", sep = "")) 

########################## DO NOT EDIT THE CODE BELOW ##########################
################################################################################

myDat <- setNames(data.frame(matrix(nrow = length(mrkTxt), ncol = 1)), c("colA")) 

for (i in 1:length(mrkTxt)) {
  myDat[i, 1] <- mrkTxt[i]
  
}

newDat <- myDat %>%
  filter(
    str_detect(colA, "Stimulus")
  )

for (i in 1:nrow(newDat)) {
  newDat$colB[i] <- gsub(".*Stimulus,", "", newDat$colA[i]) 
}

for (i in 1:nrow(newDat)) {
  newDat$colC[i] <- gsub(",1,0*.", "", newDat$colB[i]) 
}

newDat <- subset(newDat, select = -c(colA, colB)) 
proc_fileName <- "stim_timing_output.csv"
write.csv(newDat,paste(path,proc_fileName, sep = "/", collapse = NULL), row.names=FALSE)

newDat <- read.csv(paste(path,proc_fileName, sep = "/", collapse = NULL))

flanker_dat <- setNames(data.frame(matrix(ncol = 1)), c("timeDiff"))
dVal <- setNames(data.frame(matrix(ncol = 1)), c("timeDiff"))

for (i in 1:nrow(newDat)) {
  if (str_detect(newDat$colC[i], "S128")) {
    if (i != 1){
      if (str_detect(newDat$colC[i-1], "S  5")  || str_detect(newDat$colC[i-1], "S  6") || str_detect(newDat$colC[i-1], "S  7") || str_detect(newDat$colC[i-1], "S  8")){
        secondVal <- str2num(gsub(".*S128,", "", newDat$colC[i]))
        firstVal <- str2num(gsub(".*,", "", newDat$colC[i-1]))
        diffVal <- secondVal - firstVal
        dVal[1,] <- diffVal
        flanker_dat <- rbind(flanker_dat, dVal)
      }
    } else {
      next
    }
  } 
}


flanker_dat <- na.omit(flanker_dat, na.action = "omit")
time_offset_avg <- mean(flanker_dat$timeDiff)
print(paste("The average time offset is:", time_offset_avg))
time_offset_sd <- sd(flanker_dat$timeDiff)
print(paste("The standard deviation of time offse is:", time_offset_sd))
print(paste("Number of flanker task markers:", nrow(flanker_dat))) # For the stim-v1 task, this value should be 384.

############################## END OF TIMING TEST ##############################
################################################################################



