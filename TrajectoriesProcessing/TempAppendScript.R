######################################################################################################################
# Script Title: Align Trichogramma sp. Video Tracking Data with Temperature Ramps
# 
# Author: Quentin PETITJEAN
# Date Created: 18/02/2022
# Last Modified: 07/05/2025
# ============================================================================== 
# Requirements:
# - R version 4.2.3
# - Packages:
#   - MoveR v0.3.1: For reading and analyzing TRex tracking data.
#   - remotes v2.4.2: For installing packages from GitHub.
#   - data.table v1.14.8: For efficient data manipulation and file writing.
#   - lubridate v1.9.2: For handling date and time data.
# ============================================================================== 
# Script Overview:
# This script aligns video tracking data of Trichogramma sp. with temperature data recorded during the experiment.
# The workflow includes:
# 1. Installation and loading of required packages.
# 2. Setting up directories and loading cleaned tracking data and temperature ramp files.
# 3. Aligning the tracking data with temperature ramps by:
#    - Extracting strain and experiment date from file names.
#    - Identifying the appropriate temperature ramp file.
#    - Matching video time with temperature data.
#    - Handling missing temperature data by interpolating between consecutive measured points.
# 4. Generating and saving a plot of the temperature ramp across the experiment duration.
# 5. Merging the tracking data with aligned temperature data and saving the resulting dataset.
# ============================================================================== 
# Usage:
# 1. Update the 'Path' variable with the correct path to the Zenodo data repository.
# 2. Ensure the required cleaned tracking files and temperature ramp CSV files are available in the specified directories.
# 3. Run the script in an R environment.
# 4. The script will output the following files in the Results directory:
#    - Plot: 'Temp_Ramp.svg' - visual representation of temperature variation over time.
#    - CSV file: 'Clean_Temp_[AssayID].csv' - cleaned tracking data with appended temperature information.
# ============================================================================== 


######################################################################################################################
# Install and attach needed packages
######################################################################################################################

if (!require(MoveR)) {
  install.packages("remotes")
  remotes::install_github("qpetitjean/MoveR")
}
if (!require(data.table)) {
  install.packages("data.table")
}
if (!require(lubridate)) {
  install.packages("lubridate")
}

library(MoveR)

######################################################################################################################
# set parameters (only change the Path variable)
######################################################################################################################

# the path to the directory downloaded from zenodo repository
Path <- "W:/Postdoc_INRAE_SAM/Zenodo_Data_Repo"

# set the frame rate
frameRate = 25

# path of the directory containing the temperature ramp files (already processed to include video number over the timeline)
TempRampDir <- file.path(Path, "TemperatureRamp")

# the path where results are stored
GlobSavingDir <- file.path(Path, "Results")

# the list of files to process
toprocess <-
  list.files(GlobSavingDir, full.names = FALSE)

# initialize progress bar
total = length(toprocess)
pb <-
  progress::progress_bar$new(format = "Video processing [:bar] :current/:total (:percent)", total = total)
pb$tick(0)

######################################################################################################################
# Run the analysis (loop through the list of files to process)
######################################################################################################################

for (folder in toprocess) {
  # save a list of variable to keep after each loop turn
  envirList <- ls()
  
  # extract the name of the strain from folder name to retrieve temperature ramp file
  folder2 <-
    gsub("[[:digit:]]{1,4}-[[:digit:]]{1,4}-[[:digit:]]{1,4}",
         "",
         folder)
  folder3 <- gsub("^-", "", folder2)
  strainName <- gsub("[^-]*$", "", folder3)
  strainName <- gsub("-$", "", strainName)
  if (length(grep("-", strainName)) > 0) {
    strainName <- gsub("-.*$", "", strainName)
  }
  
  # path of the temperature ramp file (timeline and temperature over the video)
  possibleError0 <-
    tryCatch(
      TempRampFileTemp <-
        list.files(TempRampDir,
                   pattern = strainName,
                   full.names = TRUE),
      error = function(e)
        e
    )
  if (inherits(possibleError0, "error")) {
    print(
      paste(
        "Failed to append temperature data",
        "Problem when reading",
        folder,
        ", consider appending temperature data manually",
        "error: ",
        possibleError0,
        sep = " "
      )
    )
    rm(list = setdiff(ls(), envirList))
    gc()
    # progress bar
    pb$tick(1)
    next
  } else{
    if (length(TempRampFileTemp) > 2) {
      testRes <- c()
      for (i in seq(length(TempRampFileTemp))) {
        testResTemp <-
          grepl(regmatches(
            folder,
            regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", folder)
          ), TempRampFileTemp[[i]])
        testRes <- c(testRes, testResTemp)
      }
      
      if (!(TRUE %in% unique(testRes))) {
        testRes <- c()
        if (isTRUE(grepl("_", regmatches(
          folder,
          regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", folder)
        )))) {
          for (i in seq(length(TempRampFileTemp))) {
            testResTemp <-
              grepl(gsub("_", "-", regmatches(
                folder,
                regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", folder)
              )), TempRampFileTemp[[i]])
            testRes <- c(testRes, testResTemp)
          }
        } else if (isTRUE(grepl("-", regmatches(
          folder,
          regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", folder)
        )))) {
          for (i in seq(length(TempRampFileTemp))) {
            testResTemp <-
              grepl(gsub("-", "_", regmatches(
                folder,
                regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", folder)
              )), TempRampFileTemp[[i]])
            testRes <- c(testRes, testResTemp)
          }
        }
      }
      TempRampFileTemp <- TempRampFileTemp[which(testRes == TRUE)]
    }
    
    if (length(grep("^.*c$|^.*C$", folder)) > 0 |
        length(grep("^.*chaud$|^.*Chaud$", folder)) > 0) {
       TempFile <- TempRampFileTemp[[which(grepl("c$|C$", TempRampFileTemp))]]
    } else{
      TempFile <- TempRampFileTemp[[which(grepl("f$|F$", TempRampFileTemp))]]
    }
    TemperatureFile <- list.files(TempFile, pattern = "*Video_Time.csv", full.names = TRUE)
    
    ## create a folder named cleaned_Temp to save the results of this script
    if (length(list.dirs(paste(
      GlobSavingDir,
      folder,
      paste("Clean_Temp", folder, sep = "_"),
      sep = "/"
    ))) == 0) {
      dir.create(paste(
        GlobSavingDir,
        folder,
        paste("Clean_Temp", folder, sep = "_"),
        sep = "/"
      ))
    }
    
    # the path of the directory where everything will be saved
    SavingDir <- paste(GlobSavingDir,
                      folder,
                      paste("Clean_Temp", folder, sep = "_"),
                      sep = "/")
    # path of the cleaned tracking results (output of the CleaningScript)
    trackDatPath <- paste(paste(
      GlobSavingDir,
      folder,
      paste('cleaned', folder, sep = "_"),
      paste('cleaned', folder, sep = "_"),
      sep = "/"
    ),
    "csv",
    sep = ".")
    
    ################################################################
    #    align tracking data with real time/temperature ramp       #
    ################################################################
    # load cleaned data file
    possibleError <-
      tryCatch(
        trackDatL <-
          as.data.frame(data.table::fread(
            trackDatPath, sep = ",", dec = "."
          )),
        error = function(e)
          e
      )
    if (inherits(possibleError, "error")) {
      print(
        paste(
          "Failed to append temperature data",
          "Problem when reading",
          folder,
          ", consider appending temperature data manually",
          "error: ",
          possibleError,
          sep = " "
        )
      )
      rm(list = setdiff(ls(), envirList))
      gc()
      # progress bar
      pb$tick(1)
      next
    } else{
      # load time to temperature file
      possibleError2 <-
        tryCatch(
          Temp <-
            as.data.frame(data.table::fread(
              TemperatureFile, sep = ";", dec = "."
            )),
          error = function(e)
            e
        )
      if (inherits(possibleError2, "error")) {
        print(
          paste(
            "Failed to append temperature data",
            "Problem when reading temperature ramp file from ",
            folder,
            ", consider appending temperature data manually",
            "error: ",
            possibleError2,
            sep = " "
          )
        )
        rm(list = setdiff(ls(), envirList))
        gc()
        # progress bar
        pb$tick(1)
        next
      } else{
        Temp$RtimeRound <- lubridate::hms(Temp$Time)
        
        # retrieve starting time of the videos
        Vstart <- Temp[!duplicated(Temp[, c("Video_Numb")]), ]
        if (length(Vstart[Vstart$Video_Numb == 0]) > 0) {
          Vstart <- Vstart[-which(Vstart$Video_Numb == 0), ]
        }
        
        # retrieve the total duration and temperature of the Run in time
        ## expressed in second
        runDurs <-
          lubridate::period_to_seconds(tail(Temp$RtimeRound, 1) - Temp$RtimeRound[1]) + 1
        
        ## expressed in frames
        runDurf <- runDurs * frameRate
        
        ## reconstruct the full timeline
        ### in frame
        runTimelinef <- seq(runDurf)
        
        ### seconds
        runTimelineS <- runTimelinef / frameRate
        
        ### expressed as the elapsed time since the beginning of video acquisition
        Timeline <-
          lubridate::seconds_to_period(lubridate::period_to_seconds(Vstart$RtimeRound[1] + runTimelineS))
        Timeline <- as.data.frame(as.character(round(Timeline, 4)))
        
        ## concatenate in a dataframe including temperature and number of the corresponding video chunk
        runTimeline <- cbind(runTimelinef, runTimelineS, Timeline)
        names(runTimeline)[3] <- "Timeline"
        runTimeline$RtimeS <-
          lubridate::period_to_seconds(lubridate::hms(runTimeline$Timeline))
        runTimeline$RtimeSRound <-
          round(lubridate::period_to_seconds(lubridate::seconds_to_period(runTimeline$RtimeS)))
        Temp$RtimeSRound <-
          lubridate::period_to_seconds(Temp$RtimeRound)
        runTimeline <-
          base::merge(runTimeline, Temp[, c("Measured_Temp_Deg_C", "Video_Numb", "RtimeSRound")], by = "RtimeSRound", all = TRUE)
        
        ## remove the duplicates except NA
        runTimeline <- runTimeline[is.na(runTimeline$runTimelinef) | !duplicated(runTimeline$runTimelinef),]
        
        # sometimes temperature ramp forget 1 second
        # replace NA by averaging the last temperature measured before the "Jump"
        # and the first after the "jump" for each missing sequences
        if (NA %in% runTimeline$Measured_Temp_Deg_C) {
          MissingVal <-
            base::split(which(is.na(
              runTimeline$Measured_Temp_Deg_C
            )), base::cumsum(c(1, base::diff(
              which(is.na(
                runTimeline$Measured_Temp_Deg_C
              ))
            ) != 1)))
          
          for (i in seq(length(MissingVal))) {
            # replace NA in temperature column
            runTimeline$Measured_Temp_Deg_C[min(MissingVal[[i]]):max(MissingVal[[i]])] <-
              mean(as.numeric(
                c(
                  runTimeline$Measured_Temp_Deg_C[min(MissingVal[[i]]) - 1],
                  runTimeline$Measured_Temp_Deg_C[max(MissingVal[[i]]) + 1]
                )
              ), na.rm = TRUE)
            # replace NA in video number column
            runTimeline$Video_Numb[min(MissingVal[[i]]):max(MissingVal[[i]])] <-
              mean(as.numeric(
                c(
                  runTimeline$Video_Numb[min(MissingVal[[i]]) - 1],
                  runTimeline$Video_Numb[max(MissingVal[[i]]) + 1]
                )
              ), na.rm = TRUE)
          }
        }
         # plot Temperature ramp according to Real time
        svg(paste(SavingDir, 'Temp_Ramp.svg', sep = "/"))
        plot(
          runTimeline$Measured_Temp_Deg_C ~ lubridate::as_date(
            lubridate::hms(lubridate::seconds_to_period(runTimeline$RtimeS)),
            origin = lubridate::today()
          ),
          type = "l",
          ylab = "Measured Temperature (?C)",
          xlab = "Timeline (hh:mm)",
          main = paste("Temperature measured within the arena for", folder, sep =" ")
        )
        dev.off()
        
        # append time and temperature to tracking data
        # create the timeline of the video according to frame rate
        timeline <-
          stats::setNames(data.frame(matrix(
            ncol = 3, nrow = length(unique(trackDatL$frame))
          )), c("frame", "videotime", "Rtime"))
        timeline$frame <- sort(unique(trackDatL$frame))
        timeline$videotime <-
          lubridate::seconds_to_period(timeline$frame / frameRate)
        
        # retrieve the real timeline from starting time of the video
        timeline$Rtime <-
          round(lubridate::seconds_to_period(
            lubridate::period_to_seconds(Vstart$RtimeRound[1] + timeline$videotime)
          ), 2)
        timeline$RtimeRound <- trunc(timeline$Rtime)
        
        # concatenate tracking data and the timeline df
        trackDatLfinal <-
          base::merge(trackDatL, timeline, by = "frame", all = TRUE)
        trackDatLfinal$RtimeS <-
          lubridate::period_to_seconds(lubridate::hms(trackDatLfinal$Rtime))

        # append the tracking data to the timeline of the run
        finalDat <-
          base::merge(runTimeline,
                      trackDatLfinal,
                      by = "RtimeS",
                      all.x = TRUE)
        finalDat$trackVidId <-
          paste(finalDat$trackletId, finalDat$Video_Numb, sep = "_")
        finalDat <- finalDat[order(finalDat$runTimelinef),]
        
        # remove rows corresponding to movement before and after the start of the experiment
        finalDat <- finalDat[-c(which(is.na(finalDat$runTimelinef))), ]
        finalDat <-
          finalDat[-c(which(
            finalDat$RtimeS < min(Temp$RtimeSRound, na.rm = T) |
              finalDat$RtimeS > max(Temp$RtimeSRound, na.rm = T)
          )),]
        
        # as removing the rows corresponding to movement before and after the start of the experiment may result in 
        # tracklets shorter than 100 frames (the threshold we have defined during the cleaning process), we should remove these incomplete tracklets
        finalDatTracklets <- convert2Tracklets(finalDat, by = "trackletId")
        finalDatTrackletsCleaned <- finalDatTracklets[which(lapply(finalDatTracklets, nrow) >= 100)]
        finalDat <- convert2List(finalDatTrackletsCleaned)
        
        # save the dataset with temperature ramp appended
        data.table::fwrite(
          finalDat,
          paste(paste(
            SavingDir, paste('Clean_Temp', folder, sep = "_"), sep = "/"
          ), "csv", sep = "."),
          sep = ",",
          dec = ".",
          na = "NA"
        )
      }
    }
  }
  rm(list = setdiff(ls(), envirList))
  gc()
  # progress bar
  pb$tick(1)
}
