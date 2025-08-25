######################################################################################################################
# Script Title: Temporal Smoothing of Movement Metrics in Trichogramma sp. Video Tracking Data
# 
# Author: Quentin PETITJEAN
# Date Created: 15/03/2022
# Last Modified: 22/05/2025
# ============================================================================== 
# Requirements:
# - R version 4.2.3
# - Packages:
#   - MoveR v0.3.1: For reading and analyzing TRex tracking data.
#   - remotes v2.4.2: For installing packages from GitHub.
#   - data.table v1.14.8: For efficient data manipulation and file writing.
#   - foreach v1.5.2: For parallel computation of movement metrics.
#   - doParallel v1.0.17: For managing parallel processes.
# - Source Files:
#   - sigfigs.R: Custom function for determining significant digits in numeric values.
# ============================================================================== 
# Script Overview:
# This script performs temporal smoothing of movement metrics for each assay by:
# 1. Importing activity-classified trajectory data including other previously computed metrics.
# 3. Applying a 90-second sliding window with 5-second steps (i.e., every 125 frames with a window size of 2250 frames).
# 4. Weighting the computation by tracklet length to give more influence to longer, more reliable trajectories.
# 5. Parallelizing the smoothing operation across multiple CPU cores for efficiency.
# 6. Outputting the smoothed metrics with corresponding temperature values.
# 7. Generating SVG plots for each metric over the course of the time/temperature ramp.
# ============================================================================== 
# Usage:
# 1. Update the 'Path' variable with the path to the Zenodo data repository.
# 2. Adjust the number of cores in 'CoresToUse' to match available system resources.
# 3. Ensure each assay folder in 'Results' contains the file Clean_Temp_metrics_activ_[AssayID].csv.
# 4. Run the script to generate:
#    - CSV file: Clean_Temp_metrics_activ_smoothed_[AssayID].csv with smoothed metrics and temperature values.
#    - RDS file: Raw output of the temporalTrend R function from MoveR R package.
#    - SVG plots: One per metric, visualizing trends over time.
# ============================================================================== 

######################################################################################################################
# Install and attach needed packages
######################################################################################################################

if (!require(MoveR)) {
  install.packages("remotes")
  remotes::install_github("qpetitjean/MoveR")
}
if (!require(foreach)) {
  install.packages("foreach")
}
if (!require(data.table)) {
  install.packages("data.table")
}


library(MoveR)
library(foreach)

# import a custom function to Count Significant Digits in Numeric Values
source(
  "https://raw.githubusercontent.com/qpetitjean/TrichoG_Thermal_Biology/main/R_Func/sigfigs.R"
)

########################################################################################################################
# set parameters (only change the Path variable)
# ! CAUTION ! This script use parallel computation be careful when setting the number of core to use
########################################################################################################################

frameRate = 25 # Frame rate of the video

# determine the number of available cores on your machine
nbCores <-
  parallel::detectCores(all.tests = FALSE, logical = TRUE)

# number of cores to use for parallel computation 
# (Should be carefully chosen according to the number of core and RAM available)
CoresToUse = 4

# the path to the directory downloaded from zenodo repository
Path <- "W:/Postdoc_INRAE_SAM/Zenodo_Data_Repo"

# the path where results are stored
GlobSavingDir <- file.path(Path, "Results")

# load csv file containing number of indiv and scaling (manually counted using imageJ)
AvTrack <-
  read.csv2(file.path(Path, "VidScale_IndivDensity.csv"),
            sep = ";")
RunCol <-
  "RUN" # the name of the column containing the Run name e.g., "2022-03-01-PUG029-froid"
ScaleCol <-
  "scale.pixel.cm.1" # the name of the column containing the Run name e.g., "2022-03-01-PUG029-froid"

# the list of files to process
toprocess <-
  list.files(GlobSavingDir, full.names = FALSE)

######################################################################################################################
# Run the analysis (loop through the list of files to process)
######################################################################################################################

for (folder in toprocess) {
  # save a list of variable to keep after each loop turn
  envirList <- ls()
  
  # scaling of the video cm per pixels
  scaling <- 1 / as.numeric(AvTrack[which(AvTrack[[RunCol]] == folder), which(colnames(AvTrack) == ScaleCol)])
  # create a folder named cleaned_folder to save the results of this script
  if (length(list.dirs(paste(
    GlobSavingDir,
    folder,
    paste("Clean_Temp_metrics_activ_smoothed", folder, sep = "_"),
    sep = "/"
  ))) == 0) {
    dir.create(paste(
      GlobSavingDir,
      folder,
      paste("Clean_Temp_metrics_activ_smoothed", folder, sep = "_"),
      sep = "/"
    ))
  }
  # the path of the directory where everything will be saved
  SavingDir = paste(
    GlobSavingDir,
    folder,
    paste("Clean_Temp_metrics_activ_smoothed", folder, sep = "_"),
    sep = "/"
  )
  
  # path of the cleaned tracking results (output of the activitycompScriptBIDIME)
  trackDatPath = paste(paste(
    GlobSavingDir,
    folder,
    paste('Clean_Temp_metrics_activ', folder, sep = "_"),
    paste('Clean_Temp_metrics_activ', folder, sep = "_"),
    sep = "/"
  ),
  "csv",
  sep = ".")
  
  possibleError1 <-
  tryCatch(
     trackDatL <-
        as.data.frame(data.table::fread(
          trackDatPath, sep = ";", dec = "."
        )),
         error = function(e)
      e
     )
  if (inherits(possibleError1, "error")) {
    print(paste(possibleError1,folder, sep = " "))
    rm(list = setdiff(ls(), envirList))
    gc()
     next
  } else{
    rm(possibleError1)
    finalDatTracks <- convert2Tracklets(trackDatL, by = "trackletId")
    
    # create a list of function to compute sliding windows on the video timeline
    customFuncList <- list(
      sinuosity =
        function(x) {
          mean(x$sinuosity, na.rm = T)
        },
      Varangle_active = function(x) {
        ifelse(nrow(x[which(!is.na(x$activity2) &
                              x$activity2 == "1"),]) == 0,
               NA,
               circular::var(x$SlideMeanAngle[which(!is.na(x$activity2) &
                                                      x$activity2 == "1")], na.rm = T))
      },
      Varangle_all = function(x) {
        circular::var(x$SlideMeanAngle, na.rm = T)
      },
      speed_active =
        function(x) {
          ifelse(nrow(x[which(!is.na(x$activity2) &
                                x$activity2 == "1"),]) == 0, NA,
                 mean(x$SlidemeanSpeed[which(!is.na(x$activity2) &
                                               x$activity2 == "1")], na.rm = T))
        },
      speed_all =
        function(x) {
          mean(x$SlidemeanSpeed, na.rm = T)
        },
      activity =
        function(x) {
          if (nrow(x[!is.na(x$activity2),]) == 0) {
            NA
          } else {
            nrow(x[!is.na(x$activity2) &
                     x$activity2 == "1",]) / nrow(x[!is.na(x$activity2),])
          }
        },
      Varspeed_active =
        function(x) {
          ifelse(nrow(x[which(!is.na(x$activity2) &
                                x$activity2 == "1"),]) == 0, NA,
                 var(x$SlidemeanSpeed[which(!is.na(x$activity2) &
                                              x$activity2 == "1")], na.rm = T))
        },
      Varspeed_all =
        function(x) {
          var(x$SlidemeanSpeed, na.rm = T)
        },
      Dist2Edge =
        function(x) {
          mean(x$SlidemeanDist2Edge, na.rm = T)
        },
      TraveledDist =
        function(x) {
          mean(x$SlidemeanTraveledDist, na.rm = T)
        },
      EdgeProp = function(x) {
        nrow(x[!is.na(x$Edge) &
                 x$Edge  == "TRUE",]) / nrow(x[!is.na(x$Edge),])
      },
      ActXSpeedOnSin = function(x) {
        if (nrow(x[!is.na(x$activity2), ]) == 0) {
          NA
        } else {
          (nrow(x[!is.na(x$activity2) &
                    x$activity2 == "1", ]) / nrow(x[!is.na(x$activity2), ])) * 
            mean(x$SlidemeanSpeed, na.rm = T) / mean(x$sinuosity, na.rm = T)
        }
      }
    )
    
    ## run sliding window on the metrics included in customFuncList
    
    # Parallelizing the the computation to make analysis faster
    ## determine time duration of the video in frame
    vidT <- max(unlist(lapply(finalDatTracks, function (w)
      max(w["runTimelinef"], na.rm = T))), na.rm = T)

    # create the cluster for parallel computation (here we use the half of the total resource of the computer : 8 cores)
    myCluster <-
      parallel::makeCluster(CoresToUse, # number of cores to use
                            type = "PSOCK", outfile = paste(paste(
                              GlobSavingDir,
                              folder,
                              paste('Clean_Temp_metrics_activ_smoothed', folder, sep = "_"),
                              "logSmooth.txt",
                              sep = "/"
                            )))
    # Register the cluster
    doParallel::registerDoParallel(myCluster)
    
    # create a foreach loop to repeat the function on given time intervals
    ## specify the sampling value 
    sampling = 5 * 25
    toSampleTemp <-  seq(from = 0,
                         to = vidT,
                         by = sampling)
    toLoop <-
      toSampleTemp[c(1, cumsum(rep(ceiling(
        length(toSampleTemp) / CoresToUse
      ), CoresToUse)))]
    if (is.na(toLoop[length(toLoop)])) {
      toLoop[length(toLoop)] <- toSampleTemp[length(toSampleTemp)]
    }
    
    # import function and dataset needed for the computation
    parallel::clusterExport(
      myCluster,
      c(
        "temporalTrend",
        "customFuncList",
        "finalDatTracks",
        "scaling"
      )
    )

    # run the computation
    output.wtd_ALL <-
      foreach::foreach(i = toLoop[1:length(toLoop) - 1], .combine = 'c') %dopar% temporalTrend(
        trackDat = finalDatTracks,
        timeCol = "runTimelinef",
        customFunc = customFuncList,
        Tinterval = c(i, ifelse(i == toLoop[length(toLoop) - 1],
                                toLoop[length(toLoop)],
                                toLoop[which(toLoop == i) + 1])),
        Tstep = 25 * 90,
        sampling = sampling,
        wtd = TRUE
      )
    
    parallel::stopCluster(myCluster)
    
    # save the raw output of the smoothed function
    saveRDS(output.wtd_ALL, file = paste(paste(
      GlobSavingDir,
      folder,
      paste('Clean_Temp_metrics_activ_smoothed', folder, sep = "_"),
      paste(
        paste(
          'Clean_Temp_metrics_activ_smoothed',
          folder,
          "outputWtdALL",
          sep = "_"
        ),
        ".Rds",
        sep = ""
      ),
      sep = "/"
    )))
    
    # group the output computed with different cores
    SmoothedRes <- list()
    for (n in names(customFuncList)) {
      ResTemp <- output.wtd_ALL[which(names(output.wtd_ALL) == n)]
      SmoothedResTemp <-
        do.call("rbind", ResTemp)
      SmoothedResTemp <- SmoothedResTemp[!duplicated(SmoothedResTemp),]
      SmoothedRes[[n]] <-
        data.frame(sapply(unique(names(
          SmoothedResTemp
        )), function(x)
          unname(unlist(
            SmoothedResTemp[, names(SmoothedResTemp) == x]
          ))))
      names(SmoothedRes[[n]]) <-
        gsub("^.*\\.", "", names(SmoothedRes[[n]]))
    }
    
    # transform the output as a DF
    SmoothedResDf <- do.call("cbind", SmoothedRes)
    # simplify column names
    names(SmoothedResDf) <- gsub("^.*\\.", "", names(SmoothedResDf))
    # remove duplicated column (runtimelinef)
    SmoothedResDf <-
      SmoothedResDf[, !duplicated(colnames(SmoothedResDf))]
    
    # append temperature ramp
    finalDatList <- trackDatL
    finalDatList <- finalDatList[order(finalDatList$runTimelinef), ]
    SmoothedResDf$Temp <-
      as.numeric(finalDatList$Measured_Temp_Deg_C[match(SmoothedResDf$runTimelinef, finalDatList$runTimelinef)])
    rm(trackDatL, finalDatList)
    gc()
    data.table::fwrite(
      SmoothedResDf,
      paste(paste(
        GlobSavingDir,
        folder,
        paste('Clean_Temp_metrics_activ_smoothed', folder, sep = "_"),
        paste(
          paste('Clean_Temp_metrics_activ_smoothed', folder, sep = "_"),
          ".csv",
          sep = ""
        ),
        sep = "/"
      )),
      sep = ";",
      dec = ".",
      na = "NA"
    )
    
    # plot the smoothed results
    par(mfrow = c(1, 1))
    for (p in names(SmoothedResDf)) {
      svg(paste(
        paste(paste(
          GlobSavingDir,
          folder,
          paste('Clean_Temp_metrics_activ_smoothed', folder, sep = "_"),
          sep = "/"
        )),
        paste("SmoothedMetrics", "_", p, ".svg", sep = ""),
        sep = "/"
      ))
      plot(
        NULL,
        ylim = c(round(
          min(SmoothedResDf[[p]][!is.infinite(SmoothedResDf[[p]]) &
                                   !is.na(SmoothedResDf[[p]])], na.rm = T),
          digits = max(sigfigs(SmoothedResDf[[p]][!is.infinite(SmoothedResDf[[p]]) &
                                                    !is.na(SmoothedResDf[[p]])]), na.rm =
                         T)
        ),
        round(
          max(SmoothedResDf[[p]][!is.infinite(SmoothedResDf[[p]]) &
                                   !is.na(SmoothedResDf[[p]])], na.rm = T),
          digits = max(sigfigs(SmoothedResDf[[p]][!is.infinite(SmoothedResDf[[p]]) &
                                                    !is.na(SmoothedResDf[[p]])]), na.rm =
                         T)
        )),
        xlim = c(0, max(
          SmoothedResDf$runTimelinef, na.rm = T
        )),
        main = p
      )
      lines(SmoothedResDf[[p]] ~ SmoothedResDf$runTimelinef , col = "red")
      dev.off()
    }
  }
  rm(list = setdiff(ls(), envirList))
  gc()
  print(folder)
}
