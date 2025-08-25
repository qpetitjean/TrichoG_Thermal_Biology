######################################################################################################################
# Script Title: Compute Movement Metrics and Activity States from Trichogramma sp. Video Tracking Data
# 
# Author: Quentin PETITJEAN
# Date Created: 11/03/2022
# Last Modified: 07/05/2025
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
#   - FindPnVsFunc.r: Custom function to detect peaks and valleys in a distribution (loaded from GitHub).
# ============================================================================== 
# Script Overview:
# This script computes movement metrics from cleaned Trichogramma sp. video tracking data. 
# It uses parallel processing to compute metrics across large datasets faster. 
# The workflow includes:
# 1. Installation and loading of required packages.
# 2. Setting up directories and loading cleaned tracking data.
# 3. Automatic threshold determination for activity state and edge proximity.
# Note that edge proximity consider the perception distance of trichogramma (about 4 mm, see Wajnberg and Colazza 1998)
# 4. Generating a histogram of speed with the activity threshold.
# 5. Calculation of movement metrics, including:
#    - Sinuosity, turning angle, speed, distance traveled, and activity state.
#    - Distance to arena edge to determine edge proximity.
# 6. Smoothing of computed metrics using sliding windows.
# 7. Saving cleaned data.
# ============================================================================== 
# Usage:
# 1. Update the 'Path' variable with the correct path to the Zenodo data repository.
# 2. Set the number of cores for parallel processing ('CoresToUse') according to available resources. 4 should be enough.
# 3. Ensure the required cleaned tracking files and necessary scripts are available in the specified directories.
# 4. Run the script in an R environment.
# 5. The script will output the following files in the Results directory:
#    - CSV file: 'Clean_Temp_metrics_[AssayID].csv' - cleaned tracking data with computed metrics.
#    - Plot: 'SpeedactivTresh.svg' - speed distribution histogram with activity threshold.
#    - Log files: 'logMetricCompute.txt', 'logMetricSmooth.txt' - logs of parallel computation processes.
# ============================================================================== 
# References:
# Wajnberg, E., Colazza, S., 1998. Genetic variability in the area searched by a parasitic wasp: analysis from automatic video tracking of the walking path. Journal of Insect Physiology 44, 437–444. https://doi.org/10.1016/S0022-1910(98)00032-8
#
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
if (!require(foreach)) {
  install.packages("foreach")
}
if (!require(doParallel)) {
  install.packages("doParallel")
}

library(MoveR)
library(foreach)

# import a custom function to detect peak and valley in a distribution
source(
  "https://raw.githubusercontent.com/qpetitjean/TrichoG_Thermal_Biology/main/R_Func/FindPnVsFunc.r"
)

########################################################################################################################
# set parameters (only change the Path variable)
# ! CAUTION ! This script use parallel computation be careful when setting the number of core to use
########################################################################################################################

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

# path of the directory containing the arena border files (Grey scale gradient to the edge of the arena performed using imageJ macro: manual_arena_detect.txt available in OtherTools directory)
ArenaDir <- file.path(Path, "ArenaBorder")

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
  scaling <-
    1 / as.numeric(AvTrack[which(AvTrack[[RunCol]] == folder), which(colnames(AvTrack) == ScaleCol)]) # scaling of the video cm per pixels
  
  # path of the distance to the arena border matrix
  ArenaFile = list.files(paste(ArenaDir,
                               folder,
                               sep = "/"),
                         pattern = "*.txt",
                         full.names = TRUE)
  # create a folder named cleaned_Temp_metrics to save the results of this script
  if (length(list.dirs(paste(
    GlobSavingDir,
    folder,
    paste("Clean_Temp_metrics", folder, sep = "_"),
    sep = "/"
  ))) == 0) {
    dir.create(paste(
      GlobSavingDir,
      folder,
      paste("Clean_Temp_metrics", folder, sep = "_"),
      sep = "/"
    ))
  }
  # the path of the directory where everything will be saved
  SavingDir = paste(GlobSavingDir,
                    folder,
                    paste("Clean_Temp_metrics", folder, sep = "_"),
                    sep = "/")
  
  # path of the cleaned tracking results (output of the TempAppendScript)
  trackDatPath = paste(paste(
    GlobSavingDir,
    folder,
    paste('Clean_Temp', folder, sep = "_"),
    paste('Clean_Temp', folder, sep = "_"),
    sep = "/"
  ),
  "csv",
  sep = ".")
  
  ##########################
  #  Compute some metrics  #
  ##########################
  
  possibleError <-
    tryCatch(
      trackDatL <- as.data.frame(data.table::fread(
        trackDatPath,
        sep = ",",
        dec = "."
      )),
      error = function(e)
        e
    )
  if (inherits(possibleError, "error")) {
    print(
      paste(
        "Failed to compute metrics",
        "Problem when reading",
        folder,
        ", consider compute metrics manually",
        "error: ",
        possibleError,
        sep = " "
      )
    )
    rm(list = setdiff(ls(), envirList))
    gc()
    next
  } else{
    if (names(trackDatL[1]) == "V1") {
      trackDatL <- trackDatL[, -1]
    }
    if (length(which(names(trackDatL) == "X")) > 0) {
      trackDatL <- trackDatL[, -which(names(trackDatL) == "X")]
    }
    trackDat <- convert2Tracklets(trackDatL, by = "trackletId")
    
    # define parameter which will be used by the customFunc
    ## find the speed threshold above which individuals will be considered as active (TRUE) or not (False)
    activTresh <-
      FindPnVs(
        log10(trackDatL$speed),
        levels = 30,
        n = 3,
        minFreq = 500
      )
    
    if (length(activTresh$valleys) == 0) {
      activTresh$valleys <- rep(min(trackDatL$speed), 2)
    }
    if (length(activTresh$valleys) >= 2) {
      activTresh$valleys <- activTresh$valleys[2]
    }
    
    #### save the plot displaying the threshold
    svg(paste(SavingDir, 'SpeedactivTresh.svg', sep = "/"))
    hist(log10(trackDatL$speed),
         breaks = 100,
         main = "speed (log10)")
    abline(v = activTresh$valleys, col = "red")
    dev.off()
    ## find the distance threshold above which individuals will be considered as at the edge of the arena (TRUE)
    ## or at the center (False)
    ### load gradient to the edge of the arena
    arenaGrad <- as.data.frame(data.table::fread(ArenaFile,
                                                 sep = "\t",
                                                 dec = "."))
    
    ### find the border and draw it
    edge <- data.frame(which(arenaGrad == 1, arr.ind = T))
    names(edge)[c(1, 2)] <- c("y.pos", "x.pos")
    ### find the center and draw it
    center = c(mean(edge[, "x.pos"]), mean(edge[, "y.pos"]))
    TrichPerceptDist <-
      4 # an individual can perceive the arena (about 4 mm, see Wajnberg and Colazza 1998)
    radius <- mean(unlist(sqrt((center[1] - edge["x.pos"]) ^ 2 +
                                 (center[2] - edge["y.pos"]) ^ 2
    )), na.rm = T)
    radiusmm <-
      (round(radius, digits = -2) * scaling) * 10 # convert the arena radius in mm
    edgeTresh <-
      radius * TrichPerceptDist / radiusmm # compute the distance below which individuals perceive the arena edge
    
    # in case trackDat still contains tracklets with length below 1 frame removed them
    if (length(which(unlist(lapply(trackDat, function(x)
      nrow(x))) <= 1)) > 0) {
      trackDat <-
        trackDat[-which(unlist(lapply(trackDat, function(x)
          nrow(x))) <= 1)]
    }
    
    # Specify the batch of function to pass to the analyseTracklets function for metric computation along tracklets
    customFuncList = list(
      # compute sinuosity
      sinuosity = function(x)
        MoveR::sinuosity(x, scale = scaling, timeCol = "runTimelinef"),
      # compute turning angles
      turnAngle = function(x)
        if (nrow(x) >= 3) {
          MoveR::turnAngle(x, unit = "radians", timeCol = "runTimelinef")
        } else{
          NA
        },
      # compute simple activity
      actives = function(x)
        MoveR::activity1(
          x,
          speedCol = "speed",
          minSpeed = 10 ^ activTresh$valleys
        ),
      # compute distance traveled
      distTraveled = function(x)
        MoveR::distTraveled(x, step = 1),
      # compute proportion of individuals at the edge of the arena
      Edge = function(x)
        abs(x$dist2Edge) < edgeTresh
    )
    
    # Parallelizing analyseTracklets to make analysis faster
    ## determine total number of tracklets
    Trackn <- length(trackDat)
    
    # create the cluster for parallel computation (here we use the half of the total resource of the computer : 8 cores)
    myCluster <-
      parallel::makeCluster(
        CoresToUse,
        # number of cores to use
        type = "PSOCK",
        outfile =
          paste(SavingDir,
                "logMetricCompute.txt",
                sep = "/")
      )
    
    # Register the cluster
    doParallel::registerDoParallel(myCluster)
    # create a foreach loop to repeat the function on given time intervals
    toLoop <- seq(from = 1,
                  to = Trackn,
                  by = Trackn / CoresToUse)
    toLoop <- round(toLoop)
    if (!Trackn %in% toLoop) {
      toLoop <- c(toLoop, Trackn)
    }
    # import function and dataset needed for the computation
    parallel::clusterExport(
      myCluster,
      c(
        "analyseTracklets",
        "customFuncList",
        "edge",
        "scaling",
        "edgeTresh",
        "activTresh",
        "trackDat",
        "Trackn"
      )
    )
    
    # run the computation
    trackDat2 <-
      foreach::foreach(i = toLoop[1:length(toLoop) - 1], .combine = 'c') %dopar%
      analyseTracklets(trackDat[i:ifelse(i == toLoop[length(toLoop) - 1],
                                         toLoop[length(toLoop)],
                                         toLoop[which(toLoop == i) + 1] - 1)],
                       customFunc = customFuncList)
    
    parallel::stopCluster(myCluster)
    
    # as parallel processing aggregate lists of results, it do not return a trackletClass object, so specify it
    trackDat2 <- trackletsClass(trackDat2)
    
    ##########################
    #      Smooth metrics    #
    ##########################
    
    # Specify the batch of function to pass to the analyseTracklets function to smooth metrics along tracklets
    customFuncList = list(
      # smooth turning angles
      SlideMeanAngle = function (y)
        MoveR::slidWindow(y$turnAngle,
                      Tstep = 10, 
                      statistic = "mean", 
                      na.rm = T),
      # smooth turning angles variance
      SlideVarAngle = function (y)
        MoveR::slidWindow(y$turnAngle,
                      Tstep = 10, 
                      statistic = "circular.var", 
                      na.rm = T),
      # smooth speed
      SlidemeanSpeed = function (y)
        MoveR::slidWindow(y$speed,
                      Tstep = 10, 
                      statistic = "mean", 
                      na.rm = T),
      # smooth activity
      Slidemeanactivity = function (y)
        MoveR::slidWindow(y$actives,
                      Tstep = 10, 
                      statistic = "mean", 
                      na.rm = T),
      # smooth the proportion of time spent near the edge
      SlidemeanEdgeProp = function (y)
        MoveR::slidWindow(y$Edge,
                      Tstep = 10, 
                      statistic = "mean", 
                      na.rm = T),
      # smooth the distance to the edge
      SlidemeanDist2Edge = function (y)
        MoveR::slidWindow(y$dist2Edge,
                      Tstep = 10, 
                      statistic = "mean", 
                      na.rm = T),
      # smooth speed variance
      SlideVarSpeed = function (y)
        MoveR::slidWindow(y$speed,
                      Tstep = 10, 
                      statistic = "var", 
                      na.rm = T),
      # smooth traveled distance
      SlidemeanTraveledDist = function (y)
        MoveR::slidWindow(y$distTraveled,
                      Tstep = 10, 
                      statistic = "mean", 
                      na.rm = T)
    )
    
    # Parallelizing analyseTracklets to make analysis faster
    ## determine total number of tracklets
    Trackn <- length(trackDat2)
    # create the cluster for parallel computation (here we use the half of the total resource of the computer : 8 cores)
    myCluster <-
      parallel::makeCluster(CoresToUse, # number of cores to use
                            type = "PSOCK",
                            outfile =
                              paste(SavingDir,
                                    "logMetricSmooth.txt",
                                    sep = "/"))
    # Register the cluster
    doParallel::registerDoParallel(myCluster)
    # create a foreach loop to repeat the function on given time intervals
    toLoop <- seq(from = 1,
                  to = Trackn,
                  by = Trackn / CoresToUse)
    toLoop <- round(toLoop)
    if (!Trackn %in% toLoop) {
      toLoop <- c(toLoop, Trackn)
    }
    
    # import function and dataset needed for the computation
    parallel::clusterExport(myCluster,
                            c(
                              "analyseTracklets",
                              "customFuncList",
                              "trackDat2",
                              "Trackn"
                            ))
    
    trackDat3 <-
      foreach::foreach(i = toLoop[1:length(toLoop) - 1], .combine = 'c') %dopar% 
      analyseTracklets(trackDat2[i:ifelse(i == toLoop[length(toLoop) - 1],
                                                toLoop[length(toLoop)],
                                                toLoop[which(toLoop == i) + 1] - 1)],
                             customFunc = customFuncList)
    
    parallel::stopCluster(myCluster)
    
    # save the dataset
    trackDatL <- convert2List(trackletsClass(trackDat3))
    data.table::fwrite(
      trackDatL,
      paste(paste(
        SavingDir,
        paste('Clean_Temp_metrics', folder, sep = "_"),
        sep = "/"
      ), "csv", sep = "."),
      sep = ",",
      dec = ".",
      na = "NA"
    )
  }
  
  print(folder)
  rm(list = setdiff(ls(), envirList))
  gc()
}
