######################################################################################################################
# Script Title: Clean and Analyze Trichogramma sp. Video Tracking Data
# 
# Author: Quentin PETITJEAN
# Date Created: 11/03/2022
# Last Modified: 07/05/2022
# ============================================================================== 
# Requirements:
# - R version 4.2.3
# - Packages:
#   - MoveR v0.3.1: For reading and analyzing TRex tracking data.
#   - remotes v2.4.2: For installing packages from GitHub.
#   - data.table v1.14.8: For efficient data manipulation and file writing.
# ============================================================================== 
# Script Overview:
# This script processes video tracking data of Trichogramma sp. by cleaning and filtering raw tracking outputs. 
# The workflow includes:
# 1. Installation and loading of required packages.
# 2. Setting up directories and loading raw data.
# 3. Running the cleaning pipeline to:
#    - Filter out infinite values corresponding to lost individuals.
#    - Filter out particles with lengths outside the 95% CI of individual size.
#    - Calculate and filter individual speeds above the 999th quantile.
#    - Calculate distances to the arena edge and filter trajectories with more than half the body length outside.
# 4. Generating and saving histograms of data characteristics (individual length, speed, distance to edge).
# 5. Saving summary statistics of raw and cleaned data, including:
#    - Filter summary CSV file.
#    - Video summary files before and after filtering.
#    - Cleaned trajectory datasets.
# ============================================================================== 
# Usage:
# 1. Update the 'Path' variable with the correct path to the Zenodo data repository.
# 2. Ensure the required raw tracking files and metadata CSV are available in the specified directories.
# 3. Run the script in an R environment. 
# 4. The script will output the following files in the Results directory:
#    - Histograms: 'Infvalues.svg', 'IndLen.svg', 'IndSpeed.svg', 'IndOut.svg', 'TrackLen.svg'
#    - CSV files: 'FilterSummary.csv', 'VideoSummaryBeforeFilter.csv', 'VideoSummaryAfterFilter.csv', 'cleaned_[AssayID].csv'
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

library(MoveR)

######################################################################################################################
# set parameters (only change the Path variable)
######################################################################################################################

# the path to the directory downloaded from zenodo repository
Path <- "W:/Postdoc_INRAE_SAM/Zenodo_Data_Repo"

# path of the directory containing the files to process (tracking output)
RawDataDir <- file.path(Path, "VideoTracking")

# path of the directory containing the arena border files (Grey scale gradient to the edge of the arena performed using imageJ macro: manual_arena_detect.txt available in OtherTools directory)
ArenaDir <- file.path(Path, "ArenaBorder")

# the list of files to process
toprocess <-
  list.files(RawDataDir, full.names = FALSE)

# load csv file containing number of indiv and scaling (manually counted using imageJ)
AvTrack <-
  read.csv2(file.path(Path, "VidScale_IndivDensity.csv"),
            sep = ";")

RunCol <-
  "RUN" # the name of the column containing the Run name e.g., "2022-03-01-PUG029-froid"
ScaleCol <-
  "scale.pixel.cm.1" # the name of the column containing the Run name e.g., "2022-03-01-PUG029-froid"
CountedIndivs <-
  "manually_counted_Indiv." # the name of the column containing the Run name e.g., "2022-03-01-PUG029-froid"

# the path where each results folders are stored
GlobSavingDir <- file.path(Path, "Results")
if (length(list.dirs(GlobSavingDir)) == 0) {
  dir.create(GlobSavingDir)
}

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
  
  # Set the scaling of the video in cm per pixels
  scaling = 1 / as.numeric(AvTrack[which(AvTrack[[RunCol]] == folder), which(colnames(AvTrack) == ScaleCol)]) 
  
  # Set the manually counted number of individuals
  nbrIndcounted <- as.numeric(AvTrack[which(AvTrack[[RunCol]] == folder), which(colnames(AvTrack) == CountedIndivs)]) 
  
  ## if GlobSavingDir does not contains subfolder named as the strain (folder) create it
  if (length(list.dirs(paste(GlobSavingDir,
                             folder,
                             sep = "/"))) == 0) {
    dir.create(paste(GlobSavingDir,
                     folder,
                     sep = "/"))
  }
  ## create a folder named cleaned_folder to save the results of this script
  if (length(list.dirs(paste(
    GlobSavingDir,
    folder,
    paste("cleaned", folder, sep = "_"),
    sep = "/"
  ))) == 0) {
    dir.create(paste(
      GlobSavingDir,
      folder,
      paste("cleaned", folder, sep = "_"),
      sep = "/"
    ))
  }
  ## the following path are semi-automatically generated
  # path of the tracking results
  videoPath = paste(RawDataDir, folder, sep = "/")
  # path of the distance to the arena border matrix
  ArenaFile = list.files(paste(ArenaDir,
                               folder,
                               sep = "/"),
                         pattern = "*.txt",
                         full.names = TRUE)
  # the path of the directory where everything will be saved
  SavingDir = paste(GlobSavingDir,
                    folder,
                    paste("cleaned", folder, sep = "_"),
                    sep = "/")
  
  if (length(scaling) < 1 |
      length(nbrIndcounted) < 1 |
      length(ArenaFile) < 1 |
      length(list.files(videoPath, full.names = T)) == 0) {
    print(
      paste(
        "Failed to clean:",
        folder,
        "Problem to find scaling, nbrindcounted or arenafile, consider running the cleaning manually",
        sep = " "
      )
    )
    rm(list = setdiff(ls(), envirList))
    gc()
    # progress bar
    pb$tick(1)
    next
  }
  
  ####################
  #  Load TRex Data  #
  ####################
  
  possibleError <-
    tryCatch(
      trackDat <- readTrex(
        trexPath = videoPath,
        imgHeight = 1080,
        flipY = TRUE
      ),
      error = function(e)
        e
    )
  
  if (inherits(possibleError, "error")) {
    print(
      paste(
        "Failed to clean",
        "Problem when reading",
        folder,
        ", consider running the cleaning manually",
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
    
    ####################
    #     Cleaning     #
    ####################
    
    ## define filters based on the presence of infinite values
    ### Check the number of infinite values per time unit
    Data_Trex <- convert2List(trackDat)
    svg(paste(SavingDir, 'Infvalues.svg', sep = "/"))
    hist(Data_Trex$frame[which(is.infinite(Data_Trex$x.pos))],
         main = "freq of Infinite value across time (frame)")
    dev.off()
    
    ### define the filter
    filter.Inf <-
      filterFunc(
        trackDat,
        toFilter = "x.pos",
        customFunc = function(x)
          is.infinite(x)
      )
    
    ### filter Infinite values
    trackDat.Infilt <-
      filterTracklets(
        trackDat,
        filter = filter.Inf,
        splitCond = TRUE,
        minDur = 100
      )
    
    ## define filters based on 95% IC of the individuals length (already listed in the raw tracking data)
    trackDat.InfiltList <- convert2List(trackDat.Infilt[[2]])
    
    ### compute IC
    indLength <- log10(trackDat.InfiltList$maj.ax)
    if (length(which(is.infinite(indLength)) > 0)) {
      indLength <- indLength[-c(which(is.infinite(indLength)))]
    }
    if (length(which(is.na(indLength)) > 0)) {
      indLength <- indLength[-c(which(is.na(indLength)))]
    }
    IC <- quantile(indLength, c(0.025, 0.975))
    
    ### save the plot
    svg(paste(SavingDir, 'IndLen.svg', sep = "/"))
    hist(log10(trackDat.InfiltList$maj.ax),
         breaks = 100,
         main = "Indiv length (log10) and 95% IC")
    abline(v = c(IC[1], IC[2]))
    dev.off()
    
    ### create the filter
    filter.length <-
      filterFunc(
        trackDat.Infilt[[2]],
        toFilter = "maj.ax",
        customFunc = function(x)
          x < 10 ^ IC[1] | x > 10 ^ IC[2]
      )
    ### filter individual length
    trackDat.lenfilt <-
      filterTracklets(trackDat.Infilt[[2]],
                  filter.length,
                  splitCond = TRUE,
                  minDur = 100)
    
    ## define filters based on the 999th percentile of the individuals speed (not listed in the raw tracking data, need some computation)
    ### compute speed
    trackDat2 <- trackDat.lenfilt[[2]]
    trackDat2 <-
      analyseTracklets(trackDat2,
                   customFunc = list(
                     speed = function(x)
                       MoveR::speed(
                         x,
                         scale = scaling,
                         timeCol = "frame",
                         frameR = getInfo(trackDat2)$frameR
                       )
                   ))
    
    ### compute quantile 999th
    trackDat.speedfiltList <- convert2List(trackDat2)
    indSpeed <- log10(trackDat.speedfiltList$speed)
    if (length(which(is.infinite(indSpeed)) > 0)) {
      indSpeed <- indSpeed[-c(which(is.infinite(indSpeed)))]
    }
    if (length(which(is.na(indSpeed)) > 0)) {
      indSpeed <- indSpeed[-c(which(is.na(indSpeed)))]
    }
    quant999th <- quantile(indSpeed, c(0.999))
    ### save the plot
    svg(paste(SavingDir, 'IndSpeed.svg', sep = "/"))
    hist(indSpeed, breaks = 100, main = "Indiv speed (log10) and 999th quantile")
    abline(v = quant999th)
    dev.off()
    ### create the filter
    filter.speed <-
      filterFunc(
        trackDat2,
        toFilter = "speed",
        customFunc = function(x)
          x < 0 | x > 10 ^ quant999th
      )
    ### filter individual speed
    trackDat.speedfilt <-
      filterTracklets(trackDat2,
                  filter.speed,
                  splitCond = TRUE,
                  minDur = 100)

    ## define filter based on the presence of individual outside the arena (> to half body length detected outside)
    arenaGrad <- read.delim(ArenaFile)

    ### retrieve the border
    edge <- data.frame(which(arenaGrad == 1, arr.ind = T))
    names(edge)[c(1, 2)] <- c("y.pos", "x.pos")
    ### retrieve the center
    center <- c(mean(edge[, "x.pos"]), mean(edge[, "y.pos"]))
    
    ### compute the distance to the edge
    trackDat2 <- trackDat.speedfilt[[2]]
    trackDat2 <-
      analyseTracklets(trackDat2,
                   customFunc = list(
                     dist2Edge = function(x)
                       MoveR::dist2Edge(x, edge,
                                       customFunc = "CircularArena")
                   ))
    ### retrieve the distance to the edge
    indOut <- convert2List(trackDat2)$dist2Edge
    
    ### retrieve half the body length of individuals
    meanBodyL <-
      mean(convert2List(trackDat2)$maj.ax, na.rm = T) / 2
    
    svg(paste(SavingDir, 'IndOut.svg', sep = "/"))
    hist(indOut, breaks = 100, main = "Distance to the edge, \npositive values correspond to Indivs that are outside the arena")
    abline(v = meanBodyL)
    dev.off()
    
    ### create the filter (here we consider that an individual is truly detected
    #### outside when more than half of the mean body length is out of the arena)
    filter.out <-
      filterFunc(
        trackDat2,
        toFilter = "dist2Edge",
        customFunc = function(x)
          x > meanBodyL
      )
    
    ### filter individual outside the arena
    trackDat.outfilt <-
      filterTracklets(trackDat2,
                  filter.out,
                  splitCond = TRUE,
                  minDur = 100)

    #### check histogram of trajectory length
    svg(paste(SavingDir, 'TrackLen.svg', sep = "/"))
    hist(log10(unlist(
      lapply(trackDat.outfilt[[2]], function(x)
        dim(x)[1])
    )), breaks = 100, main = "Trajectory length (log10)")
    dev.off()
    
    ### create a summary of each filter and save it as .csv
    FilterSummary <- do.call("cbind",
                             list(
                               data.frame(Infilt = unlist(trackDat.Infilt[[1]])),
                               data.frame(lenfilt = unlist(trackDat.lenfilt[[1]])),
                               data.frame(speedfilt = unlist(trackDat.speedfilt[[1]])),
                               data.frame(outfilt = unlist(trackDat.outfilt[[1]]))
                             ))
    FilterSummary$Stats <- rownames(FilterSummary)
    rownames(FilterSummary) <- NULL
    
    data.table::fwrite(
      FilterSummary,
      paste(paste(SavingDir, 'FilterSummary', sep = "/"), "csv", sep = "."),
      sep = ",",
      dec = ".",
      na = "NA"
    )
    
    # look at the video and trajectories summary
    Data_Trex_stats_before_filter <-
      summary(trackDat)
    
    Data_Trex_stats_after_filter <-
      summary(trackDat.outfilt[[2]])
    
    # save filter stats files
    data.table::fwrite(
      as.data.frame(Data_Trex_stats_after_filter),
      paste(
        paste(SavingDir, 'VideoSummaryAfterFilter', sep = "/"),
        "csv",
        sep = "."
      ),
      sep = ",",
      dec = ".",
      na = "NA"
    )
    data.table::fwrite(
      as.data.frame(Data_Trex_stats_before_filter),
      paste(
        paste(SavingDir, 'VideoSummaryBeforeFilter', sep = "/"),
        "csv",
        sep = "."
      ),
      sep = ",",
      dec = ".",
      na = "NA"
    )
    
    # save the cleaned dataset
    data.table::fwrite(
      convert2List(trackDat.outfilt[[2]]),
      paste(paste(
        SavingDir, paste('cleaned', folder, sep = "_"), sep = "/"
      ), "csv", sep = "."),
      sep = ",",
      dec = ".",
      na = "NA"
    )
    
  }
  rm(list = setdiff(ls(), envirList))
  gc()
  # progress bar
  pb$tick(1)
}
