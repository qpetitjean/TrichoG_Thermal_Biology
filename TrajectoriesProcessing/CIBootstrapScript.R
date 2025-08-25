######################################################################################################################
# Script Title: Compute Bootstrapped 95% Confidence Intervals for Smoothed Movement Metrics in Trichogramma sp.
# 
# Author: Quentin PETITJEAN
# Date Created: 24/03/2022
# Last Modified: 22/05/2025
# ============================================================================== 
# Requirements:
# - R version 4.2.3
# - Packages:
#   - data.table v1.14.8: For efficient data manipulation and file writing.
#   - foreach v1.5.2: For parallel computation of movement metrics.
#   - doParallel v1.0.17: For managing parallel processes.
# - Source Files:
#   - sigfigs.R: Custom function for determining significant digits in numeric values.
#   - correctBootCI.R: Custom function for correcting bootstrap estimation based on minimal number of tracklets and expected range of the metric.
# ============================================================================== 
# Script Overview:
# This script computes bootstrapped 95% confidence intervals (CI) around smoothed movement metrics using 
# a 90-second sliding window sampled every 5 seconds. The steps include:
# 1. Loading activity-classified movement data.
# 2. Computing weighted (by tracklet length) means for each metric of interest.
# 3. Bootstrapping each time point (500 samples) to generate 95% Studentized confidence intervals.
# 4. Saving raw bootstrap outputs and summarizing results for each movement metric.
# 5. correcting bootstrap estimated 95%CI based on the number of tracklets and expected range of the metrics.
# 6. Plotting smoothed time series with CI "envelopes" and exporting SVG figures.
# 7. Exporting the final confidence interval dataset for downstream use.
# ============================================================================== 
# Usage:
# 1. Set the 'Path' variable to the Zenodo data repository location.
# 2. Adjust 'CoresToUse' based on available computational resources.
# 3. Ensure that each assay folder contains `Clean_Temp_metrics_activ_[AssayID].csv`.
# 4. Run the script to generate:
#    - CSV: `Clean_Temp_metrics_activ_smoothed95CI_[AssayID].csv` with bootstrapped 95% CIs per metric.
#    - RDS: Raw bootstrap outputs (`outputWtdBootsALL.rds`) for reproducibility.
#    - SVG: CI plots per metric with smoothed values and shaded intervals.
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
# import a custom function to correct bootstrap estimate
source(
  "https://raw.githubusercontent.com/qpetitjean/TrichoG_Thermal_Biology/main/R_Func/correctBootCI.R"
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
    paste("Clean_Temp_metrics_activ_smoothed95CI", folder, sep = "_"),
    sep = "/"
  ))) == 0) {
    dir.create(paste(
      GlobSavingDir,
      folder,
      paste("Clean_Temp_metrics_activ_smoothed95CI", folder, sep = "_"),
      sep = "/"
    ))
  }
  # the path of the directory where everything will be saved
  SavingDir = paste(
    GlobSavingDir,
    folder,
    paste("Clean_Temp_metrics_activ_smoothed95CI", folder, sep = "_"),
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
    print(paste(possibleError1, folder, sep = " "))
    rm(list = setdiff(ls(), envirList))
    gc()
    next
  } else{
    rm(possibleError1)
    finalDatTracks <- MoveR::convert2Tracklets(trackDatL, by = "trackletId")
    gc()

    # create a list of function to compute sliding windows on the video timeline
    customFuncList <- list(
      sinuosity =
        function(x) {
          mean(x$sinuosity, na.rm = T)
        },
      Varangle_active = function(x) {
        ifelse(nrow(x[which(!is.na(x$activity2) &
                              x$activity2 == "1"), ]) == 0,
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
                                x$activity2 == "1"), ]) == 0, NA,
                 mean(x$SlidemeanSpeed[which(!is.na(x$activity2) &
                                               x$activity2 == "1")], na.rm = T))
        },
      speed_all =
        function(x) {
          mean(x$SlidemeanSpeed, na.rm = T)
        },
      activity =
        function(x) {
          if (nrow(x[!is.na(x$activity2), ]) == 0) {
            NA
          } else {
            nrow(x[!is.na(x$activity2) &
                     x$activity2 == "1", ]) / nrow(x[!is.na(x$activity2), ])
          }
        },
      Varspeed_active =
        function(x) {
          ifelse(nrow(x[which(!is.na(x$activity2) &
                                x$activity2 == "1"), ]) == 0, NA,
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
                 x$Edge  == "TRUE", ]) / nrow(x[!is.na(x$Edge), ])
      },
      ActXSpeedOnSin = function(x) {
        if (nrow(x[!is.na(x$activity2),]) == 0) {
          NA
        } else {
          (nrow(x[!is.na(x$activity2) &
                    x$activity2 == "1",]) / nrow(x[!is.na(x$activity2),])) *
            mean(x$SlidemeanSpeed, na.rm = T) / mean(x$sinuosity, na.rm = T)
        }
      }
    )
    
    # Parallelizing the the computation to make analysis faster
    ## determine time duration of the video in frame
    vidT <- max(unlist(lapply(finalDatTracks, function (w)
      max(w["runTimelinef"], na.rm = T))), na.rm = T)
    # determine the number of available cores
    nbCores <-
      parallel::detectCores(all.tests = FALSE, logical = TRUE)
    # create the cluster for parallel computation 
    myCluster <-
      parallel::makeCluster(CoresToUse, # number of cores to use
                            type = "PSOCK", outfile = paste(
                              paste(
                                GlobSavingDir,
                                folder,
                                paste('Clean_Temp_metrics_activ_smoothed95CI', folder, sep = "_"),
                                "logSmooth.txt",
                                sep = "/"
                              )
                            ))
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
        "temporalBoot",
        "customFuncList",
        "finalDatTracks",
        "scaling"
      )
    )
    
    # run the computation
    output.wtdBoots_ALL <-
      foreach::foreach(i = toLoop[1:length(toLoop) - 1], .combine = 'c') %dopar% temporalBoot(
        trackDat = finalDatTracks,
        timeCol = "runTimelinef",
        customFunc = customFuncList,
        Tinterval = c(i, ifelse(i == toLoop[length(toLoop) - 1],
                                toLoop[length(toLoop)],
                                toLoop[which(toLoop == i) + 1])),
        Tstep = 25 * 90,
        sampling = sampling,
        bootn = 500,
        wtd = TRUE
      )
    
    parallel::stopCluster(myCluster)

    rm(finalDatTracks, 
       myCluster, 
       vidT,
       toLoop)
    gc()
    
    # save the raw output of the smoothed function
    saveRDS(output.wtdBoots_ALL, file = paste(paste(
      GlobSavingDir,
      folder,
      paste('Clean_Temp_metrics_activ_smoothed95CI', folder, sep = "_"),
      paste(
        paste(
          'Clean_Temp_metrics_activ_smoothed95CI',
          folder,
          "outputWtdBootsALL",
          sep = "_"
        ),
        ".Rds",
        sep = ""
      ),
      sep = "/"
    )))

    # group the output computed with different cores: Bootstudents results
    # (it is also possible to retrieve the bootstrap sampling and raw wtd mean and sd for each sliding windows, see the str of the result list)
    boot.ci <-
      output.wtdBoots_ALL[c(which(names(output.wtdBoots_ALL) == "BootCiStudent"))]
    BootRes <- list()
    for (n in names(customFuncList)) {
      BootResTemp <-
        do.call("rbind", lapply(boot.ci, function(x) {
          x[[n]]
        }))
      BootResTemp <- BootResTemp[!duplicated(BootResTemp),]
      BootRes[[n]] <- BootResTemp
      # append temperature data
      BootRes[[n]]$Temp <-
        as.numeric(trackDatL$Measured_Temp_Deg_C[match(BootRes[[n]]$runTimelinef, trackDatL$runTimelinef)])
    }
    
    rm(output.wtdBoots_ALL,
       trackDatL,
       boot.ci)
    gc()
    
    # CI needs some correction as computation made a low number of tracklet may result into wrong estimation
    # here we use a custom `correctBootCI` function to replace CI values performed on less than 5 tracklets 
    # and CI values falling outside the expected range for the considered metrics (see below)
    
    expectedRange <- list(
      sinuosity = list(0, function(x) { 2 * mean(x, na.rm = TRUE) }),
      Dist2Edge = list(-255, 0),
      EdgeProp = list(0, 1),
      activity = list(0, 1),
      Varangle_active = list(0, function(x) { 9.99 * mean(x, na.rm = TRUE) }),
      Varangle_all = list(0, function(x) { 9.99 * mean(x, na.rm = TRUE) }),
      Varspeed_active = list(0, function(x) { 9.99 * mean(x, na.rm = TRUE) }),
      Varspeed_all = list(0, function(x) { 9.99 * mean(x, na.rm = TRUE) }),
      ActXSpeedOnSin = list(0, function(x) { 9.99 * mean(x, na.rm = TRUE) }),
      speed_all = list(0, function(x) { 9.99 * mean(x, na.rm = TRUE) }),
      speed_active = list(0, function(x) { 9.99 * mean(x, na.rm = TRUE) })
    )
    
    BootResCorr <- correctBootCI(BootRes, minTracks = 5, expectedRange)
  
    # remove frame where NA and Inf are returned to draw the 95%IC  as an "envelope"
    SplittedIC <- list()
    for (p in names(BootResCorr)) {
      splitTemp <-
        data.frame(split(BootResCorr[[p]], is.na(BootResCorr[[p]][1]))["FALSE"])
      if (nrow(splitTemp) > 0) {
        splittedTemp <-
          data.frame(split(splitTemp, is.infinite(splitTemp$`FALSE.2.5.`))["FALSE"])
      } else{
        splittedTemp <- data.frame()
      }
      if (nrow(splittedTemp) > 0) {
        SplittedIC[[p]] <-
          data.frame(split(
            splittedTemp,
            is.infinite(splittedTemp$`FALSE.FALSE.97.5.`)
          )["FALSE"])
      } else{
        SplittedIC[[p]] <-  data.frame()
      }
      if (!nrow(SplittedIC[[p]]) > 0) {
        SplittedIC[[p]] <-
          data.frame(matrix(ncol = ncol(BootResCorr[[p]]), nrow = 1, NA))
      }
      names(SplittedIC[[p]]) <- names(BootResCorr[[p]])
    }
    
    # load smoothed metrics data to display them with the 95% CI enveloppe
    SmoothedResDf <- as.data.frame(data.table::fread(paste(
      paste(
        GlobSavingDir,
        folder,
        paste('Clean_Temp_metrics_activ_smoothed', folder, sep = "_"),
        paste(
          paste('Clean_Temp_metrics_activ_smoothed', folder, sep = "_"),
          ".csv",
          sep = ""
        ),
        sep = "/"
      )
    ), sep = ";", dec = "."))
    
    # plot IC as an "envelope", wtd mean as black dots and previously smoothed data as red line
    par(mfrow = c(1, 1))
    for (p in names(SplittedIC)) {
      svg(paste(
        paste(paste(
          GlobSavingDir,
          folder,
          paste('Clean_Temp_metrics_activ_smoothed95CI', folder, sep = "_"),
          sep = "/"
        )),
        paste("SmoothedMetrics95CI", "_", p, ".svg", sep = ""),
        sep = "/"
      ))
      
      if (nrow(SplittedIC[[p]]) <= 1 && is.na(SplittedIC[[p]][1,])) {
        plot(
          NULL,
          ylim = c(0, 1),
          xlim = c(0, 1) ,
          main = p
        )
        mtext("No IC", 1, -2)
      } else{
        plot(
          NULL,
          ylim = c(round(
            min(c(
              SplittedIC[[p]]$`2.5%`, SplittedIC[[p]]$`97.5%`
            ) , na.rm = T), digits = max(sigfigs(SplittedIC[[p]]$`97.5%`), na.rm =
                                           T)
          ),
          round(
            max(c(
              SplittedIC[[p]]$`2.5%`, SplittedIC[[p]]$`97.5%`
            ), na.rm = T), digits = max(sigfigs(SplittedIC[[p]]$`97.5%`), na.rm =
                                          T)
          )),
          xlim = c(0, max(
            SplittedIC[[p]]$runTimelinef, na.rm = T
          )),
          main = p
        )
        lines(SmoothedResDf[[p]] ~ SmoothedResDf$runTimelinef , col = "red")
        points(
          SplittedIC[[p]]$mean ~ SplittedIC[[p]]$runTimelinef,
          pch = 19,
          cex = 0.1
        )
        lines(SplittedIC[[p]]$`97.5%` ~ SplittedIC[[p]]$runTimelinef)
        lines(SplittedIC[[p]]$`2.5%` ~ SplittedIC[[p]]$runTimelinef)
        polygon(
          x = c(
            SplittedIC[[p]]$runTimelinef,
            rev(SplittedIC[[p]]$runTimelinef)
          ),
          y = c(SplittedIC[[p]]$`2.5%`, rev(SplittedIC[[p]]$`97.5%`)),
          col = rgb(1, 0, 0, 0.1),
          border = NA
          ,
          density = NA
        )
      }
      dev.off()
    }
    # transform the output as a DF
    SplittedICDf <- do.call("rbind", SplittedIC)
    SplittedIClength <-
      lapply(names(SplittedIC), function(x)
        nrow(SplittedIC[[x]]))
    names(SplittedIClength) <- names(SplittedIC)
    SplittedICDf <-
      data.frame(SplittedICDf, customfunc = unlist(lapply(names(SplittedIClength), function(x)
        rep(x, SplittedIClength[[x]]))))
    rownames(SplittedICDf) <- NULL
    
    #Save the dataset as csv
    data.table::fwrite(
      SplittedICDf,
      paste(paste(
        GlobSavingDir,
        folder,
        paste('Clean_Temp_metrics_activ_smoothed95CI', folder, sep = "_"),
        paste(
          paste('Clean_Temp_metrics_activ_smoothed95CI', folder, sep = "_"),
          ".csv",
          sep = ""
        ),
        sep = "/"
      )),
      sep = ";",
      dec = ".",
      na = "NA"
    )
    rm(list = setdiff(ls(), envirList))
    gc()
    print(folder)
  }
}
