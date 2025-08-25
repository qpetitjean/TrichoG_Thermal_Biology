######################################################################################################################
# Script Title: Activity State Classification from Trichogramma sp. Video Tracking Data
# 
# Author: Quentin PETITJEAN
# Date Created: 11/03/2022
# Last Modified: 16/05/2025
# ============================================================================== 
# Requirements:
# - R version 4.2.3
# - Packages:
#   - MoveR v0.3.1: For reading and analyzing TRex tracking data.
#   - remotes v2.4.2: For installing packages from GitHub.
#   - data.table v1.14.8: For efficient data manipulation and file writing.
#   - plotly v4.10.1: For interactive 3D visualization.
#   - htmlwidgets v1.6.2: For saving interactive plots as HTML files.
#   - progress v1.2.2: For displaying progress during data processing.
# ============================================================================== 
# Script Overview:
# This script classifies the activity states (active/inactive) of Trichogramma sp. from video tracking data using density based clustering  (Henning 2020; Ester et al., 1996). 
# The workflow includes:
# 1. Installation and loading of required packages.
# 2. Setting up directories and loading paired warm and cold ramps for each assay.
# 3. Merging data from paired ramps to create a single dataset.
# 4. Classifying activity states using a density-based clustering method:
#    - Activity is determined based on the relationship between speed (log10-transformed) and turning angle variance.
#    - The clustering is performed in a 3D space using speed and turning angle variance.
# 5. Visualizing the classified states using an interactive 3D density plot.
# 6. Saving the 3D plot as an HTML file for interactive viewing.
# 7. Splitting and saving the warm and cold datasets with activity classification.
# ============================================================================== 
# Usage:
# 1. Update the 'Path' variable with the correct path to the Zenodo data repository.
# 2. Ensure the 'Pairfiles.csv' file is available in the specified directory.
# 3. Run the script in an R environment.
# 4. The script will output the following files in the Results directory:
#    - HTML file: '3d_activityClustPlot.html' - Interactive 3D density plot of activity states.
#    - CSV files: 'Clean_Temp_metrics_activ_[AssayID].csv' - Warm and cold datasets with activity classification.
# ============================================================================== 
# References:
# Christian, H., (2020). fpc: Flexible Procedures for Clustering. R package version 2.2-9.
# Ester, M., Kriegel H.P., Sander, J., Xu X., (1996). A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise. Institute for Computer Science, University of Munich. Proceedings of 2nd International Conference on Knowledge Discovery and Data Mining (KDD-96)

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
if (!require(plotly)) {
  install.packages("plotly")
}
if (!require(htmlwidgets)) {
  install.packages("htmlwidgets")
}

library(MoveR)

######################################################################################################################
# set parameters
######################################################################################################################

# the path to the directory downloaded from zenodo repository
Path <- "W:/Postdoc_INRAE_SAM/Zenodo_Data_Repo"

# the path where results are stored
GlobSavingDir <- file.path(Path, "Results")

# load the file indicating which ramps (cold and warm) are paired together
PairList <- data.table::fread(file.path(Path, "Pairfiles.csv"),
                              sep = ";",
                              dec = ".")
toprocess <-
  list.files(GlobSavingDir, full.names = FALSE)

# initialize progress bar
total <- length(toprocess)
pb <-
  progress::progress_bar$new(format = "File processing [:bar] :current/:total (:percent)", total = total)
pb$tick(0)

######################################################################################################################
# Run the analysis (loop through the list of files to process)
######################################################################################################################

for (folder in toprocess) {
  # save a list of variable to keep after each loop turn
  envirList <- ls()
  
  # retrieve paired ramps
  if (!is.null(PairList)) {
    if (!length(which(PairList[["Pair1"]] == folder)) == 0) {
      pair1 <- PairList[which(PairList[["Pair1"]] == folder), ][["Pair1"]]
      pair2 <-
        PairList[which(PairList[["Pair1"]] == folder), ][["Pair2"]]
    } else {
      pair1 <- PairList[which(PairList[["Pair2"]] == folder), ][["Pair2"]]
      pair2 <-
        PairList[which(PairList[["Pair2"]] == folder), ][["Pair1"]]
    }
    FileTemp <- list(
      list.files(GlobSavingDir,
                 pattern = pair1,
                 full.names = TRUE),
      list.files(GlobSavingDir,
                 pattern = pair2,
                 full.names = TRUE)
    )
    
  }
  # extract the name of the strain from file
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
  # find pair ramps from the same strain
  if (is.null(PairList)) {
    possibleError0 <-
      tryCatch(
        FileTemp <-
          list.files(GlobSavingDir,
                     pattern = strainName,
                     full.names = TRUE),
        error = function(e)
          e
      )
    if (inherits(possibleError0, "error")) {
      print(paste(folder, "failed" , possibleError0, sep = " "))
      rm(list = setdiff(ls(), envirList))
      gc()
      pb$tick(1)
      next
    }
  }
  # in case there is more than 2 folder corresponding to the same strain (more than 2 ramps has been made)
  if (length(FileTemp) > 2) {
    testRes <- c()
    for (i in seq(length(FileTemp))) {
      testResTemp <-
        grepl(regmatches(folder,
                         regexpr(
                           "[0-9]{4}-[0-9]{2}-[0-9]{2}", folder
                         )), FileTemp[[i]])
      testRes <- c(testRes, testResTemp)
    }
    
    if (!(TRUE %in% unique(testRes))) {
      testRes <- c()
      if (isTRUE(grepl("_", regmatches(
        folder,
        regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", folder)
      )))) {
        for (i in seq(length(FileTemp))) {
          testResTemp <-
            grepl(gsub("_", "-", regmatches(
              folder,
              regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", folder)
            )), FileTemp[[i]])
          testRes <- c(testRes, testResTemp)
        }
      } else if (isTRUE(grepl("-", regmatches(
        folder,
        regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", folder)
      )))) {
        for (i in seq(length(FileTemp))) {
          testResTemp <-
            grepl(gsub("-", "_", regmatches(
              folder,
              regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", folder)
            )), FileTemp[[i]])
          testRes <- c(testRes, testResTemp)
        }
      }
    }
    FileTemp <- FileTemp[which(testRes == TRUE)]
  }
  if (length(FileTemp) != 2) {
    print(paste(folder, "failed" , "There is a more than two pairfile", sep = " "))
    rm(list = setdiff(ls(), envirList))
    gc()
    pb$tick(1)
    next
  }
  FileTempOrigin <- NULL
  if (length(grep("^.*c$|^.*C$", FileTemp[1])) > 0 |
      length(grep("^.*chaud$|^.*Chaud$", FileTemp[1])) > 0 |
      length(grep("^.*f$|^.*F$", FileTemp[1])) > 0 |
      length(grep("^.*froid$|^.*Froid$", FileTemp[1])) > 0) {
    FileTemp[1] <- FileTemp[1]
  } else {
    FileTempOrigin <- FileTemp
    if (length(grep("^.*froid|^.*Froid", FileTemp[1])) > 0 |
        length(grep("^.*chaud|^.*Chaud", FileTemp[1])) > 0) {
      FileTemp[1] <-
        gsub(
          paste(
            paste("(?<=", strainName, "froid).*", sep = "_"),
            paste("(?<=", strainName, "froid).*", sep = "-"),
            paste("(?<=", strainName, "Froid).*", sep = "_"),
            paste("(?<=", strainName, "Froid).*", sep = "-"),
            paste("(?<=", strainName, "chaud).*", sep = "_"),
            paste("(?<=", strainName, "chaud).*", sep = "-"),
            paste("(?<=", strainName, "Chaud).*", sep = "_"),
            paste("(?<=", strainName, "Chaud).*", sep = "-"),
            sep = "|"
          ),
          "",
          FileTemp[1],
          perl = TRUE
        )
    } else{
      FileTemp[1] <- gsub(
        paste(
          paste("(?<=", strainName, "f).*", sep = "_"),
          paste("(?<=", strainName, "f).*", sep = "-"),
          paste("(?<=", strainName, "F).*", sep = "_"),
          paste("(?<=", strainName, "F).*", sep = "-"),
          paste("(?<=", strainName, "c).*", sep = "_"),
          paste("(?<=", strainName, "c).*", sep = "-"),
          paste("(?<=", strainName, "C).*", sep = "_"),
          paste("(?<=", strainName, "C).*", sep = "-"),
          sep = "|"
        ),
        "",
        FileTemp[1],
        perl = TRUE
      )
    }
  }
  
  if (length(grep("^.*c$|^.*C$", FileTemp[2])) > 0 |
      length(grep("^.*chaud$|^.*Chaud$", FileTemp[2])) > 0 |
      length(grep("^.*f$|^.*F$", FileTemp[2])) > 0 |
      length(grep("^.*froid$|^.*Froid$", FileTemp[2])) > 0) {
    FileTemp[2] <- FileTemp[2]
  } else {
    FileTempOrigin <- FileTemp
    if (length(grep("^.*froid|^.*Froid", FileTemp[2])) > 0 |
        length(grep("^.*chaud|^.*Chaud", FileTemp[2])) > 0) {
      FileTemp[2] <-
        gsub(
          paste(
            paste("(?<=", strainName, "froid).*", sep = "_"),
            paste("(?<=", strainName, "froid).*", sep = "-"),
            paste("(?<=", strainName, "Froid).*", sep = "_"),
            paste("(?<=", strainName, "Froid).*", sep = "-"),
            paste("(?<=", strainName, "chaud).*", sep = "_"),
            paste("(?<=", strainName, "chaud).*", sep = "-"),
            paste("(?<=", strainName, "Chaud).*", sep = "_"),
            paste("(?<=", strainName, "Chaud).*", sep = "-"),
            sep = "|"
          ),
          "",
          FileTemp[2],
          perl = TRUE
        )
    } else{
      FileTemp[2] <- gsub(
        paste(
          paste("(?<=", strainName, "f).*", sep = "_"),
          paste("(?<=", strainName, "f).*", sep = "-"),
          paste("(?<=", strainName, "F).*", sep = "_"),
          paste("(?<=", strainName, "F).*", sep = "-"),
          paste("(?<=", strainName, "c).*", sep = "_"),
          paste("(?<=", strainName, "c).*", sep = "-"),
          paste("(?<=", strainName, "C).*", sep = "_"),
          paste("(?<=", strainName, "C).*", sep = "-"),
          sep = "|"
        ),
        "",
        FileTemp[2],
        perl = TRUE
      )
    }
  }
  
  
  if (!is.null(FileTempOrigin)) {
    Newfolder <- sub(".*/", "", FileTemp[grepl(folder, FileTempOrigin)])
  } else{
    Newfolder <- folder
  }
  
  
  if (length(grep("^.*c$|^.*C$", Newfolder)) > 0 |
      length(grep("^.*chaud$|^.*Chaud$", Newfolder)) > 0) {
    if (length(which(grepl("^.*froid$|^.*Froid$", FileTemp))) == 0) {
      pairFile <- FileTemp[[which(grepl("^.*f$|^.*F$", FileTemp))]]
    } else{
      pairFile <-
        FileTemp[[which(grepl("^.*froid$|^.*Froid$", FileTemp))]]
    }
  } else if (length(grep("^.*f$|^.*F$", Newfolder)) > 0 |
             length(grep("^.*froid$|^.*Froid$", Newfolder)) > 0) {
    if (length(which(grepl("^.*chaud$|^.*Chaud$", FileTemp))) == 0) {
      pairFile <- FileTemp[[which(grepl("^.*c$|^.*C$", FileTemp))]]
    } else{
      pairFile <-
        FileTemp[[which(grepl("^.*chaud$|^.*Chaud$", FileTemp))]]
    }
  }
  if (!is.null(FileTempOrigin)) {
    toTest <- FileTempOrigin[grepl(pairFile, FileTempOrigin)]
  } else{
    toTest <- pairFile
  }
  if (length(which(toprocess %in% sub(".*/", "", toTest))) == 0) {
    print(paste(
      folder,
      "failed" ,
      "Pairfile is not found in toprocess list",
      sep = " "
    ))
    rm(list = setdiff(ls(), envirList))
    gc()
    pb$tick(1)
    next
  }
  # in case the current ramp has already been analysed with its twin ramp go to the next iteration
  if (which(toprocess %in% folder) > which(toprocess %in% sub(".*/", "", toTest))) {
    rm(list = setdiff(ls(), envirList))
    gc()
    pb$tick(1)
    next
  }
  # load data of various metrics use for activity determination from both ramps
  pairFile1 <-
    list.files(
      paste(
        GlobSavingDir,
        folder,
        paste("Clean_Temp_metrics", folder, sep = "_"),
        sep = "/"
      ),
      pattern = "*.csv",
      full.names = TRUE
    )
  pairFile2 <-
    list.files(paste(toTest,
                     paste(
                       "Clean_Temp_metrics", sub(".*/", "", toTest), sep = "_"
                     ),
                     sep = "/"),
               pattern = "*.csv",
               full.names = TRUE)
  if (length(pairFile1) == 0 | length(pairFile2) == 0) {
    print(paste(folder, "failed" , "missing or mispelled pairfile", sep = " "))
    rm(list = setdiff(ls(), envirList))
    gc()
    pb$tick(1)
    next
  }
  
  ## create a folder named cleaned_folder to save the results of this script
  if (length(list.dirs(paste(
    GlobSavingDir,
    Newfolder,
    paste("Clean_Temp_metrics_activ", Newfolder, sep = "_"),
    sep = "/"
  ))) == 0) {
    dir.create(paste(
      GlobSavingDir,
      Newfolder,
      paste("Clean_Temp_metrics_activ", Newfolder, sep = "_"),
      sep = "/"
    ))
  }
  if (length(list.dirs(paste(
    GlobSavingDir,
    sub(".*/", "", toTest),
    paste("Clean_Temp_metrics_activ", sub(".*/", "", toTest), sep = "_"),
    sep = "/"
  ))) == 0) {
    dir.create(paste(
      GlobSavingDir,
      sub(".*/", "", toTest),
      paste("Clean_Temp_metrics_activ",
            sub(".*/", "", toTest),
            sep = "_"),
      sep = "/"
    ))
  }
  # the path of the directory where everything will be saved
  SavingDir <- list(paste(
    GlobSavingDir,
    Newfolder,
    paste("Clean_Temp_metrics_activ", Newfolder, sep = "_"),
    sep = "/"
  ),
  paste(
    toTest,
    paste("Clean_Temp_metrics_activ", sub(".*/", "", toTest), sep = "_"),
    sep = "/"
  ))
  
  ## import the warm and cold ramps
  pairFile1 <-
    list.files(
      paste(
        GlobSavingDir,
        folder,
        paste("Clean_Temp_metrics", folder, sep = "_"),
        sep = "/"
      ),
      pattern = "Clean_Temp_metrics",
      full.names = TRUE
    )
  pairFile2 <-
    list.files(paste(toTest,
                     paste(
                       "Clean_Temp_metrics", sub(".*/", "", toTest), sep = "_"
                     ),
                     sep = "/"),
               pattern = "Clean_Temp_metrics",
               full.names = TRUE)
  
  # merge the cold and warm ramps
  trackDatWarm <- as.data.frame(data.table::fread(pairFile1, sep = ",", dec = "."))
  
  trackDatCold <- as.data.frame(data.table::fread(pairFile2, sep = ",", dec = "."))
  
  trackDatCold$ramp <- "cold"
  trackDatWarm$ramp <- "warm"
  
  trackDatAllList <- rbind(trackDatWarm, trackDatCold)
  trackDatAll <- MoveR::convert2Tracklets(trackDatAllList, by = "trackletId")
  
  # use density based clustering (Henning 2020; Ester et al., 1996) on a two dimensions array to classify particles’ between active and inactive states.
  # here as speed decrease and variance of turning angle increase simultaneously under temperature stress it allow identify inactive (0) vs active (1) momentum in a 3d space 
  trackDatAct <-
    MoveR::activity2(
      trackDatAll,
      var1 = "SlidemeanSpeed",
      var2 = "SlideVarAngle",
      var1T = function(x) {
        log10(x + 1)
      },
      nbins = 100,
      eps = 0.15,
      minPts = 5
    )
  
  # Display the 3D plot showing the distribution of active and inactive states within the specified space
  ## compute the count matrix in the 2d space with the same parameter than activity2 function
  trackDatL <- MoveR::convert2List(trackDatAct)
  outP <-
    MoveR::countMat(
      x = log10(trackDatL[["SlidemeanSpeed"]]),
      y = trackDatL[["SlideVarAngle"]],
      groups = trackDatL[["activity2"]],
      nbins = 100,
      output = "matrix"
    )
  
  ## display the interactive 3d plot using plotly
  ### retrieve max and min value of sqrt-transformed count over the groups for plotting
  maxVal <-
    max(unlist(lapply(outP, function(z)
      max(sqrt(
        z
      )))))
  minVal <-
    min(unlist(lapply(outP, function(z)
      min(sqrt(
        z
      )))))
  
  ### draw the plot using plotly R package (3d interactive plot)
  fig <-
    plotly::plot_ly(
      x =  ~ colnames(outP[[1]]),
      y =  ~ rownames(outP[[1]]),
      contours = list(
        z = list(
          show = TRUE,
          start = round(minVal, -2),
          project = list(z = TRUE),
          end = round(maxVal, -2),
          size = max(maxVal) / 10,
          color = "white"
        )
      )
    )
  ### add inactive layer
  fig <- plotly::add_surface(
    p = fig,
    z = sqrt(outP[[2]]),
    opacity = 0.8,
    colorscale = "Hot",
    cmin = min(sqrt(outP[[1]])),
    cmax = max(sqrt(outP[[1]])),
    colorbar = list(title = "inactive\ncounts (sqrt)")
  )
  ### add active layer
  fig <- plotly::add_surface(
    p = fig,
    z = sqrt(outP[[1]]),
    opacity = 1,
    colorscale = list(
      c(0, 0.25, 1),
      c("rgb(20,20,20)", "rgb(58,139,44)", "rgb(234,239,226)")
    ),
    colorbar = list(title = "active\ncounts (sqrt)")
  )
  fig <- plotly::layout(fig,
                        title = '3D density plot of activity states',
                        scene1 = list(
                          xaxis = list(title = "Speed (log10)"),
                          yaxis = list(title = "Angle variance"),
                          zaxis = list(title = "counts (sqrt)")
                        ))
  fig
  
  ### save the plot
  htmlwidgets::saveWidget(fig,
                          file = paste(SavingDir[[1]], "3d_activityClustPlot.html", sep = "/"))
  
  ### copy the files (plots) that are common for both ramp in the second ramp folder result
  file.copy(
    from = list.files(SavingDir[[1]], pattern = "*.html", full.names = TRUE),
    to = SavingDir[[2]],
    overwrite = TRUE,
    recursive = FALSE,
    copy.mode = TRUE
  )
  
  # split the cold and warm dataset including activity2 classification and save them
  trackDatL <- as.data.frame(trackDatL)
  trackDatLWarm <- trackDatL[which(trackDatL$ramp == "warm"), ]
  trackDatLCold <- trackDatL[which(trackDatL$ramp == "cold"), ]
  
  for (i in 1:2) {
    # save the data
    data.table::fwrite(
      if (i == 1) trackDatLWarm else trackDatLCold,
      paste(
        SavingDir[[i]],
        paste(
          'Clean_Temp_metrics_activ_',
          gsub(".*Clean_Temp_metrics_",
               "",
               sub(".*/", "", if(i == 1) pairFile1 else pairFile2)),
          sep = ""
        ),
        sep = "/"
      ),
      sep = ";",
      dec = ".",
      na = "NA"
    )
  }

rm(list = setdiff(ls(), envirList))
gc()
# progress bar
pb$tick(1)
}
