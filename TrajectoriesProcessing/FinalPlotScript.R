######################################################################################################################
# Script Title: Generate Final Plots of Trichogramma Movement Metrics Over Temperature Ramps
# 
# Author: Quentin PETITJEAN
# Date Created: 24/03/2022
# Last Modified: 23/05/2025
# ============================================================================== 
# Requirements:
# - R version 4.2.3
# - Packages:
#   - MoveR v0.3.1: For movement metrics and sliding window operations.
#   - data.table v1.14.8: For efficient data manipulation and file writing.
# - Source Files:
#   - Strains_plotFunc.r: Custom plotting function for overlaying behavioral and thermal trends.
# ============================================================================== 
# Script Overview:
# This script generates the final plot for each Trichogramma assay. Each plot overlays:
# - the temperature ramp,
# - smoothed behavioral metrics (e.g., activity, speed, sinuosity),
# - 95% confidence intervals (CI),
# - and the relative number of detected individuals over time.
# ============================================================================== 
# Usage:
# 1. Set the 'Path' variable to the base directory containing the Zenodo data repository.
# 2. Ensure the required input files exist in the corresponding assay-specific folders.
# 3. Run the script to generate:
#    - SVG plots: `FinalPlot_[AssayID].svg` with temperature, movement metrics, CI, and detection rate.
#    - Results are saved in `Results/FinalPlot_[AssayID]`.
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

# import a custom function to plot the temporal trends of trichogramma movement over temperature changes
source(
  "https://raw.githubusercontent.com/qpetitjean/TrichoG_Thermal_Biology/main/R_Func/Strains_plotFunc.r"
)

########################################################################################################################
# set parameters (only change the Path variable)
########################################################################################################################

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
IndCol <-
  "manually_counted_Indiv." # the name of the column containing the number of individuals for a given run

# the list of files to process
toprocess <-
  list.files(GlobSavingDir, full.names = FALSE)

# initialize progress bar
total = length(toprocess)
pb <-
  progress::progress_bar$new(format = "Plots processing [:bar] :current/:total (:percent)", total = total)
pb$tick(0)

######################################################################################################################
# Run the analysis (loop through the list of files to process)
######################################################################################################################

for (folder in toprocess) {
  # save a list of variable to keep after each loop turn
  envirList <- ls()
  
  nbrIndcounted = as.numeric(AvTrack[which(AvTrack[[RunCol]] == folder), which(colnames(AvTrack) == IndCol)])
  # create a folder named cleaned_folder to save the results of this script
  if (length(list.dirs(paste(
    GlobSavingDir,
    folder,
    paste("FinalPlot", folder, sep = "_"),
    sep = "/"
  ))) == 0) {
    dir.create(paste(
      GlobSavingDir,
      folder,
      paste("FinalPlot", folder, sep = "_"),
      sep = "/"
    ))
  }
  # the path of the directory where everything will be saved
  SavingDir = paste(GlobSavingDir,
                    folder,
                    paste("FinalPlot", folder, sep = "_"),
                    sep = "/")
  
  # path of the smoothed tracking results (output of the SmoothMetricScriptBIDIME)
  smoothDatPath = paste(paste(
    GlobSavingDir,
    folder,
    paste('Clean_Temp_metrics_activ_smoothed', folder, sep = "_"),
    paste('Clean_Temp_metrics_activ_smoothed', folder, sep = "_"),
    sep = "/"
  ),
  "csv",
  sep = ".")
  # path of the 95%CI tracking results (output of the SmoothMetricScriptBIDIME)
  CIDatPath = paste(paste(
    GlobSavingDir,
    folder,
    paste('Clean_Temp_metrics_activ_smoothed95CI', folder, sep = "_"),
    paste('Clean_Temp_metrics_activ_smoothed95CI', folder, sep = "_"),
    sep = "/"
  ),
  "csv",
  sep = ".")
  
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
  possibleError0 <-
    tryCatch(
      trackDatL <-
        as.data.frame(data.table::fread(
          trackDatPath, sep = ";", dec = "."
        )),
      error = function(e)
        e
    )
  if (inherits(possibleError0, "error")) {
    print(paste(possibleError0, folder, sep = " "))
    rm(list = setdiff(ls(), envirList))
    gc()
    # progress bar
    pb$tick(1)
    next
  } else{
    possibleError1 <-
      tryCatch(
        SmoothedResDf <-
          as.data.frame(data.table::fread(
            smoothDatPath, sep = ";", dec = "."
          )),
        error = function(e)
          e
      )
    if (inherits(possibleError1, "error")) {
     print(paste(possibleError1, folder, sep = " "))
      rm(list = setdiff(ls(), envirList))
      gc()
      # progress bar
      pb$tick(1)
      next
    } else{
      possibleError2 <-
        tryCatch(
          SplittedIC <-
            as.data.frame(data.table::fread(
              CIDatPath, sep = ";", dec = "."
            )),
          error = function(e)
            e
        )
      if (inherits(possibleError2, "error")) {
        print(paste(possibleError2, folder, sep = " "))
        rm(list = setdiff(ls(), envirList))
        gc()
        # progress bar
        pb$tick(1)
        next
      } else{
        rm(
          possibleError0,
          possibleError1,
          possibleError2,
          trackDatPath,
          smoothDatPath,
          CIDatPath
        )
        gc()
        
        # rename the CI column and split the CI data into sublists to work with the strains_plot function
        names(SplittedIC)[[1]] <- "97.5%"
        names(SplittedIC)[[2]] <-"2.5%"
        SplittedICL <- split(SplittedIC, SplittedIC$customfunc)
        
        # save the plot displaying speed, sinuosity and activity
        svg(paste(
          SavingDir,
          paste("FinalPlot", "_", folder, ".svg",  sep = ""),
          sep = "/"
        ),
        width = 9.375,
        height = 5.7291666667)
        strains_plot(
          finalDatList = trackDatL,
          smoothDf = SmoothedResDf,
          TempAxis = "y",
          TempCol = "Temp",
          TimeLimit = c(0, (120*60*25)), 
          TimeCol = "runTimelinef",
          TempScale = if(length(grep("chaud", folder) == 1)){c(18,53,5)} else {c(-8,20,5)},
          colTemp = "#333333",
          nbrIndcounted = nbrIndcounted,
          IC = SplittedICL,
          variables = c("activity" , "speed_active", "sinuosity"),
          colvariable =  c("#339999", "#6600CC", "#669933"),
          scaleVariable = list(c(0, max(SplittedICL[["activity"]]["2.5%"], na.rm=T)), 
                               c(0, max(SplittedICL[["speed_active"]]["2.5%"], na.rm=T)), 
                               c(0, max(SplittedICL[["sinuosity"]]["2.5%"], na.rm=T))),
          Lab = list(
            "Activity (%)",
            expression("Speed (" * list(cm.s ^ -1) * ")"),
            "Sinuosity"
          ),
          xTick = 20,
          InnerSpace = 0.57
        )
        dev.off()
      }
    }
  }
  rm(list = setdiff(ls(), envirList))
  gc()
  # progress bar
  pb$tick(1)
} 
