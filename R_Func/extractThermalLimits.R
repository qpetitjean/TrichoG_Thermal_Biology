#' @title Extract Thermal Limits from Smoothed Time Series
#'
#' @description
#' This function identifies thermal limits (e.g., CTmin, CTmax, recovery temperature)
#' for a specified behavioral variable by detecting plateau and deviations in smoothed
#' movement metric over time/temperature ramp. It uses segmented regression
#' to define breakpoints and allows either interactive or non-interactive detection of treshold
#' for reproducibility.
#'
#' In interactive mode, the user defines the number of breakpoints and selects the
#' plateau of interest via the R console. In non-interactive mode, this information
#' is loaded from a previously saved configuration file (.Rds) allowing full reproducibility.
#'
#' @param metricName Character. Name of the metric to analyze (e.g., `"activity"`, `"sinuosity"`).
#' @param folder Character. Name of the folder (assay ID) being processed.
#' @param strainName Character. Strain identifier for the current data.
#' @param Ramp Character. Type of temperature ramp ("Warm" or "Cold").
#' @param Colorvec Character vector. Colors used for plotting each metric.
#' @param metricList Character vector. List of all metric names considered in the analysis.
#' @param SmoothedResCI Data frame. Combined data frame of smoothed means and confidence intervals for all metrics (i.e., output from `temporalBoot()` + `temporalTrend()` merged).
#' @param frameRate Numeric. Frame rate of the video used in the experiment (default = 25 fps).
#' @param GlobSavingDir Character. Path to the root directory where results are saved.
#' @param interactive Logical. Whether to run the function in interactive mode (default = TRUE). If FALSE, a configuration file must be available.
#' @param configFilePath Optional character. If provided, path to a custom configuration `.Rds` file. Otherwise, uses a default path based on `folder`.
#'
#' @return
#' A data frame containing the following extracted metrics:
#' \describe{
#'   \item{CTMax}{Critical thermal maximum temperature (Warm ramps only)}
#'   \item{CTMin}{Critical thermal minimum temperature (Cold ramps only)}
#'   \item{CTMaxTime_frame}{Frame at which CTMax occurred}
#'   \item{CTMinTime_frame}{Frame at which CTMin occurred}
#'   \item{CTMaxValue / CTMinValue}{Metric value at thermal limit}
#'   \item{RecoveryTemp / RecoveryTime / RecoveryValue}{Temperature, time, and value where the metric shows recovery (Cold ramps only)}
#'   \item{meanValue}{Mean value of the metric over the ramp, excluding first 15 min}
#'   \item{maxValue}{Maximum observed metric value before end of ramp}
#'   \item{MetricBased}{Name of the analyzed metric}
#' }
#'
#' Additionally, an `.Rds` file is saved during interactive runs to allow full reproducibility (`Segmentation_Config_[AssayID].Rds`).
#'
#' @details
#' - Plateaus are detected using the `segmented` package on a linear model fitted to the smoothed data.
#' - Users interactively select the plateau of interest (e.g., initial behavioral stability).
#' - Deviations from the plateau are identified based on overlap between the 95% confidence intervals of the smoothed signal and the plateau edges.
#'
#' @seealso
#' \code{\link[segmented]{segmented}}, \code{\link{getDeviation}}
#'
#' @importFrom segmented
#' @export

extractThermalLimits <- function(metricName,
                                 folder,
                                 strainName,
                                 Ramp,
                                 Colorvec,
                                 metricList,
                                 SmoothedResCI,
                                 frameRate = 25,
                                 GlobSavingDir = getwd(),
                                 interactive = TRUE,
                                 configFilePath = NULL) {
  # Prepare file to store or read user input (interactive)
  configName <- paste0("Segmentation_Config_", folder, ".Rds")
  configPath <- file.path(GlobSavingDir,
                          folder,
                          paste("FindMetrics", folder, sep = "_"),
                          configName)
  if (!interactive && !file.exists(configPath)) {
    stop("Non-interactive mode requested but no configuration file was found.")
  }
  
  # Prepare data
  SmoothedResCIMetric <- SmoothedResCI[SmoothedResCI$customfunc == metricName, ]
  SmoothedResCIMetric <- SmoothedResCIMetric[!duplicated(SmoothedResCIMetric$runTimelinef), ]
  
  # specify the model from which breakpoints should be evaluated
  y <- SmoothedResCIMetric$SmoothedMean[!is.na(SmoothedResCIMetric$SmoothedMean)]
  x <- SmoothedResCIMetric$runTimelinef[!is.na(SmoothedResCIMetric$SmoothedMean)]
  mod <- glm(y ~ x)
  
  # retrieve max and mean metric values
  meanValue <- mean(SmoothedResCIMetric$SmoothedMean[SmoothedResCIMetric$runTimelinef >= 15 * 60 * frameRate &
                                                       SmoothedResCIMetric$runTimelinef <= 110 * 60 * frameRate], na.rm = TRUE)
  
  maxValue <- max(SmoothedResCIMetric$SmoothedMean[SmoothedResCIMetric$runTimelinef <= 110 * 60 * frameRate], na.rm = TRUE)
  
  # iterate to find best segmented model based on expected number of breakpoint
  if (interactive) {
    for (i in 1:1000) {
      par(mfrow = c(1, 1))
      plot(
        SmoothedResCIMetric$SmoothedMean ~ SmoothedResCIMetric$runTimelinef,
        ylim = c(0, max(SmoothedResCIMetric$SmoothedMean, na.rm = TRUE) + max(SmoothedResCIMetric$SmoothedMean, na.rm = TRUE)*20/100),
        main = paste(folder, ":", metricName),
        type = "l",
        col = Colorvec[[which(metricList == metricName)]]
      )
      if (all(c("X2.5.", "X97.5.") %in% names(SmoothedResCIMetric))) {
        valid <- complete.cases(SmoothedResCIMetric[, c("X2.5.", "X97.5.")])
        polygon(
          x = c(
            SmoothedResCIMetric$runTimelinef[valid],
            rev(SmoothedResCIMetric$runTimelinef[valid])
          ),
          y = c(
            SmoothedResCIMetric$X2.5.[valid],
            rev(SmoothedResCIMetric$X97.5.[valid])
          ),
          col = adjustcolor(Colorvec[[which(metricList == metricName)]], alpha = 0.15),
          border = NA
        )
      }
      
      BP <- readline("Expected number of breakpoints: ")
      if (is.na(BP)|grepl("^na$|0", BP, ignore.case = TRUE))
        break
      BP <- as.integer(BP)
      Segmod <- segmented::segmented(mod, seg.Z = ~ x, npsi = BP)
      plot(Segmod,
           add = TRUE,
           col = "firebrick3",
           lwd = 2)
      
      satisfied <- readline("Is the result satisfactory? (Y/N): ")
      if (tolower(satisfied) == "y")
        break
    }
    if (interactive & (is.na(BP)|grepl("^na$|0", BP, ignore.case = TRUE))) {
      BP <- NA
      SelectedPlateau <- NA
    } else{
      breaks <- Segmod$psi[, 2]
      breakpoints <- do.call("rbind", lapply(breaks, function(j) {
        SmoothedResCIMetric[which.min(abs(SmoothedResCIMetric$runTimelinef - round(j))), ]
      }))
      
      Plateaux <- c(list(c(0, breakpoints$runTimelinef[1])), lapply(seq_len(nrow(breakpoints)), function(j) {
        if (j == nrow(breakpoints)) {
          c(breakpoints$runTimelinef[j],
            max(SmoothedResCIMetric$runTimelinef))
        } else {
          c(breakpoints$runTimelinef[j],
            breakpoints$runTimelinef[j + 1])
        }
      }))
      print(Plateaux)
      SelectedPlateau <- as.numeric(readline("Select plateau (list index): "))
    }
    # Save user choices to RDS for reproducibility
    userChoices <- list()
    if (file.exists(configPath)) {
      userChoices <- readRDS(configPath)
    }
    userChoices[[strainName]][[folder]][[metricName]] <- list(nBreakpoints = BP, SelectedPlateau = SelectedPlateau)
    saveRDS(userChoices, configPath)
  } else {
    # Load saved config
    userChoices <- readRDS(configPath)
    BP <- userChoices[[strainName]][[folder]][[metricName]]$nBreakpoints
    SelectedPlateau <- userChoices[[strainName]][[folder]][[metricName]]$SelectedPlateau
  }
  if ((is.na(BP)|grepl("^na$|0", BP, ignore.case = TRUE)) |
      is.na(SelectedPlateau)|grepl("^na$|0", SelectedPlateau, ignore.case = TRUE)) {
    outDf <- data.frame(
      Folder = folder,
      Strain = strainName,
      Ramp = Ramp,
      CTMax = NA,
      CTMin = NA,
      CTMaxTime_frame = NA,
      CTMinTime_frame = NA,
      CTMaxValue = NA,
      CTMinValue = NA,
      RecoveryTemp = NA,
      RecoveryTime = NA,
      RecoveryValue = NA,
      meanValue = meanValue,
      maxValue = maxValue,
      MetricBased = metricName
    )
  } else{
    # Segment with stored number of breakpoints
    Segmod <- segmented::segmented(mod, seg.Z = ~ x, npsi = BP)
    breaks <- Segmod$psi[, 2]
    breakpoints <- do.call("rbind", lapply(breaks, function(j) {
      SmoothedResCIMetric[which.min(abs(SmoothedResCIMetric$runTimelinef - round(j))), ]
    }))
    
    Plateaux <- c(list(c(0, breakpoints$runTimelinef[1])), lapply(seq_len(nrow(breakpoints)), function(j) {
      if (j == nrow(breakpoints)) {
        c(breakpoints$runTimelinef[j],
          max(SmoothedResCIMetric$runTimelinef))
      } else {
        c(breakpoints$runTimelinef[j],
          breakpoints$runTimelinef[j + 1])
      }
    }))
    # Extract plateau limits
    PlateEdges <- breakpoints[breakpoints$runTimelinef %in% Plateaux[[SelectedPlateau]], ]
    
    if (all(c("X2.5.", "X97.5.") %in% names(SmoothedResCIMetric))) {
      deviation <- getDeviation(SmoothedResCIMetric, PlateEdges, Ramp, frameRate)
      plateauTemp <- deviation$plateauTemp
      plateauTime <- deviation$plateauTime
      plateauMetric <- deviation$plateauMetric
    } else {
      plateauTemp <- PlateEdges$Temp
      plateauTime <- PlateEdges$runTimelinef
      plateauMetric <- PlateEdges$SmoothedMean
    }
    # store results into a df
    outDf <- data.frame(
      Folder = folder,
      Strain = strainName,
      Ramp = Ramp,
      CTMax = ifelse(Ramp == "Warm", plateauTemp[1], NA),
      CTMin = ifelse(Ramp == "Cold", plateauTemp[1], NA),
      CTMaxTime_frame = ifelse(Ramp == "Warm", plateauTime[1], NA),
      CTMinTime_frame = ifelse(Ramp == "Cold", plateauTime[1], NA),
      CTMaxValue = ifelse(Ramp == "Warm", plateauMetric[1], NA),
      CTMinValue = ifelse(Ramp == "Cold", plateauMetric[1], NA),
      RecoveryTemp = ifelse(Ramp == "Cold", plateauTemp[2], NA),
      RecoveryTime = ifelse(Ramp == "Cold", plateauTime[2], NA),
      RecoveryValue = ifelse(Ramp == "Cold", plateauMetric[2], NA),
      meanValue = meanValue,
      maxValue = maxValue,
      MetricBased = metricName
    )
  }
  
  return(outDf)
}
