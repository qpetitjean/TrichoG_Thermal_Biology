#' @title Identify Deviation Points from a Plateau
#'
#' @description 
#' Given a selected plateau interval from smoothed time series, this function identifies the 
#' first time point(s) where the metric significantly deviates from the plateau, based on 
#' non-overlapping 95% confidence intervals. This helps define CTmin/CTmax or recovery 
#' transitions in Cold or Warm ramps.
#'
#' @param SmoothedResCIMetric Data frame. Subset of smoothed and bootstrapped confidence interval data for a single metric.
#' @param PlateEdges Data frame. Two-row data frame containing the beginning and end of the selected plateau.
#' @param Ramp Character. Ramp type: either "Cold" or "Warm".
#' @param frameRate Numeric. Frame rate of the video used in the experiment (default = 25 fps).
#'
#' @return 
#' A named list containing:
#' \describe{
#'   \item{plateauTemp}{Vector of temperature(s) where deviation from plateau is detected. Two values for Cold ramps, one for Warm ramps.}
#'   \item{plateauTime}{Corresponding video frame(s) where deviation occurs.}
#'   \item{plateauMetric}{Smoothed metric value(s) at the deviation point(s).}
#' }
#' 
#' If no significant deviation is found or no confidence interval data is available, 
#' returns the original plateau edge(s).
#'
#' @seealso 
#' \code{\link{extractThermalLimits}}, \code{\link[segmented]{segmented}}
#'
#' @export
#' 
getDeviation <- function(SmoothedResCIMetric, PlateEdges, Ramp, frameRate = 25) {
  # Default: assume no deviation detected
  plateauTemp <- PlateEdges[["Temp"]]
  plateauTime <- PlateEdges[["runTimelinef"]]
  plateauMetric <- PlateEdges[["SmoothedMean"]]
  
  hasCI <- all(c("X97.5.", "X2.5.") %in% names(SmoothedResCIMetric)) &&
    all(!is.na(SmoothedResCIMetric$X97.5.)) &&
    all(!is.na(SmoothedResCIMetric$X2.5.))
  
  if (!hasCI) {
    return(list(
      plateauTemp = plateauTemp,
      plateauTime = plateauTime,
      plateauMetric = plateauMetric
    ))
  }
  
  # Helper to find the deviation from a given edge
  findDeviation <- function(row, side = c("before", "after")) {
    side <- match.arg(side)
    if (side == "before") {
      subset <- SmoothedResCIMetric[which(
        SmoothedResCIMetric$runTimelinef <= row$runTimelinef &
          SmoothedResCIMetric$SmoothedMean < row$SmoothedMean &
          SmoothedResCIMetric$X97.5. < row$X97.5. &
          SmoothedResCIMetric$X2.5. < row$X97.5.
      ), ]
      if (nrow(subset) > 0) {
        return(subset[which.max(subset$runTimelinef), ])
      }
    } else if (side == "after") {
      subset <- SmoothedResCIMetric[which(
        SmoothedResCIMetric$runTimelinef >= row$runTimelinef &
          SmoothedResCIMetric$SmoothedMean < row$SmoothedMean &
          SmoothedResCIMetric$X97.5. < row$X97.5. &
          SmoothedResCIMetric$X2.5. < row$X97.5.
      ), ]
      if (nrow(subset) > 0) {
        return(subset[which.min(subset$runTimelinef), ])
      }
    }
    return(row) # fallback
  }
  
  if (Ramp == "Cold") {
    below1 <- findDeviation(PlateEdges[1, ], "before")
    below2 <- findDeviation(PlateEdges[2, ], "after")
    return(list(
      plateauTemp = c(below1$Temp, below2$Temp),
      plateauTime = c(below1$runTimelinef, below2$runTimelinef),
      plateauMetric = c(below1$SmoothedMean, below2$SmoothedMean)
    ))
  } else if (Ramp == "Warm") {
    below1 <- findDeviation(PlateEdges[1, ], "before")
    return(list(
      plateauTemp = c(below1$Temp, NA),
      plateauTime = c(below1$runTimelinef, NA),
      plateauMetric = c(below1$SmoothedMean, NA)
    ))
  } else {
    # Unknown ramp type, return plate edges
    return(list(
      plateauTemp = plateauTemp,
      plateauTime = plateauTime,
      plateauMetric = plateauMetric
    ))
  }
}