#' @title Correct Bootstrapped Confidence Intervals
#'
#' @description 
#' This function filters and corrects bootstrapped confidence interval (CI) outputs (temporalBoot function from MoveR) based on two criteria:
#' (1) the minimum number of tracklets used to estimate the CI, and (2) an expected valid range for the metric values. 
#' CI values that do not meet these criteria are replaced with NA.
#'
#' @param BootRes A named list of data frames, each corresponding to a movement metric, containing 
#' bootstrapped confidence intervals (`2.5%` and `97.5%`) and a column `nbTracklets` indicating 
#' the number of tracklets used in each time window.
#'
#' @param minTracks An integer specifying the minimum number of tracklets required to retain CI estimates.
#' If the number of tracklets for a given time window is below this threshold, both CI bounds are set to NA.
#'
#' @param expectedRange A named list with one element per metric, where each element is a vector of length two
#' specifying the lower and upper expected limits. These can be either numeric values or functions that
#' dynamically compute limits based on CI values.
#'
#' @return A list of data frames with corrected `2.5%` and `97.5%` CI values. Values that fail the
#' tracklet count or expected range criteria are set to NA.
#'
#' @author Quentin Petitjean
#'
#' @export


correctBootCI <- function(BootRes, minTracks, expectedRange) {
  
  for (metric in names(BootRes)) {
    df <- BootRes[[metric]]
    
    # Replace CI values with NA if nbTracklets < minTracks
    lowTracks <- which(df$nbTracklets < minTracks)
    if (length(lowTracks) > 0) {
      df[lowTracks, c("97.5%", "2.5%")] <- NA
    }
    
    # Correct CI values outside expected range
    rangeLim <- expectedRange[[metric]]
    
    if (!is.null(rangeLim)) {
      minLim <- if (is.function(rangeLim[[1]])) rangeLim[[1]](df$`97.5%`) else rangeLim[[1]]
      maxLim <- if (is.function(rangeLim[[2]])) rangeLim[[2]](df$`2.5%`) else rangeLim[[2]]
      
      df$`97.5%`[which(df$`97.5%` < minLim)] <- NA
      df$`2.5%`[which(df$`2.5%` > maxLim)] <- NA
    }
    
    BootRes[[metric]] <- df
  }
  
  return(BootRes)
}
