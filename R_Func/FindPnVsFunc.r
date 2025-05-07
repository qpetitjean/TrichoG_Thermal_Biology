#' @title Find peaks and valley in data
#'
#'
#' @description Given a vector of data, the function create the frequency histogram and identify and return
#' local minima and maxima.
#'
#'
#' @param img A numeric vector where to look for minima and maxima
#'
#' @param levels the number of levels from which histogram breaks will be define, high value might introduce noise in
#' minima and maxima detection, on the contrary low value might make difficult to identify them
#'
#' @param n the number of maxima and minima to identify
#'
#' @return this function returns a list containing minima and maxima values and plot the frequency histogram with
#' minima and maxima values as red and blue vertical lines respectively
#'
#'
#' @authors Quentin Petitjean (Based on Evan Friedland solution, see https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima))
#'
#' @export


FindPnVs <- function(img,
                     levels = 100,
                     n = 1, 
                     minFreq = 0) {
  # remove NA and infinite values
  if (length(which(is.na(img))) > 0) {
    img = img[-c(which(is.na(img)))]
  }
  if (length(which(is.infinite(img))) > 0) {
    img = img[-c(which(is.infinite(img)))]
  }
  # draw frequency histogram
  range = c(min(img), max(img))
  breaks = seq(range[1], range[2], length.out = levels + 1)
  h = hist.default(img, breaks = breaks, plot = TRUE)
  h$breaks <- h$breaks[which(h$counts >= minFreq)]
  h$counts <- h$counts[which(h$counts >= minFreq)]
  # find peaks and valley (solution adapted from Evan Friedland: 
  # https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima)
  inflect <- function(x, threshold = 1) {
    up   <- sapply(1:threshold, function(y)
      c(x[-(seq(y))], rep(NA, y)))
    down <-
      sapply(-1:-threshold, function(y)
        c(rep(NA, abs(y)), x[-seq(length(x), length(x) - abs(y) + 1)]))
    a    <- cbind(x, up, down)
    list(minima = which(apply(a, 1, min) == a[, 1]),
         maxima = which(apply(a, 1, max) == a[, 1]))
  }
  bottoms <-
    lapply(seq(n), function(x)
      inflect(h$counts, threshold = x)$minima)
  tops <-
    lapply(seq(n), function(x)
      inflect(h$counts, threshold = x)$maxima)
  # store the results in vector
  valleys <- vector()
  peaks <- vector()
  for (i in seq(n)) {
    bottomsDf <- cbind(bottoms[[i]], h$counts[bottoms[[i]]])
    topsDf <- cbind(tops[[i]], h$counts[tops[[i]]])
    bottom <- bottomsDf[which(bottomsDf[,2] == max(bottomsDf[,2])), 1]
    top <- topsDf[which(topsDf[,2] == max(topsDf[,2])), 1]
    valleys <- c(valleys, h$breaks[bottom])
    peaks <- c(peaks, h$breaks[top])
  }
  # draw a red vertical line at valleys and a blue vertical line at peaks
  abline(v = valleys, col = "red")
  abline(v = peaks, col = "blue")
  return(list(valleys = valleys, peaks = peaks))
}

