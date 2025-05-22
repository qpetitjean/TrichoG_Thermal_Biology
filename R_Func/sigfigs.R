#' @title Count Significant Digits in Numeric Values
#' 
#'
#' @description Calculates the number of significant digits in a numeric. This function removes leading and trailing zeros 
#' as well as the decimal point and returns the number of remaining digits.
#'
#'
#' @param x A numeric vector or object coercible to character. Typically a vector of numeric values for which 
#'
#'
#' @return An integer vector indicating the number of significant digits for each input element.
#'
#'
#' @authors Based on Benjamin's solution found on Stack Overflow: https://stackoverflow.com/questions/49730224/finding-the-number-of-significant-digits-versus-digits-after-decimal-in-r
#'
#' @export
#' 

sigfigs <- function(x) {
  orig_scipen <- getOption("scipen")
  options(scipen = 999)
  on.exit(options(scipen = orig_scipen))
  
  x <- as.character(x)
  x <- sub("\\.", "", x)
  x <- gsub("(^0+|0+$)", "", x)
  nchar(x)
}