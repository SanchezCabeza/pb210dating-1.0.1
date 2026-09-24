#' @title Draws a color-filled ellipse for model validation
#'
#' @description
#' This is an internal routine, not called directly by the user, used to plot 
#' validation ellipses on the age model for visual inspection of 
#' independent chronological markers.
#' 
#' @importFrom graphics polygon
#' @importFrom grDevices colorRampPalette
#'
#' @param x0 x-coordinate of the ellipse center (e.g., estimated year).
#' @param y0 y-coordinate of the ellipse center (e.g., estimated depth).
#' @param a Length of the larger axis.
#' @param b Length of the shorter axis.
#' @param col Color to fill the ellipse. Default is \code{"red"}.
#' @param intensity Color intensity, from 0 (white) to 1 (full color). Default is \code{1}.
#' @param ... Extra arguments passed to \code{\link[graphics]{polygon}} or other plotting functions.
#'
#' @details
#' The non-rotated ellipse equation is used to render a color-filled polygon. 
#' This is particularly useful in the \eqn{^{210}}Pb dating workflow for 
#' representing 1-\eqn{\sigma} and 2-\eqn{\sigma} uncertainties of validation 
#' events (such as \eqn{^{137}}Cs peaks or exotic markers) against the 
#' calculated age model.
#'
#' @returns An ellipse is drawn on top of the active plot (typically the age model). 
#' This helps to visually assess the validity of the geochronology.
#'
#' @references
#' Monteiro, A. (2006). [R] Add ellipse to plot. 
#' \url{https://stat.ethz.ch/pipermail/r-help/2006-October/114652.html}. 
#' Accessed 2025-10-14.
#' 
#' @examples
#' \dontrun{
#'   Intensity <- 1
#'   EllipseFull(x0 = YearRef, y0 = DepthRef, 
#'               a = abs(YearUpper - YearLower), b = abs(DepthLower - DepthUpper),
#'               col = "red", intensity = Intensity / 4)}
#' 
#' @seealso \code{\link{AgeModel}}
#' 
#' @keywords internal
#' @noRd
#' 
EllipseFull <- function(x0, y0, a, b, col = "red", intensity = 1, ...) {
  
  theta <- seq(0, 2 * pi, length = 100)
  x <- x0 + a * cos(theta)
  y <- y0 + b * sin(theta)
  
  # Graded colour
  f.Colour  <- colorRampPalette(c("white", col)) # Colour function
  Colour    <- f.Colour(255)                     # Colours from white to red [1, 255]
  
  # Ensure index is an integer between 1 and 255
  Intensity <- max(1, min(255, round(intensity * 255)))
  
  # Fill the ellipse
  polygon(x, y, col = Colour[Intensity], border = FALSE, ...)
}

#' @title Calculate rolling mean of 2 adjacent values in a vector
#' 
#' @description 
#' Internal helper to calculate the midpoints between adjacent values in a vector.
#' 
#' @param x A numeric vector.
#' @returns A numeric vector of length `length(x) - 1`.
#' 
#' @keywords internal
#' @noRd

RollMean2V <- function(x) {
  if (length(x) < 2) return(numeric(0))
  return((x[-1] + x[-length(x)]) / 2)
}

#' @title Calculate the rolling mean of 2 adjacent columns in a matrix
#' 
#' @description 
#' Internal helper to calculate the midpoints between adjacent columns 
#' of a numeric matrix or data frame.
#' 
#' @param dat A numeric matrix or data frame.
#' @returns A numeric matrix or data frame with \code{ncol(dat) - 1} columns.
#' 
#' @keywords internal
#' @noRd

RollMean2M <- function(dat) {
  n <- ncol(dat)
  if (is.null(n) || n < 2) return(dat)
  return((dat[, 2:n, drop = FALSE] + dat[, 1:(n - 1), drop = FALSE]) / 2)
}

