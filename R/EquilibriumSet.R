#' @title Manually set the equilibrium layer for Pb-210 dating
#'
#' @description
#' Manually assigns the deepest layer below which there is no measurable excess \eqn{^{210}}Pb
#' (the "equilibrium layer"). This value is required for the dating models. 
#'
#' Use this function when you prefer to override the automatic proposal or the interactive 
#' selection performed by \code{\link{Equilibrium}}.
#'
#' @param List Dating list created by \code{\link{ReadData}}.
#' @param LayerDating Numeric. Bottom depth (cm) of the equilibrium layer. 
#' It \strong{must} be one of the values present in \code{Val$DepthLayerBottom}.
#' @param verbose Logical. If TRUE, prints progress logs to the console (default \code{TRUE}).
#'
#' @details
#' The function performs a strict validation: \code{LayerDating} must match exactly 
#' one of the section bottom depths. No interpolation or nearest-neighbour matching is done.
#'
#' @returns Updated dating list with the value of \code{LayerDating} stored in 
#'   \code{$Val$LayerDating}.
#' 
#' @examples
#' \dontrun{
#' List <- EquilibriumSet(List, LayerDating = "your_dating_layer_(cm)")}
#' 
#' @seealso \code{\link{ConstantRa}}, \code{\link{Equilibrium}}, \code{\link{MissingInventory}}, 
#' \code{\link{Pb210CF}}, \code{\link{Pb210CFCS}}
#'
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

EquilibriumSet <- function(List, LayerDating, verbose = TRUE) {
  
  cat("Start the  function to assign the layer below which there is no excess Pb-210. \n")
  
  # Safety checks
  if (!"DepthLayerBottom" %in% names(List$Val)) {
    stop("   - List does not contain DepthLayerBottom. Run ReadData() first.")
  }
  
  if (!is.numeric(LayerDating) || length(LayerDating) != 1 || is.na(LayerDating)) {
    stop("   - LayerDating must be a single numeric value (in cm).")
  }
  
  if (!LayerDating %in% List$Val$DepthLayerBottom) {
    stop("   - LayerDating (", LayerDating, " cm) is not a core layer (cm): ",
         paste(sort(List$Val$DepthLayerBottom), collapse = ", "))
  }
  
  # Assign
  List$Val$LayerDating <- LayerDating
  
  # Official ordering (same as all other package functions)
  List$Val <- List$Val[order(names(List$Val))]
  cat("   - LayerDating manually set to", LayerDating, "cm.\n")
  invisible(List)
}
