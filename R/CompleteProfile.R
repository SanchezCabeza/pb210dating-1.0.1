#' @title Complete missing data of core profiles for the CF model
#' 
#' @description
#' The function completes missing values in sediment core profiles used for 
#' \eqn{^{210}}Pb radiochronology.
#' Only interpolation (not extrapolation) is performed. It leaves original data intact. 
#'
#' @param List Dating list created by \code{ReadData} and processed by other package functions.
#' @param verbose Logical.Print progress messages. Default is \code{TRUE}.
#'
#' @details
#' Required surface values: \code{MassSection} and \code{Pb210ExcessDecay}.
#' \code{Ra226} is optional. If no \eqn{^{226}}Ra data is provided, it is skipped with a warning.
#'
#' @details
#' Interpolation rules:
#' \itemize{
#'   \item \code{MassSection}: Linear interpolation vs. mean depth.
#'   \item \code{Ra226}: Linear interpolation vs. cumulative mass.
#'   \item \code{Pb210Total}: Exponential interpolation vs. cumulative mass.
#'   \item \code{Pb210ExcessDecay}: Exponential interpolation vs. cumulative mass.
#' }
#'
#' @returns Updated dating list with complete profiles.
#' 
#' @examples
#' \dontrun{
#' List <- ReadData()
#' List <- CompleteProfile(List)}
#' 
#' @seealso \code{\link{Pb210CF}}. 
#' 
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

CompleteProfile <- function(List, verbose = TRUE) {
  
  if (verbose) cat("Complete profiles for CF modelling. \n")
  Val <- List$Val
  Unc <- List$Unc
  
  #### Extract Base Variables ####
  MassSection      <- Val$MassSection
  MassSectionUnc   <- Unc$MassSection
  DepthLayerBottom <- Val$DepthLayerBottom
  
  #### checks ####
  # Masses
  if (is.na(MassSection[1]))     stop("   - Top section MassSection is missing. Please revise.")
  if (is.na(Unc$MassSection[1])) stop("   - Uncertainty for top MassSection is missing. Please revise.")
  # Check for Pb210ExcessDecay
  if (is.null(List[["Val"]][["Pb210ExcessDecay"]]) || all(is.na(List[["Val"]][["Pb210ExcessDecay"]]))) {
    stop("Variable 'Pb210ExcessDecay' is missing or empty.\n",
      "Please run 'DecayCorrection()' on your dating list before proceeding.",
      call. = FALSE
    )
  }
  
  # Read Ra226
  if (!"Ra226" %in% names(Val) || all(is.na(Val$Ra226))) {
    stop("   - Ra226 is required to compute excess Pb-210. You may need to use ConstantRa().")
  }
  
  # pb_source id the variable used in this script
  pb_source     <- Val$Pb210ExcessDecay
  pb_unc_source <- Unc$Pb210ExcessDecay
  if (is.na(pb_source[1]))     stop("   - Top section of excess Pb210 is missing. Please revise.")
  if (is.na(pb_unc_source[1])) stop("   - Uncertainty for top excess Pb210 is missing. Please revise.")
  
  #### Complete MassSection ####
  DepthLayerExtended <- c(0, DepthLayerBottom)
  DepthSectionMean   <- RollMean2V(DepthLayerExtended)
  
  # approx() with rule = 2 ensures extreme values are clamped to the nearest neighbor (no NA generated)
  MassSection    <- approx(x = DepthSectionMean, y = MassSection, xout = DepthSectionMean, rule = 2)$y
  MassSectionUnc <- approx(x = DepthSectionMean, y = MassSectionUnc, xout = DepthSectionMean, rule = 2)$y
  
  # Update list with completed masses
  Val$MassSection <- MassSection
  Unc$MassSection <- MassSectionUnc
  
  # Cumulative Mass
  MassLayerCumulative   <- c(0, cumsum(MassSection))
  MassSectionCumulative <- RollMean2V(MassLayerCumulative)
  
  # Helper to perform the interpolation
  .interp <- function(x, y, method = "linear") {
    if (method == "exponential") {
      valid_idx <- !is.na(x) & !is.na(y) & y > 0
      # If we have fewer than 2 valid points, we cannot interpolate exponentially
      if (sum(valid_idx) < 2) {
        # Fallback to linear if exponential is mathematically impossible
        return(approx(x = x, y = y, xout = x, rule = 2)$y)
      }
      # loginterpolation 
      log_interp <- approx(x = x[valid_idx], y = log(y[valid_idx]), xout = x, rule = 2)$y
      return(exp(log_interp))
    } else {
      # Standard linear path remains unchanged
      return(approx(x = x, y = y, xout = x, rule = 2)$y)
    }
  }
  
  ## Complete Ra226 (Linear)
  hasRa226 <- !is.null(Val$Ra226) && any(!is.na(Val$Ra226))
  if (hasRa226) {
    Val$Ra226Complete <- .interp(MassSectionCumulative, Val$Ra226, "linear")
    Unc$Ra226Complete <- .interp(MassSectionCumulative, Unc$Ra226, "linear")
  } else {
    if (verbose) cat("   - Warning: No 226Ra data provided. Ra226Complete interpolation will be skipped.\n")
  }
  
  ## Complete Pb210 Excess (Exponential vs. Cumulative Mass)
  Val$Pb210ExcessDecayComplete <- .interp(MassSectionCumulative, pb_source, "exponential")
  Unc$Pb210ExcessDecayComplete <- .interp(MassSectionCumulative, pb_unc_source, "exponential")
  
  # Finalize
  List$Val <- Val[order(names(Val))]
  List$Unc <- Unc[order(names(Unc))]
  if (verbose) cat("   - Complete() successful. \n")  
  invisible(List)
}
