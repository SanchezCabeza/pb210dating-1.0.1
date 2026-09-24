#' @title Radioactive decay correction of Pb-210 and Po-210 to sampling date
#'
#' @description
#' Corrects the \strong{measured (not-decay-corrected) total} activities 
#' (`Pb210Total` for gamma spectrometry and `Po210Total` for alpha spectrometry) 
#' from the counting/plating period back to the core sampling date.
#' 
#' This function \strong{must} be called \strong{immediately after} \code{\link{ReadData}} and
#' \strong{before} \code{\link{CompleteProfile}}, \code{\link{ConstantRa}}, or any dating model.
#' 
#' The function uses the per-section dates supplied in the input file
#' (`DatePlating`, `MeasurementAlphaStart`, `MeasurementAlphaEnd`,
#' `MeasurementGammaStart`, `MeasurementGammaEnd`).
#' If dates are not provided, no decay correction is performed. 
#' If no \eqn{^{226}}Ra is present (or missing where Pb/Po data exist) it interpolates linearly
#' in depth only the needed positions.
#' 
#' \strong{Correction rules for \eqn{^{210}}Pb\eqn{_{ex}}} (the activity that will be used for dating):
#' \itemize{
#' \item \strong{Gamma (\eqn{^{210}}Pb)}: corrected with \eqn{\lambda_{Pb210}} from the sampling date to the mid-point of the gamma measurement interval.
#' \item \strong{Alpha (\eqn{^{210}}Po)}: corrected with \eqn{\lambda_{Pb210}} from the sampling date to the plating date
#' + \eqn{\lambda_{Po210}} from the plating date to the mid-point of the alpha measurement interval
#' (secular equilibrium assumption).
#' }
#'
#' Uncertainties are scaled by the same correction factor(s). The function works whether
#' the core file contains only gamma, only alpha, or both datasets.
#' 
#' @param List Dating list returned by \code{\link{ReadData}}.
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#'
#' @returns Updated dating list with Decay, Excess, and ExcessDecay fields 
#' in \code{$Val} and \code{$Unc}.
#'
#' @examples
#' \dontrun{
#'   a <- ReadData()
#'   a <- DecayCorrection(a)}
#' 
#' @seealso \code{\link{ReadData}}.
#'
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

DecayCorrection <- function(List, verbose = TRUE) {
  
  cat("Decay Correction started:\n")
  
  Val <- List$Val
  Unc <- List$Unc
  
  #### Safety checks ####
  ## Handle Dates: if blank, do not perform decay correction (all dates = sampling date) 
  # use character dates, later properly transformed
  SampleDate <- format(Val$DateSampling, "%Y-%m-%d")
  if (is.null(Val$DateSampling) | is.na(Val$DateSampling)) stop("Sampling date (DateSampling) is missing. Please revise the input file.\n")
  Val$DatePlating[is.na(Val$DatePlating)]                     <- SampleDate
  Val$MeasurementAlphaStart[is.na(Val$MeasurementAlphaStart)] <- SampleDate
  Val$MeasurementAlphaEnd  [is.na(Val$MeasurementAlphaEnd)  ] <- SampleDate
  Val$MeasurementGammaStart[is.na(Val$MeasurementGammaStart)] <- SampleDate
  Val$MeasurementGammaEnd  [is.na(Val$MeasurementGammaEnd)  ] <- SampleDate
  
  ## Check if we have at least one isotope to process
  has_Pb_data <- !is.null(Val$Pb210Total) && any(!is.na(Val$Pb210Total))
  has_Po_data <- !is.null(Val$Po210Total) && any(!is.na(Val$Po210Total))
  
  if (!has_Pb_data && !has_Po_data)
    stop("No Pb-210 or Po-210 data found to correct.\n")
  
  #### Ra-226 Handling ####
  # Check if we have Ra226 data. If not, ConstantRa() is run
  has_Ra <- !is.null(Val$Ra226) && any(!is.na(Val$Ra226))
  
  if (!has_Ra) {
    cat("   - Ra-226 is required to calculate excess activities but is missing. 
        I run ConstantRa(). This value should be carefully revised before proceeding with dating. \n")
    
    List <- ConstantRa(List = List, mode = "auto", verbose = TRUE)
    Val <- List$Val
    Unc <- List$Unc
  }
  
  # Interpolate Ra-226 only where needed
  idx_need_ra <- (!is.na(Val$Pb210Total) | !is.na(Val$Po210Total)) & is.na(Val$Ra226)
  
  if (any(idx_need_ra)) {
    if (verbose) cat("   - Interpolating missing Ra-226 values.\n")
    known <- !is.na(Val$Ra226)
    if (sum(known) < 2) stop("Cannot interpolate Ra-226: at least two known values are required.")
    
    Val$Ra226[idx_need_ra] <- approx(x = Val$DepthLayerBottom[known], y = Val$Ra226[known], 
                                     xout = Val$DepthLayerBottom[idx_need_ra], method = "linear", rule = 2)$y
    Unc$Ra226[idx_need_ra] <- approx(x = Val$DepthLayerBottom[known], y = Unc$Ra226[known], 
                                     xout = Val$DepthLayerBottom[idx_need_ra], method = "linear", rule = 2)$y
  }
  
  t_sample <- as.Date(Val$DateSampling)
  
  if (verbose) {
    cat("   - Sampling date = ", format(t_sample), "\n")
  }
  
  #### Decay constants ####
  Val$HalfLifePb210 <- if (is.null(Val$HalfLifePb210) || is.na(Val$HalfLifePb210)) 22.23 else Val$HalfLifePb210
  Val$HalfLifePo210 <- if (is.null(Val$HalfLifePo210) || is.na(Val$HalfLifePo210)) (138.3763 / 365.25) else Val$HalfLifePo210
  
  lambdaPb210 <- log(2) / Val$HalfLifePb210
  lambdaPo210 <- log(2) / Val$HalfLifePo210
  
  #### 1. Gamma (Pb-210) Correction ####
  if (has_Pb_data) {
    t_start_gamma <- as.Date(Val$MeasurementGammaStart)
    t_end_gamma   <- as.Date(Val$MeasurementGammaEnd)
    t_end_gamma[is.na(t_end_gamma)] <- t_start_gamma[is.na(t_end_gamma)]
    t_meas_gamma  <- t_start_gamma + (t_end_gamma - t_start_gamma) / 2
    delta_gamma   <- as.numeric(difftime(t_meas_gamma, t_sample, units = "days")) / 365.25
    corr_gamma    <- exp(lambdaPb210 * delta_gamma)
    if (verbose) cat("   - Gamma correction factor range =", round(range(corr_gamma, na.rm = TRUE), 4), "\n")
  }
  
  #### 2. Alpha (Po-210) Correction ####
  if (has_Po_data) {
    t_plated      <- as.Date(Val$DatePlating)
    t_start_alpha <- as.Date(Val$MeasurementAlphaStart)
    t_end_alpha   <- as.Date(Val$MeasurementAlphaEnd)
    t_end_alpha[is.na(t_end_alpha)] <- t_start_alpha[is.na(t_end_alpha)]
    t_meas_alpha  <- t_start_alpha + (t_end_alpha - t_start_alpha) / 2
    
    delta_alpha_plating     <- as.numeric(difftime(t_plated, t_sample, units = "days")) / 365.25
    delta_alpha_measurement <- as.numeric(difftime(t_meas_alpha, t_plated, units = "days")) / 365.25
    
    corr_alpha_plating      <- exp(lambdaPb210 * delta_alpha_plating)
    corr_alpha_measurement  <- exp(lambdaPo210 * delta_alpha_measurement)
    corr_alpha              <- corr_alpha_plating * corr_alpha_measurement
    if (verbose) cat("   - Alpha correction factor range =", round(range(corr_alpha, na.rm = TRUE), 4), "\n")
  }
  
  #### 4. Apply Corrections ####
  
  if (has_Po_data) {
    Val$Po210Excess      <- Val$Po210Total - Val$Ra226
    Unc$Po210Excess      <- sqrt(Unc$Po210Total^2 + Unc$Ra226^2)
    Val$Po210ExcessDecay <- Val$Po210Excess * corr_alpha
    Unc$Po210ExcessDecay <- Unc$Po210Excess * corr_alpha
    if (verbose) cat("   - Excess Po210 (alpha) corrected.\n")
  }

  if (has_Pb_data) {
    Val$Pb210Excess      <- Val$Pb210Total - Val$Ra226
    Unc$Pb210Excess      <- sqrt(Unc$Pb210Total^2 + Unc$Ra226^2)
    Val$Pb210ExcessDecay <- Val$Pb210Excess * corr_gamma
    Unc$Pb210ExcessDecay <- Unc$Pb210Excess * corr_gamma
    if (verbose) cat("   - Excess Pb210 (gamma) corrected.\n")
  } else {
    # if no Pb data are present, copy Po210 data 
    Val$Pb210Excess      <- Val$Po210Excess
    Unc$Pb210Excess      <- Unc$Po210Excess
    Val$Pb210ExcessDecay <- Val$Po210ExcessDecay
    Unc$Pb210ExcessDecay <- Unc$Po210ExcessDecay
    if (verbose) cat("   - Po210 data have been copied to Pb210 for dating.\n")
  }
  
  # Finalize
  List$Val <- Val[order(names(Val))]
  List$Unc <- Unc[order(names(Unc))]
  
  if (verbose) cat("   - Decay Correction finished successfully.\n")
  invisible(List)
}
