#' @title Pb-210 dating of sediment cores with the CFCS model
#' 
#' @description
#' Performs \eqn{^{210}}Pb dating using the Constant Flux and Constant Sedimentation (CFCS) model,
#' with optional piece-wise segments, although this option is not yet fully tested. 
#' The user should consider if the profile can meet simultaneoulsy 
#' the hypothesis of constant excess \eqn{^{210}}Pb flux and constant sedimentation.
#' In principal, the user must indicate to which depth the model can be used.
#' Although debatable, some authors use a piece-wise model (by segments),
#' and the user can also provide the proposed break points for the model.
#' However, one must take into account that, in this case,
#' the constant flux hypothesis is only maintained in each profile segment. 
#' 
#' When selecting the break sections, the user should be aware of the limitations of
#' linear regression analysis. This very much depends on the dataset and the expected robustness,
#' but we think that a linear regression of less than 5 points may usually be not statitiscally significant.
#'
#' Please note the required contents of the Breaks() vector,
#' which also defines the maximum layer to which CFCS fitting will be done.
#'
#' @param List Dating list created by \code{ReadData} and processed by other package functions.
#' @param Breaks In order to be as clear as possible about the user intention,
#' the vector must contain the following layers (section bottom depth):
#' \itemize{
#' \item{The surface layer (0).}
#' \item{The following intermediate break points (if any).}
#' \item{The depth until which fitting will take place.
#' Its maximum value is the layer down to which excess \eqn{^{210}}Pb is observed. 
#' It  can be last break point.}
#' The minimum number of elements of Breaks is 2, corresponding to the surface layer and
#' the end of the segment to be used for CFCS fitting.}
#' @param MonteCarloRunsCFCS Number of Monte Carlo iterations for uncertainty. Default = 10000.
#' @param verbose Logical. If TRUE, prints progress logs to the console (default \code{TRUE}).
#'
#' @returns Updated dating list with CFCS results in $Val and uncertainties in $Unc.
#' 
#' @examples
#' \dontrun{
#' List <- Pb210CFCS(List, Breaks = c(0, "your_core_ dating_layer"))}
#' 
#' @seealso \code{\link{Pb210CF}}, \code{\link{AgeModel}}, 
#' \code{\link{MissingInventory}}, \code{\link{Output}}. 
#'
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

Pb210CFCS <- function(List, 
                      Breaks = NA, 
                      MonteCarloRunsCFCS = 1e4, 
                      verbose = TRUE) {
  
  if(verbose) cat("Start the CFCS Pb-210 dating model.\n")
  
  Val <- List$Val
  Unc <- List$Unc
  Val$MonteCarloRunsCFCS <- MonteCarloRunsCFCS
  
  # Initialization of missing inventories
  if (is.null(Val$MissingInventoryBottom) || is.na(Val$MissingInventoryBottom)) Val$MissingInventoryBottom <- 0
  if (is.null(Unc$MissingInventoryBottom) || is.na(Unc$MissingInventoryBottom)) Unc$MissingInventoryBottom <- 0
  if (is.null(Val$MissingInventoryTop)    || is.na(Val$MissingInventoryTop))    Val$MissingInventoryTop    <- 0
  if (is.null(Unc$MissingInventoryTop)    || is.na(Unc$MissingInventoryTop))    Unc$MissingInventoryTop    <- 0
  
  # If LayerDating not supplied, use whole profile
  if (!"LayerDating" %in% names(Val) || is.na(Val$LayerDating) || is.null(Val$LayerDating)) {
    if (verbose) cat("   - No LayerDating is found. The whole profile is used.\n")
    Val$LayerDating <- tail(Val$DepthLayerBottom, 1)
  }
  
  # Handle Breaks
  if (all(is.na(Breaks))) {
    BreakPoints <- c(0, Val$LayerDating)
  } else {
    BreakPoints <- sort(Breaks)
    if (BreakPoints[1] != 0) {
      if (verbose) cat("   - Break = 0 cm has been added.\n")
      BreakPoints <- c(0, BreakPoints)
    }
  }
  
  if (length(BreakPoints) < 2) stop("   - At least 2 break points (including 0) must be provided.")
  if (max(BreakPoints) > Val$LayerDating) stop("   - Last break is beyond the datable core segment.")
  if (sum(BreakPoints %in% Val$DepthLayerBottom) < length(BreakPoints) - 1) 
    stop("   - Breaks do not correspond to bottom layer depths.")
  Val$BreakPoints <- BreakPoints
  
  # Index of the last datable section (LayerDating is a depth value)
  idx <- which(Val$DepthLayerBottom == Val$LayerDating)
  if (length(idx) == 0) stop("LayerDating (", Val$LayerDating, " cm) not found in DepthLayerBottom.")
  Val$NumberOfDatedSections <- idx[1]
  keep_idx <- 1:Val$NumberOfDatedSections
  
  SectionIndicesAtBreaks <- c(1, which(Val$DepthLayerBottom %in% BreakPoints))
  NumberOfCFCSsegments   <- length(BreakPoints) - 1
  
  # Excess Pb210 
  if ("Pb210ExcessDecay" %in% names(Val) && !all(is.na(Val$Pb210ExcessDecay[keep_idx]))) {
    Val$Pb210 <- Val$Pb210ExcessDecay[keep_idx]
    Unc$Pb210 <- Unc$Pb210ExcessDecay[keep_idx]
    if (verbose) cat("   - Using decay corrected excess complete Pb-210.\n")
  } else {
    stop("The CFCS model requires Pb210ExcessDecay. Please use DecayCorrection() first.")
  }
  
  #### Core dating calculations ####
  Val$DecayConstantPb210 <- log(2) / Val$HalfLifePb210
  Unc$DecayConstantPb210 <- Val$DecayConstantPb210 * Unc$HalfLifePb210 / Val$HalfLifePb210
  
  Val$CoreSurfaceArea <- pi * (Val$Diameter / 2)^2 / 1e4
  Unc$CoreSurfaceArea <- pi / 2 * Val$Diameter * Unc$Diameter / 1e4
  
  Val$DepthLayerExtended <- c(0, Val$DepthLayerBottom)
  Unc$DepthLayerExtended <- c(0, Unc$DepthLayerBottom)
  
  Val$SectionThickness <- diff(Val$DepthLayerExtended)
  Unc$SectionThickness <- sqrt(Unc$DepthLayerExtended[2:(Val$NSections + 1)]^2 +
                                 Unc$DepthLayerExtended[1:Val$NSections]^2)
  
  Val$DepthSectionMean <- RollMean2V(Val$DepthLayerExtended)
  Unc$DepthSectionMean <- sqrt(Unc$DepthLayerExtended[2:(Val$NSections + 1)]^2 +
                                 Unc$DepthLayerExtended[1:Val$NSections]^2)
  
  Val$DensityDry <- Val$MassSection / Val$CoreSurfaceArea / Val$SectionThickness / 1e4
  Unc$DensityDry <- Val$DensityDry * sqrt((Unc$MassSection / Val$MassSection)^2 +
                                            (Unc$CoreSurfaceArea / Val$CoreSurfaceArea)^2 +
                                            (Unc$SectionThickness / Val$SectionThickness)^2)
  
  Val$MassLayerCumulative   <- c(0, cumsum(Val$MassSection))
  Val$MassSectionCumulative <- RollMean2V(Val$MassLayerCumulative)
  Unc$MassSectionCumulative <- Unc$MassSection
  Unc$MassSectionCumulative[1] <- Unc$MassSection[1]
  for (i in 2:Val$NSections) {
    Unc$MassSectionCumulative[i] <- sqrt(Unc$MassSectionCumulative[i - 1]^2 +
                                           Unc$MassSection[i]^2)
  }
  
  Val$MassDepthSectionCumulative <- Val$MassSectionCumulative /
    Val$CoreSurfaceArea / 1e4
  Unc$MassDepthSectionCumulative <- Val$MassDepthSectionCumulative *
    sqrt((Unc$MassSectionCumulative / Val$MassSectionCumulative)^2 +
           (Unc$CoreSurfaceArea / Val$CoreSurfaceArea)^2)
  
  Val$AgeCFCS  <- rep(NA, Val$NSections)
  Val$YearCFCS <- rep(NA, Val$NSections)
  
  Time0 <- 0
  MassDepth0 <- 0
  
  for (i in 1:NumberOfCFCSsegments) {
    idx <- SectionIndicesAtBreaks[i]:SectionIndicesAtBreaks[i + 1]
    C <- Val$Pb210[idx]                     # canonical excess
    MassDepth <- Val$MassDepthSectionCumulative[idx]
    
    MarReg <- lm(log(C) ~ MassDepth, na.action = na.exclude)
    MarSum <- summary(MarReg)
    MarSlope <- MarSum$coefficients[2, 1]
    Val$MassAccumulationRateMean[i] <- -Val$DecayConstantPb210 / MarSlope
    Val$PValueMassAccumulationRate[i] <- MarSum$coefficients[2, 4]
    
    SarReg <- lm(log(C) ~ Val$DepthSectionMean[idx], na.action = na.exclude)
    SarSum <- summary(SarReg)
    SarSlope <- SarSum$coefficients[2, 1]
    Val$SedimentAccumulationRateMean[i] <- -Val$DecayConstantPb210 / SarSlope
    Val$PValueSedimentAccumulationRate[i] <- SarSum$coefficients[2, 4]
    
    Val$AgeCFCS[idx] <- (MassDepth - MassDepth0) /
      Val$MassAccumulationRateMean[i] + Time0
    Val$YearCFCS[idx] <- decimal_date(Val$DateSampling) - Val$AgeCFCS[idx]
    
    Time0 <- Val$AgeCFCS[SectionIndicesAtBreaks[i + 1]]
    MassDepth0 <- Val$MassDepthSectionCumulative[SectionIndicesAtBreaks[i + 1]]
  }
  
  #Val$AgeSectionMean <- Val$AgeCFCS
  
  #### Monte Carlo uncertainty ####
  if (Val$MonteCarloRunsCFCS <= 0) {
    Val <- Val[order(names(Val))]
    Unc <- Unc[order(names(Unc))]
    if (verbose) cat("   - CFCS model finished (Monte Carlo skipped).\n")
    return(list(Val = Val, Unc = Unc))
  }
  
  if (verbose) cat("   - Running Monte Carlo uncertainty estimation (", 
                   Val$MonteCarloRunsCFCS, " iterations)...\n")
  
  # suppress spurious MC warnings
  old_warn <- getOption("warn")
  options(warn = -1)
  
  # Preallocate
  n <- Val$MonteCarloRunsCFCS
  n_seg <- NumberOfCFCSsegments
  
  Ran <- vector("list", 10)
  names(Ran) <- c("MassAccumulationRateMean", "SedimentAccumulationRateMean",
                  "PValueMassAccumulationRate", "PValueSedimentAccumulationRate",
                  "AgeCFCS", "HalfLifePb210", "DecayConstantPb210",
                  "MassDepthSectionCumulative", "Pb210", "DepthSectionMean")
  
  Ran$MassAccumulationRateMean <- matrix(NA_real_, nrow = n, ncol = n_seg)
  Ran$SedimentAccumulationRateMean <- matrix(NA_real_, nrow = n, ncol = n_seg)
  Ran$PValueMassAccumulationRate <- matrix(NA_real_, nrow = n, ncol = n_seg)
  Ran$PValueSedimentAccumulationRate <- matrix(NA_real_, nrow = n, ncol = n_seg)
  Ran$AgeCFCS <- matrix(NA_real_, nrow = n, ncol = Val$NSections)
  
  # Safe rnorm
  safe_rnorm <- function(mean_val, sd_val, n) {
    sd_val[is.na(sd_val) | sd_val <= 0] <- 1e-6   # small positive sd
    rnorm(n, mean_val, sd_val)
  }
  
  # Randomize inputs
  Ran$Pb210 <- mapply(safe_rnorm, Val$Pb210, Unc$Pb210, n)
  Ran$MassDepthSectionCumulative <- mapply(safe_rnorm, 
                                           Val$MassDepthSectionCumulative,
                                           Unc$MassDepthSectionCumulative, n)
  Ran$DepthSectionMean <- mapply(safe_rnorm, 
                                 Val$DepthSectionMean,
                                 Unc$DepthSectionMean, n)
  
  Ran$HalfLifePb210 <- safe_rnorm(Val$HalfLifePb210, Unc$HalfLifePb210, n)
  Ran$DecayConstantPb210 <- log(2) / Ran$HalfLifePb210
  
  Ran$Pb210[Ran$Pb210 <= 0] <- NA
  
  # Segment-wise Monte Carlo
  Time0_mean <- 0
  Time0_unc  <- 0
  MassDepth0_mean <- 0
  MassDepth0_unc  <- 0
  
  for (i in 1:n_seg) {
    idx <- SectionIndicesAtBreaks[i]:SectionIndicesAtBreaks[i + 1]
    
    Ran$Time0      <- safe_rnorm(Time0_mean, Time0_unc, n)
    Ran$MassDepth0 <- safe_rnorm(MassDepth0_mean, MassDepth0_unc, n)
    
    for (j in 1:n) {
      Cj <- Ran$Pb210[j, idx]
      Mj <- Ran$MassDepthSectionCumulative[j, idx]
      
      # Skip bad iterations (not enough valid points for regression)
      if (sum(!is.na(Cj) & Cj > 0) < 2) next
      
      # MAR regression
      MarReg <- tryCatch(lm(log(Cj) ~ Mj, na.action = na.exclude), 
                         error = function(e) NULL)
      if (!is.null(MarReg) && nrow(summary(MarReg)$coefficients) >= 2) {
        MarSum <- summary(MarReg)
        Ran$MassAccumulationRateMean[j, i] <- -Ran$DecayConstantPb210[j] / 
          MarSum$coefficients[2, 1]
        Ran$PValueMassAccumulationRate[j, i] <- MarSum$coefficients[2, 4]
      }
      
      # SAR regression
      Sj <- Ran$DepthSectionMean[j, idx]
      SarReg <- tryCatch(lm(log(Cj) ~ Sj, na.action = na.exclude), 
                         error = function(e) NULL)
      if (!is.null(SarReg) && nrow(summary(SarReg)$coefficients) >= 2) {
        SarSum <- summary(SarReg)
        Ran$SedimentAccumulationRateMean[j, i] <- -Ran$DecayConstantPb210[j] / 
          SarSum$coefficients[2, 1]
        Ran$PValueSedimentAccumulationRate[j, i] <- SarSum$coefficients[2, 4]
      }
    }
    
    # Age calculation - protect against division by zero / NA
    valid_mar <- !is.na(Ran$MassAccumulationRateMean[, i]) & 
      Ran$MassAccumulationRateMean[, i] > 0
    
    if (any(valid_mar)) {
      Ran$AgeCFCS[valid_mar, idx] <- (Ran$MassDepthSectionCumulative[valid_mar, idx] - 
                                        Ran$MassDepth0[valid_mar]) / 
        Ran$MassAccumulationRateMean[valid_mar, i] + Ran$Time0[valid_mar]
    }
    
    # Update for next segment
    Time0_mean <- Val$AgeCFCS[SectionIndicesAtBreaks[i + 1]]
    Time0_unc  <- sd(Ran$AgeCFCS[, SectionIndicesAtBreaks[i + 1]], na.rm = TRUE)
    MassDepth0_mean <- Val$MassDepthSectionCumulative[SectionIndicesAtBreaks[i + 1]]
    MassDepth0_unc  <- sd(Ran$MassDepthSectionCumulative[, SectionIndicesAtBreaks[i + 1]], 
                          na.rm = TRUE)
  }
  
  # Final summaries
  Unc$MassAccumulationRateMean <- apply(Ran$MassAccumulationRateMean, 2, sd, na.rm = TRUE)
  Unc$SedimentAccumulationRateMean <- apply(Ran$SedimentAccumulationRateMean, 2, sd, na.rm = TRUE)
  Unc$AgeCFCS <- apply(Ran$AgeCFCS, 2, sd, na.rm = TRUE)
  Unc$YearCFCS <- Unc$AgeCFCS
  
  #Unc$AgeSectionMean <- Unc$AgeCFCS
  
  Val$PValueMassAccumulationRate <- apply(Ran$PValueMassAccumulationRate, 2, mean, na.rm = TRUE)
  Val$PValueSedimentAccumulationRate <- apply(Ran$PValueSedimentAccumulationRate, 2, mean, na.rm = TRUE)
  
  # Finalize
  List$Val <- Val[order(names(Val))]
  List$Unc <- Unc[order(names(Unc))]
  if (verbose) cat("   - CFCS model finished successfully.\n")
  options(warn = old_warn)
  
  invisible(List)
}
