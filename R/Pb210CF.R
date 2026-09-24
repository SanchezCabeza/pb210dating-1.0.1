#' @title Pb-210 dating of sediment cores with the Constant Flux (CF) model
#' 
#' @description
#' Performs \eqn{^{210}}Pb dating using the Constant Flux (CF) model.
#' Based on Sanchez-Cabeza and Ruiz-Fernandez (2012).
#'
#' @param List Dating list created by \code{ReadData} and processed by other package functions.
#' @param MonteCarloRunsCF Number of Monte Carlo iterations for uncertainty estimation.
#'        Default is 20000. Use \code{MonteCarloRunsCF = 0} to skip uncertainty calculation.
#' @param verbose Logical.Print progress messages. Default is \code{TRUE}.
#'
#' @details
#' Missing inventories at the top (\code{MissingInventoryTop}) and bottom
#' (\code{MissingInventoryBottom}) of the core are supported but their use requires
#' careful justification. The core is assumed to be cylindrical.
#' Uncertainties are propagated via Monte Carlo simulation
#' (Sanchez-Cabeza et al., 2014).
#'
#' @returns Updated dating list. Results in $Val, uncertainties in $Unc.
#' 
#' @examples
#' \dontrun{
#'   List <- ReadData()
#'   List <- Pb210CF(List)}
#' 
#' @seealso \code{\link{CompleteProfile}}, \code{\link{Pb210CFCS}}, 
#' \code{\link{AgeModel}}, \code{\link{Equilibrium}}, \code{\link{EquilibriumSet}}, 
#' \code{\link{MissingInventory}}, \code{\link{Output}}, \code{\link{PlotMARSAR}}.
#'
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

Pb210CF <- function(List, 
                    MonteCarloRunsCF = 2e5, 
                    verbose = TRUE) {
  
  if(verbose) cat("Start CF Pb-210 dating model.\n")
  
  Val <- List$Val
  Unc <- List$Unc
  Val$MonteCarloRunsCF <- MonteCarloRunsCF
  
  # Initialization
  if (is.null(Val$MissingInventoryBottom) || is.na(Val$MissingInventoryBottom)) Val$MissingInventoryBottom <- 0
  if (is.null(Unc$MissingInventoryBottom) || is.na(Unc$MissingInventoryBottom)) Unc$MissingInventoryBottom <- 0
  if (is.null(Val$MissingInventoryTop)    || is.na(Val$MissingInventoryTop))    Val$MissingInventoryTop    <- 0
  if (is.null(Unc$MissingInventoryTop)    || is.na(Unc$MissingInventoryTop))    Unc$MissingInventoryTop    <- 0
  
  RedFactor <- 100
  
  idx <- which(Val$DepthLayerBottom == Val$LayerDating)
  if (length(idx) == 0) stop("LayerDating (", Val$LayerDating, " cm) not found in DepthLayerBottom.")
  Val$NumberOfDatedSections <- idx[1]
  keep_idx <- 1:Val$NumberOfDatedSections
  
  Val$YearSampling <- decimal_date(Val$DateSampling)
  
  # excess Pb210 
  if ("Pb210ExcessDecayComplete" %in% names(Val) && !all(is.na(Val$Pb210ExcessDecayComplete[keep_idx]))) {
    Val$Pb210 <- Val$Pb210ExcessDecayComplete[keep_idx]
    Unc$Pb210 <- Unc$Pb210ExcessDecayComplete[keep_idx]
    if (verbose) cat("   - Using decay corrected excess complete Pb-210.\n")
  } else {
    stop("The CF model requires Pb210ExcessDecayComplete. Please use DecayCorrection() and Complete().")
  }
  
  #### Core dating calculations ####
  Val$DecayConstantPb210 <- log(2) / Val$HalfLifePb210
  Val$CoreSurfaceArea    <- pi * (Val$Diameter / 2)^2 / 1E4
  
  Val$MissingInventoryTopBq    <- Val$MissingInventoryTop    * Val$CoreSurfaceArea
  Unc$MissingInventoryTopBq    <- Unc$MissingInventoryTop    * Val$CoreSurfaceArea
  Val$MissingInventoryBottomBq <- Val$MissingInventoryBottom * Val$CoreSurfaceArea
  Unc$MissingInventoryBottomBq <- Unc$MissingInventoryBottom * Val$CoreSurfaceArea
  
  Val$DepthLayerExtended <- c(0, Val$DepthLayerBottom)
  Unc$DepthLayerExtended <- c(0, Unc$DepthLayerBottom)
  Val$SectionThickness   <- diff(Val$DepthLayerExtended)
  Val$DepthSectionMean   <- RollMean2V(Val$DepthLayerExtended)
  Val$DensityDry         <- Val$MassSection / Val$CoreSurfaceArea / Val$SectionThickness / 1e4
  Val$MassDepthSection   <- Val$MassSection / Val$CoreSurfaceArea
  Val$MassDepthLayerCumulative   <- c(0, cumsum(Val$MassDepthSection))
  Val$MassDepthSectionCumulative <- RollMean2V(Val$MassDepthLayerCumulative)
  
  # Inventory calculations
  A <- Val$Pb210[keep_idx] * Val$MassSection[keep_idx] / 1000
  Val$InventorySection <- c(Val$MissingInventoryTopBq, A, Val$MissingInventoryBottomBq)
  Val$CumulativeInventoryFromTop <- c(0, cumsum(Val$InventorySection))
  Val$CumulativeInventoryFromBottom <- c(rev(cumsum(rev(Val$InventorySection))), 0)
  Val$TotalInventoryBq <- Val$CumulativeInventoryFromBottom[1]
  Val$TotalInventory   <- Val$TotalInventoryBq / Val$CoreSurfaceArea
  Val$FluxPb210        <- Val$DecayConstantPb210 * Val$TotalInventoryBq / Val$CoreSurfaceArea
  
  # Age calculation
  A <- log(1 + Val$CumulativeInventoryFromTop / Val$CumulativeInventoryFromBottom) / Val$DecayConstantPb210
  Val$AgeLayer        <- A[2:(length(A) - 1)]
  Val$YearLayer       <- Val$YearSampling - Val$AgeLayer
  Val$AgeSectionMean  <- RollMean2V(Val$AgeLayer)
  Val$YearSectionMean <- Val$YearSampling - Val$AgeSectionMean
  Val$AgeDelta        <- diff(Val$AgeLayer)
  
  # Accumulation rates
  Val$MassAccumulationRate     <- Val$MassSection[keep_idx] / Val$AgeDelta / Val$CoreSurfaceArea / 1e4
  Val$SedimentAccumulationRate <- Val$SectionThickness[keep_idx] / Val$AgeDelta
  
  # Return if no MC
  if (Val$MonteCarloRunsCF <= 0) {
    Val <- Val[order(names(Val))]
    Unc <- Unc[order(names(Unc))]
    return(list(Val = Val, Unc = Unc))
  }
  
  #### Monte Carlo uncertainty calculation ####
  if (verbose) cat("   - Start Monte Carlo uncertatinty estimation. \n")
  
  # Preallocate Ran$
  Ran <- vector("list", length(c("Diameter", "HalfLifePb210", "MissingInventoryTop", 
                                 "MissingInventoryBottom", "DepthLayerExtended", 
                                 "MassSection", "Pb210")))
  names(Ran) <- c("Diameter", "HalfLifePb210", "MissingInventoryTop", 
                  "MissingInventoryBottom", "DepthLayerExtended", 
                  "MassSection", "Pb210")
  
  # Randomize constants
  Ran$Diameter              <- rnorm(Val$MonteCarloRunsCF, Val$Diameter,              Unc$Diameter / RedFactor)
  Ran$HalfLifePb210         <- rnorm(Val$MonteCarloRunsCF, Val$HalfLifePb210,         Unc$HalfLifePb210 / RedFactor)
  Ran$MissingInventoryTop   <- rnorm(Val$MonteCarloRunsCF, Val$MissingInventoryTop,   Unc$MissingInventoryTop / RedFactor)
  Ran$MissingInventoryBottom<- rnorm(Val$MonteCarloRunsCF, Val$MissingInventoryBottom,Unc$MissingInventoryBottom / RedFactor)
  
  # Randomize vectors
  Ran$DepthLayerExtended <- mapply(rnorm, Val$MonteCarloRunsCF, Val$DepthLayerExtended, Unc$DepthLayerExtended / RedFactor)
  Ran$MassSection        <- mapply(rnorm, Val$MonteCarloRunsCF, Val$MassSection,        Unc$MassSection / RedFactor)
  Ran$Pb210              <- mapply(rnorm, Val$MonteCarloRunsCF, Val$Pb210, Unc$Pb210 / RedFactor)
  
  # Calculated variables in Monte Carlo
  Ran$DecayConstantPb210 <- log(2) / Ran$HalfLifePb210
  Ran$CoreSurfaceArea    <- pi * (Ran$Diameter / 2)^2 / 1E4
  
  Ran$MissingInventoryTopBq    <- Ran$MissingInventoryTop    * Ran$CoreSurfaceArea
  Ran$MissingInventoryBottomBq <- Ran$MissingInventoryBottom * Ran$CoreSurfaceArea
  
  Ran$DepthDelta         <- t(diff(t(Ran$DepthLayerExtended)))
  Ran$DepthSectionMean   <- RollMean2M(Ran$DepthLayerExtended)
  Ran$DensityDry         <- Ran$MassSection / Ran$CoreSurfaceArea / Ran$DepthDelta / 1e4
  Ran$MassDepthSection   <- Ran$MassSection / Ran$CoreSurfaceArea
  
  A <- t(apply(Ran$MassDepthSection, 1, cumsum))
  Ran$MassDepthLayerCumulative <- cbind(0, A)
  Ran$MassDepthSectionCumulative <- RollMean2M(Ran$MassDepthLayerCumulative)
  
  # excess Pb-210 is Ran$Pb210
  A <- Ran$Pb210[, keep_idx] * Ran$MassSection[, keep_idx] / 1000
  Ran$InventorySection <- cbind(Ran$MissingInventoryTopBq, A, Ran$MissingInventoryBottomBq)
  
  A <- t(apply(Ran$InventorySection, 1, cumsum))
  Ran$CumulativeInventoryFromTop <- cbind(0, A)
  
  A <- t(apply(Ran$InventorySection[, ncol(Ran$InventorySection):1], 1, cumsum))
  Ran$CumulativeInventoryFromBottom <- cbind(A[, ncol(A):1], 0)
  
  Ran$TotalInventoryBq <- Ran$CumulativeInventoryFromBottom[,1]
  Ran$TotalInventory   <- Ran$TotalInventoryBq / Ran$CoreSurfaceArea
  Ran$FluxPb210        <- Ran$DecayConstantPb210 * Ran$TotalInventoryBq / Ran$CoreSurfaceArea
  
  A <- log(1 + Ran$CumulativeInventoryFromTop / Ran$CumulativeInventoryFromBottom) / Ran$DecayConstantPb210
  Ran$AgeLayer <- A[, 2:(ncol(A) - 1)]
  Ran$AgeSectionMean <- RollMean2M(Ran$AgeLayer)
  Ran$AgeDelta <- t(diff(t(Ran$AgeLayer)))
  
  Ran$MassAccumulationRate     <- Ran$MassSection[, keep_idx] / 
    Ran$AgeDelta / Ran$CoreSurfaceArea / 1e4
  Ran$SedimentAccumulationRate <- Ran$DepthDelta[, keep_idx] / Ran$AgeDelta
  
  #### Uncertainties ####
  Uncertainty <- function(x) sd(x) * RedFactor
  
  Unc$TotalInventoryBq         <- Uncertainty(Ran$TotalInventoryBq)
  Unc$TotalInventory           <- Uncertainty(Ran$TotalInventory)
  Unc$FluxPb210                <- Uncertainty(Ran$FluxPb210)
  Unc$MissingInventoryTopBq    <- Uncertainty(Ran$MissingInventoryTopBq)
  Unc$MissingInventoryBottomBq <- Uncertainty(Ran$MissingInventoryBottomBq)
  Unc$CoreSurfaceArea          <- Uncertainty(Ran$CoreSurfaceArea)
  
  Unc$Pb210                    <- apply(Ran$Pb210, 2, Uncertainty)
  Unc$DensityDry               <- apply(Ran$DensityDry, 2, Uncertainty)
  Unc$InventorySection         <- apply(Ran$InventorySection, 2, Uncertainty)
  Unc$CumulativeInventoryFromTop  <- apply(Ran$CumulativeInventoryFromTop, 2, Uncertainty)
  Unc$CumulativeInventoryFromBottom <- apply(Ran$CumulativeInventoryFromBottom, 2, Uncertainty)
  Unc$SectionThickness         <- apply(Ran$DepthDelta, 2, Uncertainty)
  Unc$DepthSectionMean         <- apply(Ran$DepthSectionMean, 2, Uncertainty)
  Unc$MassDepthSection         <- apply(Ran$MassDepthSection, 2, Uncertainty)
  Unc$MassDepthLayerCumulative <- apply(Ran$MassDepthLayerCumulative, 2, Uncertainty)
  Unc$MassDepthSectionCumulative <- apply(Ran$MassDepthSectionCumulative, 2, Uncertainty)
  Unc$MassAccumulationRate     <- apply(Ran$MassAccumulationRate, 2, Uncertainty)
  Unc$SedimentAccumulationRate <- apply(Ran$SedimentAccumulationRate, 2, Uncertainty)
  Unc$AgeDelta                 <- apply(Ran$AgeDelta, 2, Uncertainty)
  Unc$AgeLayer                 <- apply(Ran$AgeLayer, 2, Uncertainty)
  Unc$AgeSectionMean           <- apply(Ran$AgeSectionMean, 2, Uncertainty)
  
  # Singularities in accumulation rates
  Unc$MassAccumulationRate[is.infinite(Val$MassAccumulationRate)] <- Inf
  Unc$SedimentAccumulationRate[is.infinite(Val$SedimentAccumulationRate)] <- Inf
  
  # Non-calculated age variables
  Unc$YearLayer      <- Unc$AgeLayer
  Unc$YearSectionMean<- Unc$AgeSectionMean
  
  # Shorten cumulative inventories for final output
  Val$CumulativeInventoryFromTop    <- Val$CumulativeInventoryFromTop[2:(length(Val$CumulativeInventoryFromTop)-1)]
  Unc$CumulativeInventoryFromTop    <- Unc$CumulativeInventoryFromTop[2:(length(Unc$CumulativeInventoryFromTop)-1)]
  Val$CumulativeInventoryFromBottom <- Val$CumulativeInventoryFromBottom[2:(length(Val$CumulativeInventoryFromBottom)-1)]
  Unc$CumulativeInventoryFromBottom <- Unc$CumulativeInventoryFromBottom[2:(length(Unc$CumulativeInventoryFromBottom)-1)]
  
  # Finalize
  List$Val <- Val[order(names(Val))]
  List$Unc <- Unc[order(names(Unc))]
  if (verbose) cat("   - CF model finished succesfully.\n")  
  invisible(List)
}
