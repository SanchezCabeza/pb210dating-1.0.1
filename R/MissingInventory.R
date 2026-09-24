#' @title Calculate missing inventory using a reference depth or
#' sediment accumulation rate 
#'
#' @description
#' In a sediment core dated with the \code{\link{Pb210CF}} model,
#' equilibrium may not be reached at the core bottom. 
#' This function estimates the missing \eqn{^{210}}Pb\eqn{_{ex}} inventory
#' (below the deepest layer) so that the Constant Flux (CF) model
#' can be applied.
#'
#' Three estimation modes are available:
#'\itemize{
#' \item{\strong{reference_date} (or any string containing "ref" or "dat"): Uses an
#'   independent chronomarker from the validation file.}
#' \item{\strong{mean_mar} (or any string containing "mea"): Uses the mean mass
#'   accumulation rate of the deepest \code{Bottom_n} sections excluding
#'   the last section (which always has MAR = 0 when no missing
#'   inventory is yet assumed).}
#' \item{\strong{fitted_mar} (or any string containing "fit"): Fits \eqn{\log(C)}
#'   versus mass depth on the deepest \code{Bottom_n} sections.}}
#'
#' At the end, the function re-runs \code{\link{Pb210CF}(List)}
#' with the updated \code{MissingInventoryBottom} (or \code{MissingInventoryTop})
#' to provide consistent ages and accumulation rates.
#'
#' @param List Dating list returned by \code{\link{ReadData}} (can be raw
#'   or already processed by \code{\link{Pb210CF}}).
#' @param Mode Character. One of \code{reference_date} (default),
#'   \code{mean_mar}, or \code{fitted_mar}. Partial matching is
#'   supported (e.g. \code{ref}, \code{dat}, \code{mea}, \code{fit}).
#' @param ValidationNumber Integer. Which validation point (row in the
#'   validation table) to use when \code{Mode = "reference_date"}.
#'   Ignored for the other modes. Default is 1.
#' @param Bottom_n Integer. Number of deepest dated sections to use for
#'   MAR calculation when \code{Mode = "mean_mar"} or \code{"fitted_mar"}.
#'   Must be \eqn{\ge 1} for \code{"mean_mar"} and \eqn{\ge 2} for
#'   \code{"fitted_mar"}. Default is 3.
#' @param MonteCarloRunsCF Number of Monte Carlo iterations for uncertainty.
#'   Default is 20000. Use 0 to skip uncertainty calculations.
#' @param verbose Logical. If TRUE, prints progress logs (default TRUE).
#'
#' @returns Updated dating list with the missing inventory in
#'   \code{Val$MissingInventoryBottom} (or \code{Val$MissingInventoryTop})
#'   and the estimated Monte Carlo uncertainty in the corresponding \code{Unc$} entry.
#'  
#' @examples  
#' \dontrun{
#' List <- MissingInventory(List, Mode = "ref")}
#'
#' @seealso \code{\link{Equilibrium}}, \code{\link{Pb210CF}}, \code{\link{Pb210CFCS}}
#'
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

MissingInventory <- function(List, Mode = "reference_date", 
                             ValidationNumber = 1, Bottom_n = 3, MonteCarloRunsCF = 2e5,
                             verbose = TRUE) {
  
  if(verbose) cat("Estimation of missing inventory for incomplete profiles. \n")
  
  # ----------------------------------------------------------------------
  # Normalise Mode
  # ----------------------------------------------------------------------
  Mode <- tolower(trimws(Mode))
  if (grepl("ref|dat", Mode)) Mode <- "reference_date"
  if (grepl("mea", Mode)) Mode <- "mean_mar"
  if (grepl("fit", Mode)) Mode <- "fitted_mar"
  Mode <- match.arg(Mode, c("reference_date", "mean_mar", "fitted_mar"))
  
  # ----------------------------------------------------------------------
  # If not yet CF dated, run Pb210CF silently
  # ----------------------------------------------------------------------
  Val <- List$Val
  Unc <- List$Unc
  
  if (!"FluxPb210" %in% names(Val)) {
    if (verbose) cat("   - List has not been dated yet: running Pb210CF(List) internally...\n")
    List <- Pb210CF(List, verbose = FALSE)
    Val  <- List$Val
    Unc  <- List$Unc
  }
  
  # Initialize missing inventories
  Val$MissingInventoryBottom <- Val$MissingInventoryTop <- 
    Val$MissingInventoryBottomBq <- Val$MissingInventoryTopBq <- 0

  # excess Pb210 
  if ("Pb210ExcessDecayComplete" %in% names(Val) && !all(is.na(Val$Pb210ExcessDecayComplete))) {
    Val$Pb210 <- Val$Pb210ExcessDecayComplete
    Unc$Pb210 <- Unc$Pb210ExcessDecayComplete
    if (verbose) cat("   - Using decay corrected excess complete Pb-210.\n")
  } else {
    stop("The CF model requires Pb210ExcessDecayComplete. Please use DecayCorrection() and Complete().")
  }
  
  # ----------------------------------------------------------------------
  # Input validation
  # ----------------------------------------------------------------------
  if (Mode %in% c("mean_mar", "fitted_mar")) {
    if (!"NumberOfDatedSections" %in% names(Val) || !"Pb210Excess" %in% names(Val)) {
      stop("For MAR-based modes List must have been processed by Pb210CF(...) first.")
    }
    if (Mode == "mean_mar" && Bottom_n < 1) stop("For mean_mar, Bottom_n must be at least 1.")
    if (Mode == "fitted_mar" && Bottom_n < 2) stop("For fitted_mar, Bottom_n must be at least 2.")
    
    N_dated <- Val$NumberOfDatedSections
    if (Bottom_n >= N_dated) stop("Bottom_n cannot be equal to or larger than the number of dated sections.")
  }
  
  # ----------------------------------------------------------------------
  # REFERENCE DATE
  # ----------------------------------------------------------------------
  if (Mode == "reference_date") {
    
    DepthRef   <- Val$ValidationDepthMean [ValidationNumber]
    DepthLower <- Val$ValidationDepthLower[ValidationNumber]
    DepthUpper <- Val$ValidationDepthUpper[ValidationNumber]
    YearRef    <- Val$ValidationYearMean  [ValidationNumber]
    YearLower  <- Val$ValidationYearLower [ValidationNumber]
    YearUpper  <- Val$ValidationYearUpper [ValidationNumber]
    
    # Check values
    if (YearRef <= 1900 & verbose)  cat ("   - The reference date is too old; the estimation of the missing inventory is not reliable.\n")
    if (DepthRef > Val$LayerDating) stop("   - The reference depth is beyond the dating range.")
    if (YearUpper - YearLower < 0) {
      if (verbose) cat("   - Reference dates are inverted and have been swapped.\n")
      A <- YearUpper; YearUpper <- YearLower; YearLower <- A
    }
    if (DepthLower - DepthUpper < 0) {
      if (verbose) cat("   - Depths are inverted and have been swapped.\n")
      A <- DepthUpper; DepthUpper <- DepthLower; DepthLower <- A
    }
    
    NLayer     <- min(which(Val$DepthLayerBottom >= DepthRef))
    Difference <- Val$YearSectionMean[NLayer] - YearRef
    if (Difference == 0) stop("   - Validation is exact. Missing inventory is not needed.\n")
    Type       <- YearRef - Val$YearSectionMean[NLayer] > 0
    
    if (Type) {
      if (verbose) cat("   - Bottom inventory is calculated.\n")
      YearCF <- function(x) {
        List$Val$MissingInventoryBottom <- x
        Pb210CF(List, MonteCarloRunsCF = 0, verbose = FALSE)$Val$YearSectionMean[NLayer]
      }
      f.MissingInventoryBottom <- seq(0, 3 * Val$TotalInventory, length.out = 1000)
      f.YearCF   <- lapply(f.MissingInventoryBottom, YearCF)
      f.YearDiff <- unlist(f.YearCF) - YearRef
      
      fun.MissingInventoryBottom <- approxfun(f.YearDiff, f.MissingInventoryBottom)
      Val$MissingInventoryBottom <- fun.MissingInventoryBottom(0)
      
    } else {
      cat("   - Top inventory is needed. This needs a strong justification.\n")
      YearCF <- function(x) {
        List$Val$MissingInventoryTop <- x
        Pb210CF(List, MonteCarloRunsCF = 0, verbose = FALSE)$Val$YearSectionMean[NLayer]
      }
      f.MissingInventoryTop <- seq(0, 3 * Val$TotalInventory, length.out = 1000)
      f.YearCF   <- lapply(f.MissingInventoryTop, YearCF)
      f.YearDiff <- unlist(f.YearCF) - YearRef
      fun.MissingInventoryTop <- approxfun(f.YearDiff, f.MissingInventoryTop)
      Val$MissingInventoryTop <- fun.MissingInventoryTop(0)
    }
  }
  
  # ----------------------------------------------------------------------
  # USING MEAN_MAR
  # ----------------------------------------------------------------------
  if (Mode == "mean_mar") {
    idx_bottom <- (Val$NumberOfDatedSections - Bottom_n) : (Val$NumberOfDatedSections - 1)
    C_bottom   <- Val$Pb210Excess[idx_bottom]
    C_last     <- Val$Pb210Excess[Val$NumberOfDatedSections]
    C_last_u   <- Unc$Pb210Excess[Val$NumberOfDatedSections]
    lambda     <- log(2) / Val$HalfLifePb210
    lambda_u   <- lambda * Unc$HalfLifePb210 / Val$HalfLifePb210
    if (any(C_bottom <= 0, na.rm = TRUE)) warning("Some bottom concentrations are <= 0; results may be unreliable.")
    
    # iteration for mar convergence
    delta <- 1e6
    MIBot <- numeric()
    for (i in 1:100) {
      mar_values   <- Val$MassAccumulationRate[idx_bottom]
      mar_values_u <- Unc$MassAccumulationRate[idx_bottom]
      valid_mar    <- mar_values  [is.finite(mar_values) & mar_values > 0]
      valid_mar_u  <- mar_values_u[is.finite(mar_values) & mar_values > 0]
      MAR_gcm2yr <- if (length(valid_mar) == 0) {
        Val$MassAccumulationRate[Val$NumberOfDatedSections - 1]  
      } else {
        mean(valid_mar, na.rm = TRUE) 
      }
      MAR_gm2yr  <- MAR_gcm2yr * 10000
      MIBot[i] <- Val$MissingInventoryBottom <- 
        (C_last / 1000) * (MAR_gm2yr / lambda)
      # New dating
      List <- list(Val = Val, Unc = Unc)
      List <- Pb210CF(List, MonteCarloRunsCF = 0, verbose = FALSE)
      Val <- List$Val
      Unc <- List$Unc
      if (i > 1) delta <- MIBot[i] - MIBot[i-1]
      if (delta < 1e-3 * MIBot[i]) break
    }
    if (verbose) cat("   - Bottom missing inventory calculated from the mean MAR of the", Bottom_n, "sections before the last one.\n")
  }
  
  # ----------------------------------------------------------------------
  # USING FITTED_MAR
  # ----------------------------------------------------------------------
  if (Mode == "fitted_mar") {
    idx_bottom <- (Val$NumberOfDatedSections - Bottom_n + 1) : (Val$NumberOfDatedSections)
    C_bottom   <- Val$Pb210Excess[idx_bottom]
    C_last     <- Val$Pb210Excess[Val$NumberOfDatedSections]
    C_last_u   <- Unc$Pb210Excess[Val$NumberOfDatedSections]
    lambda     <- log(2) / Val$HalfLifePb210
    lambda_u   <- lambda * Unc$HalfLifePb210 / Val$HalfLifePb210
    
    if (any(C_bottom <= 0, na.rm = TRUE)) warning("Some bottom concentrations are <= 0; results may be unreliable.")
    
    Mass_bottom <- Val$MassDepthSectionCumulative[idx_bottom]
    
    fit       <- lm(log(C_bottom) ~ Mass_bottom, na.action = na.exclude)
    slope     <- summary(fit)$coefficients[2,1]
    slope_u   <- summary(fit)$coefficients[2,2]
    slope_p   <- summary(fit)$coefficients[2,4]
    
    # SAFE handling for NA/NaN p-value (e.g. perfect fit or singular matrix)
    if (is.na(slope_p) || is.nan(slope_p)) {
      slope_p <- 1.0   # treat as non-significant
      if (verbose) cat("   - Package warning: p-value was NA/NaN (likely singular fit). Proceeding with calculation.\n")
    }
    
    if (slope_p > 0.05) {
      if (verbose) cat("   - Package warning: p =", round(slope_p,3), ": The fit is not statistically significant, but calculation proceeds.\n")
    }
    
    
    if (slope_p > 0.05) {
      if (verbose) cat("   - Package warning: p =", round(slope_p,3), ": The fit is not statistically significative, but calculation proceeds.\n")
    }
    MAR_gm2yr <- -lambda / slope
    Val$MissingInventoryBottom <- (C_last / 1000) * (MAR_gm2yr / lambda)
    if (verbose) cat("   - Bottom missing inventory calculated from fitted MAR on the last", Bottom_n, "sections.\n")
  }
  if (Val$MissingInventoryBottom < 0) {
    if (verbose) cat("   - Calculated MissingInventoryBottom is negative; it has been set to 0.")
    Val$MissingInventoryBottom <- 0
  }
  
  # ======================================================================
  # MONTE CARLO UNCERTAINTY
  # ======================================================================
  if (MonteCarloRunsCF > 0) {
    if (verbose) cat("   - Running Monte Carlo uncertainty estimation (",MonteCarloRunsCF, " iterations)...\n")
    isBottom  <- !is.null(Val$MissingInventoryBottom) && !is.na(Val$MissingInventoryBottom)
    
    # Pre-allocate Ran$ (for Random)
    Ran <- list(MissingInventoryBottom = numeric(MonteCarloRunsCF))
    if (Mode == "reference_date") {
      Ran <- list(YearRef                = numeric(MonteCarloRunsCF),
                  YearCF                 = numeric(MonteCarloRunsCF),
                  Difference             = numeric(MonteCarloRunsCF))
    } else {
      Ran <- list(MAR_gm2yr              = numeric(MonteCarloRunsCF),
                  C_last                 = numeric(MonteCarloRunsCF),
                  lambda                 = numeric(MonteCarloRunsCF), 
                  slope                  = numeric(MonteCarloRunsCF), 
                  MissingInventoryBottom = numeric(MonteCarloRunsCF))
      
      if (Mode == "mean_mar") Ran$valid_mar <- matrix(NA, nrow = MonteCarloRunsCF, ncol = length(idx_bottom))
    }
    
    # ======================================================================
    # MC REFERENCE DATE
    # ======================================================================
    if (Mode == "reference_date") {
      YearCF_point   <- Val$YearSectionMean[NLayer]
      YearCFU        <- Unc$YearSectionMean[NLayer]
      Ran$YearRef    <- runif(MonteCarloRunsCF, Val$ValidationYearLower, Val$ValidationYearUpper)  # square distribution
      Ran$YearCF     <- rnorm(MonteCarloRunsCF, YearCF_point, YearCFU) # normal distribution
      Ran$Difference <- Ran$YearCF - Ran$YearRef
      
      if (isBottom) {
        Ran$MissingInventoryBottom      <- fun.MissingInventoryBottom(Ran$Difference)
        Unc$MissingInventoryBottom <- sd(Ran$MissingInventoryBottom, na.rm = TRUE)
      } else {
        Ran$MissingInventoryTop      <- fun.MissingInventoryTop(Ran$Difference)
        Unc$MissingInventoryTop <- sd(Ran$MissingInventoryTop, na.rm = TRUE)
      }
      
    } else {
      
      # ======================================================================
      # MC MEAN MAR
      # ======================================================================

      # Common to both mean_mar and fitted_mar
      RedFactor <- 1
      Ran$C_last <- rnorm(MonteCarloRunsCF, C_last, C_last_u / RedFactor)
      Ran$lambda <- rnorm(MonteCarloRunsCF, lambda, lambda_u / RedFactor)
      
      if (Mode == "mean_mar") {
        for (i in seq_along(idx_bottom)) {
          Ran$valid_mar[, i] <- rnorm(MonteCarloRunsCF,
                                      valid_mar[i], valid_mar_u[i] / RedFactor)
        }
        Ran$MAR_mean  <- apply(Ran$valid_mar, 1, mean, na.rm = TRUE)
        Ran$MAR_gm2yr <- Ran$MAR_mean * 10000
        
      } else {   
        # from fitted_mar
        Ran$slope <- rnorm(MonteCarloRunsCF, slope, slope_u / RedFactor) 
        Ran$MAR_gm2yr <- - Ran$lambda / Ran$slope
      }
      
      Unc$MissingInventoryBottom <-  sd(Ran$C_last * Ran$MAR_gm2yr / Ran$lambda / 1000 * RedFactor, na.rm = TRUE)
    }
  }
  
  # ======================================================================
  # END
  # ======================================================================

  # Finalize
  List$Val <- Val[order(names(Val))]
  List$Unc <- Unc[order(names(Unc))]
  if (verbose) cat("   - Monte Carlo uncertainty estimation successful. \n")  
  invisible(List)
}
