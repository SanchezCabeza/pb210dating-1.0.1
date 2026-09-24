#' @title Calculate the bottom layer with no excess Pb-210 below
#' 
#' @description
#' Estimates the layer below which there is no measurable excess \eqn{^{210}}Pb 
#' (the "equilibrium layer"). This value is critical for the CF dating model.
#'
#' The function works after \code{\link{DecayCorrection}}, 
#' so the presence of \code{Pb210ExcessDecay} is expected. 
#' 
#' @importFrom graphics barplot polygon text axis par
#'
#' @param List Dating list created by \code{ReadData} and processed by other package functions.
#' @param Interactive Logical. If \code{TRUE} (default), allows interactive visual selection 
#'        of the equilibrium layer (highly recommended). If \code{FALSE}, it uses the automatic rule.
#' @param verbose Logical. If TRUE, prints progress logs to the console (default \code{TRUE}).
#'
#' @details
#' Excess \eqn{^{210}}Pb is shown with 1-sigma uncertainty.
#' In interactive mode the x-axis shows tick marks and labels exactly centered 
#' under each red-blue bar pair. You can zoom in/out.
#'
#' @returns Updated dating list with the value of \code{Val$LayerDating}.
#' 
#' @seealso \code{\link{ConstantRa}}, \code{\link{EquilibriumSet}}, \code{\link{MissingInventory}}, 
#' \code{\link{Pb210CF}}, \code{\link{Pb210CFCS}}
#' 
#' @examples 
#' \dontrun{
#' List <- Equilibrium(List)}
#' 
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

Equilibrium <- function(List, Interactive = TRUE, verbose = TRUE) {
  
  if (verbose) cat("Start searching for Pb-210 - Ra-226 equilibrium. \n")
  
  # Read data
  Val <- List$Val
  Unc <- List$Unc
  
  NSections        <- Val$NSections
  DepthLayerBottom <- Val$DepthLayerBottom
  # excess Pb-210
  if ("Pb210ExcessDecayComplete" %in% names(Val) && !all(is.na(Val$Pb210ExcessDecayComplete))) {
    Val$Pb210 <- Val$Pb210ExcessDecayComplete
    Unc$Pb210 <- Unc$Pb210ExcessDecayComplete
    if (verbose) cat("   - Using decay corrected excess complete Pb-210.\n")
  } else if ("Pb210ExcessDecay" %in% names(Val) && !all(is.na(Val$Pb210ExcessDecay))) {
    Val$Pb210 <- Val$Pb210ExcessDecay
    Unc$Pb210 <- Unc$Pb210ExcessDecay
    if (verbose) cat("   - Using excess decay Pb-210.\n")
  } else {
    # Excess Decay Pb-210 is needed
    stop("  - Pb210ExcessDecay or Pb210ExcessDecayComplete are required. ")
  }
  
  # Matrix for barplot: excess (red) and uncertainty (blue)
  Pb210 <-    Val$Pb210
  Pb210Unc <- Unc$Pb210
  D <- rbind(Pb210, Pb210Unc)
  
  # Equilibrium conditions
  .find_equilibrium <- function() {
    L <- max(DepthLayerBottom, na.rm = TRUE)
    if (all(Pb210 - Pb210Unc > 0, na.rm = TRUE)) {
      if (verbose) cat("   - All excess activities are positive (1 sigma).\n")
    } else if (any(Pb210 - Pb210Unc <= 0, na.rm = TRUE)) {
      L <- min(DepthLayerBottom[which(Pb210 - Pb210Unc <= 0)], na.rm = TRUE)
      if (verbose) cat("   - At least one excess activity is compatible with zero (1 sigma).\n")
    } else if (any(Pb210 == 0, na.rm = TRUE)) {
      L <- min(DepthLayerBottom[which(Pb210 == 0)], na.rm = TRUE)
      if (verbose) cat("   - At least one excess activity is exactly zero.\n")
    } else if (any(Pb210 < 0, na.rm = TRUE)) {
      L <- min(DepthLayerBottom[which(Pb210 < 0)], na.rm = TRUE)
      if (verbose) cat("   - Negative excess activities detected.\n")
    }
    return(L)
  }
  
  LayerDating_proposed <- .find_equilibrium()
  
  # Non-interactive mode
  if (!Interactive) {
    if (verbose) cat("   - The automatically proposed equilibrium layer is", LayerDating_proposed, "cm.\n")
    List$Val$LayerDating <- LayerDating_proposed
    List$Val <- List$Val[order(names(List$Val))]
    return(List)
  }
  
  # Interactive mode
  E1 <- 1
  E2 <- NSections
  
  repeat {
    cat("   - Exploring excess Pb-210 profile (see barplot)... \n")
    par(mar = c(5, 4, 4, 1))   # extra bottom margin for rotated labels
    barplot(height = D[, E1:E2], beside = TRUE, xaxs = "i", yaxs = "i", las = 1,
            col = c("red", "blue"), 
            names.arg = rep("", length(DepthLayerBottom[E1:E2])),
            main = expression({}^{210}*Pb[ex] ~ "(red) and 1" * sigma ~ "uncertainty (blue)"),
            xlab = "Bottom layer depth (cm)", 
            ylab = "Activity (Bq/kg)")
    n_bars <- length(DepthLayerBottom[E1:E2])
    tick_pos <- seq(2, by = 3, length.out = n_bars)
    axis(1, at = tick_pos, labels = DepthLayerBottom[E1:E2], 
         las = 1, cex.axis = 1, tick = TRUE)
    cat("   -> The proposed equilibrium layer is:", LayerDating_proposed, "cm.\n")
    choice <- readline("   Choose action: [a]ccept proposed, [m]anual entry, [z]oom plot: ")
    choice <- tolower(trimws(choice))
    
    if (choice %in% c("a", "accept")) {
      LayerDating <- LayerDating_proposed
      break
    } else if (choice %in% c("m", "manual")) {
      LayerDating <- as.numeric(readline("   - Enter your desired LayerDating value (cm): "))
      break
    } else if (choice %in% c("z", "zoom")) {
      z_min <- as.numeric(readline("     Enter minimum depth to display (cm): "))
      z_max <- as.numeric(readline("     Enter maximum depth to display (cm): "))
      
      i_min <- which.min(abs(DepthLayerBottom - z_min))
      i_max <- which.min(abs(DepthLayerBottom - z_max))
      
      E1 <- max(1, min(i_min, i_max))
      E2 <- min(NSections, max(i_min, i_max))
    } else {
      cat("   - Invalid choice. Please enter 'a', 'm', or 'z'.\n")
    }
  }
  
  List$Val$LayerDating <- LayerDating
  List$Val <- List$Val[order(names(List$Val))]
  if (verbose) cat("   - Equilibrium layer saved as:", LayerDating, "cm.\n")
  return(List)
}
