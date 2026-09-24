#' @title Constant Ra-226 activity
#' 
#' @description
#' If no \eqn{^{226}}Ra values are available, this function estimates a constant 
#' supported (base) \eqn{^{210}}Pb activity from the bottom-most 
#' \eqn{^{210}}Po or \eqn{^{210}}Pb total activities.
#'
#' @param List Dating list created by \code{ReadData} and processed by other package functions.
#' @param mode Character. Total activity used to determine the constant Ra-226 activity.: 
#' \code{"auto"} (default), \code{"Pb"} or \code{"Po"}. 
#' If both \eqn{^{210}}Pb and \eqn{^{210}}Po are present, the user \strong{must} specify which one to use.
#' @param verbose Logical. Shows estimation details (default \code{TRUE}).
#'
#' @details
#' At least three non-NA values are required.
#' 
#' The algorithm starts from the core bottom and, if needed, progressively includes 
#' more sections to calculate the mean and standard deviation. It stops when 
#' the activity of the section immediately above is significantly higher 
#' than the current mean + 2 SD (indicating the equilibrium zone has been found).
#' 
#' The resulting mean is assigned as constant \code{Ra226} for the entire profile.
#'
#' \strong{Important:} The user should carefully review the result. 
#' This function only handles simple cases. If the profile does not reach 
#' equilibrium, the last value is used as fallback but it is likely wrong.
#' 
#' @returns Updated dating list with constant \eqn{^{226}}Ra and its uncertainty.
#' 
#' @examples
#' \dontrun{
#'   List <- ReadData()
#'   List <- ConstantRa(List, mode = "Po")}
#' 
#' @seealso \code{\link{ReadData}}, \code{\link{DecayCorrection}}, 
#' \code{\link{Equilibrium}}, \code{\link{EquilibriumSet}}
#' 
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

ConstantRa <- function(List, mode = "auto", verbose = TRUE) {
  
  if (verbose) cat("Estimate constant Ra-226 activity: \n")
  
  if (is.null(List$Val$Pb210Total) & is.null(List$Val$Po210Total)) {
    stop("   - Pb-210 and Po-210 activities are missing. 
         I cannot estimate a constant activity of Ra-226.")
  }
  
  # ------------------------------------------------------------------
  # Autodetect or enforce mode
  # ------------------------------------------------------------------
  hasPb <- !all(is.na(List$Val$Pb210Total))
  hasPo <- !all(is.na(List$Val$Po210Total))
  mode  <- match.arg(tolower(mode), c("auto", "pb", "po"))
  
  if (mode == "auto") {
    if (hasPb && !hasPo) { 
      mode <- "pb"
      if (verbose) cat("   - Detected total Pb-210 for constant Ra-226 estimation.\n")
    } else if (hasPo && !hasPb) {
      mode <- "po"
      if (verbose) cat("   - Detected total Po-210 for constant Ra-226 estimation.\n")
    } else if (hasPb && hasPo) {
      stop("   - Both total Pb-210 and Po-210 activities are present.\n",
           "Please specify mode = 'Pb' or mode = 'Po'.")
    } else {
      stop("   - No usable Pb-210 or Po-210 data found.")
    }
  }
  
  # Select the chosen activity
  if (mode == "pb") {
    Tot   <- List$Val$Pb210Total
    TotU  <- List$Unc$Pb210Total
    nucl  <- "Pb-210"
  } else {
    Tot   <- List$Val$Po210Total
    TotU  <- List$Unc$Po210Total
    nucl  <- "Po-210"
  }
  
  # ------------------------------------------------------------------
  # Find equilibrium zone
  # ------------------------------------------------------------------
  # Use only non-NA values
  valid_idx  <- which(!is.na(Tot))
  Tot_clean  <- Tot[valid_idx]
  TotU_clean <- TotU[valid_idx]
  N          <- length(Tot_clean)
  
  if (N < 3) {
    stop("Not enough ", nucl, " data (need at least 3 values) to estimate constant Ra226.\n",
         "Please provide Ra226 values manually in the constants or profiles file.")
  }
  
  # Search for equilibrium zone from the bottom
  found <- FALSE
  for (i in (N-1):2) {
    Ra  <- mean(Tot_clean[i:N])
    RaU <- sqrt(sd(Tot_clean[i:N])^2 + sqrt(sum(TotU_clean[i:N]^2))^2)
    
    # Check if section above is significantly higher
    if (Ra + 2 * RaU < Tot_clean[i - 1]) {
      found <- TRUE
      break
    }
  }
  
  # if not found, use bottom value
  if (!found) {
    Ra  <- Tot_clean[N]
    RaU <- TotU_clean[N]
    warning("No clear equilibrium zone found. Using the last ", nucl, " value as Ra-226.")
  }
  
  # Assign to List
  List$Val$Ra226 <- rep(Ra, List$Val$NSections)
  List$Unc$Ra226 <- rep(RaU, List$Val$NSections)
  
  # Initialize also the Complete field
  List$Val$Ra226Complete <- List$Val$Ra226
  List$Unc$Ra226Complete <- List$Unc$Ra226
  
  if (verbose) {    
    cat("   - Constant Ra-226 estimated:", round(Ra, 2), "+-", round(RaU, 2), "Bq/kg.\n")
  }
  
  # Re-sort lists
  List$Val <- List$Val[order(names(List$Val))]
  List$Unc <- List$Unc[order(names(List$Unc))]
  
  invisible(List)
}
