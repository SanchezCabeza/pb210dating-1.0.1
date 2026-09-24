#' @title Read Data for Pb-210 Dating of Sediment Cores
#'
#' @description
#' This function serves as the primary data entry point for the \code{pb210dating} package. 
#' It allows users to either load the built-in example dataset (\code{TehuaII}) or import 
#' user files in CSV format. The function performs automatic cleaning of 
#' empty rows and validates that the 21-column structure meets the requirements 
#' for radiochronological modeling. 
#'
#' @details 
#' The input file must follow a strict \strong{21-column structure}:
#' \itemize{
#'   \item \strong{Columns 1-18 (Profiles):} Sample Code, Depth Bottom (cm), Depth Uncertainty, 
#'   Section Mass (g), Mass Uncertainty, Total Po-210 Activity (Bq/kg), Po-210 Uncertainty, 
#'   Po plating date, Po measurement start, Po measurement end, 
#'   Total Pb-210 Activity (Bq/kg), Pb-210 Uncertainty, Ra-226 Activity, Ra-226 Uncertainty, 
#'   Proxy Value, Proxy Uncertainty, gamma measurement start, gamma measurement end. 
#'   The column name of the Proxy column is used as the Proxy name (e.g., Cs-137). 
#'   \item \strong{Columns 19-21 (Metadata block, 15 rows):} 
#'   \itemize{
#'     \item \strong{Core} Core name (e.g. TEHUAII).
#'     \item \strong{Alpha spectrometry} Analyst(s) who performed alpha spectrometry (full name or initials, e.g. LHPB).
#'     \item \strong{Gamma spectrometry} Analyst(s) who performed gamma spectrometry (full name or initials, e.g. LHPB).
#'     \item \strong{Dating} Analyst(s) responsible for the \eqn{^{210}}Pb dating (full name or initials, e.g. JASC).
#'     \item \strong{Sampling date} yyyy-mm-dd
#'     \item \strong{Internal diameter (cm)} and uncertainty. 
#'     \item \strong{Half life Pb210 (yr)} and uncertainty.
#'     \item \strong{Half life Po210 (yr)} and uncertainty.
#'     \item \strong{Validation label} e.g., \eqn{^{137}}Cs maximum 
#'     \item \strong{Top layer (cm)} up to where the validation point can be present. 
#'     \item \strong{Mean depth (cm)} where the validation point is observed. 
#'     \item \strong{Bottom layer (cm)} down to where the validation point is observed. 
#'     \item \strong{Year upper} up to where the validation point can be present. 
#'     \item \strong{Mean year} where the validation point is observed. 
#'     \item \strong{Year lower} down to where the validation point is observed. 
#'   }
#' }
#'
#' @param NameFile Character. Path to the .csv file. 
#' If omitted or set to NA (default), 
#' the function will load the built-in \code{TehuaII} dataset from the package.
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @returns A list containing two elements: 
#' \item{Val}{A list of numerical variables, dates, and activities.}
#' \item{Unc}{A list containing the uncertainties associated with numeric variables in `Val`.}
#'
#' @examples
#' \dontrun{
#'   List <- ReadData()
#'   List <- ReadData("my_core_data.csv")}
#' 
#' @seealso \code{\link{DecayCorrection}}, 
#' \code{\link{CompleteProfile}}, \code{\link{ConstantRa}}
#' 
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

ReadData <- function(NameFile = NA, verbose = TRUE){    
  
  if (verbose == TRUE) cat("Start reading data. \n")
  
  #### Resolve file path ####
  if (length(NameFile) != 1 || is.na(NameFile) || identical(NameFile, "")) {
    # Load built-in example from inst/extdata
    NameFile <- system.file("extdata", "TehuaII.csv", package = "pb210dating")
    if (!nzchar(NameFile))
      stop("   - Built-in TehuaII.csv not found. Re-install the package or supply a file path.\n")
    if (verbose) cat("   - Using built-in example dataset: TehuaII.csv\n")
  }
  
  if (grepl("\\.csv$", NameFile, ignore.case = TRUE)) {
    T <- read.csv(NameFile, stringsAsFactors = FALSE,
                  na.strings = c("NA", "", " "), check.names = FALSE)
    if (verbose) cat("   - Read user CSV file:", basename(NameFile), "\n")
  } else {
    stop("   - Only CSV files are supported.")
  }

  #### Trim input data file ####
  if (ncol(T) > 21) {
    T <- T[, 1:21, drop = FALSE]
    if (verbose) cat("   - Columns after 21 have been deleted.\n")
  }
  
  # Convert NAs
  T[T == "NA" | T == ""] <- NA
  
  # Remove completely empty rows
  T <- T[rowSums(is.na(T) | T == "") < ncol(T), , drop = FALSE]
  
  # Less columns than expected
  if (ncol(T) < 21) stop("I need at least 16 columns. Please, revise.")
  
  #### Extract profile data (rows with content in column 2) ####
  profile_rows <- which(!is.na(T[[2]]) & trimws(T[[2]]) != "")
  if (length(profile_rows) < 3) 
    stop("I need at least 3 sections. Please check the file.")
  T_profiles <- T[profile_rows, 1:18, drop = FALSE]
  N <- nrow(T_profiles) # Number of sections
  ProxyLabel <- colnames(T)[15]
  
  # Internal column names
  colnames(T_profiles) <- c("SampleCode", "DepthLayerBottom", "DepthLayerU",
                            "MassSection", "MassSectionU", 
                            "Po210Total", "Po210TotalU", 
                            "DatePlating", "MeasurementAlphaStart", "MeasurementAlphaEnd",
                            "Pb210Total", "Pb210TotalU",
                            "Ra226", "Ra226U", 
                            "ProxyValue", "ProxyU", 
                            "MeasurementGammaStart", "MeasurementGammaEnd") 
  
  #### Extract info block (rows 1-15) ####
  T_info <- T[1:15, 19:21, drop = FALSE]
  Val <- list()
  Unc <- list()
  #colnames(T_info) <- c("Name", "Value", "Unc")
  
  ## Assign values by exact row position
  # METADATA FIELDS
  Val$Core                 <- T_info$Value[1]
  Val$AlphaSpectrometry    <- T_info$Value[2]
  Val$GammaSpectrometry    <- T_info$Value[3]
  Val$Dating               <- T_info$Value[4]
  if (is.na(Val$Core) | is.null(Val$Core)) Val$Core <- "Core"
  if (is.null(Val$AlphaSpectrometry))      Val$AlphaSpectrometry <- NA
  if (is.null(Val$GammaSpectrometry))      Val$GammaSpectrometry <- NA
  if (is.null(Val$Dating))                 Val$Dating <- NA
  
  # Row 5: Sampling date
  Val$DateSampling <- as.Date(T_info$Value[5])
  if (is.na(Val$DateSampling) | is.null(Val$DateSampling)) 
    stop("The sampling date (row 5 of metadata block) is needed. Please revise. \n")
  
  # Row 6: Internal diameter
  Val$Diameter <- as.numeric(T_info$Value[6])
  Unc$Diameter <- as.numeric(T_info$Unc[6])
  if (is.na(Val$Diameter) | is.null(Val$Diameter) | is.na(Unc$Diameter) | is.null(Unc$Diameter)) 
    stop("The internal diameter and its uncertainty (row 6) are needed. Please revise. \n")
  
  # Row 7: HalflifePb210
  Val$HalfLifePb210 <- as.numeric(T_info$Value[7])
  Unc$HalfLifePb210 <- as.numeric(T_info$Unc[7])
  if (is.na(Val$HalfLifePb210)) Val$HalfLifePb210 <- 22.23
  if (is.na(Unc$HalfLifePb210)) Unc$HalfLifePb210 <- 0.12
  
  # Row 8: HalflifePo210
  Val$HalfLifePo210 <- as.numeric(T_info$Value[8])
  Unc$HalfLifePo210 <- as.numeric(T_info$Unc[8])
  if (is.na(Val$HalfLifePo210)) Val$HalfLifePo210 <- 0.3788537
  if (is.na(Unc$HalfLifePo210)) Unc$HalfLifePo210 <- 0.0000047
  
  # Rows 9-15: Validation point
  Val$ValidationLabel      <-            T_info$Value[9]
  Val$ValidationDepthUpper <- as.numeric(T_info$Value[10])
  Val$ValidationDepthMean  <- as.numeric(T_info$Value[11])
  Val$ValidationDepthLower <- as.numeric(T_info$Value[12])
  Val$ValidationYearUpper  <- as.numeric(T_info$Value[13])
  Val$ValidationYearMean   <- as.numeric(T_info$Value[14])
  Val$ValidationYearLower  <- as.numeric(T_info$Value[15])
  
  # LayerDating and Missing Inventories left as NA
  Val$LayerDating            <- NA
  Val$MissingInventoryTop    <- NA
  Val$MissingInventoryBottom <- NA
  Unc$MissingInventoryTop    <- NA
  Unc$MissingInventoryBottom <- NA
  
  Val$MonteCarloRuns <- 2e5
  
  #### Fill profile data ####
  Val$SampleCode         <- T_profiles$SampleCode
  Val$DepthLayerBottom   <- T_profiles$DepthLayerBottom
  Unc$DepthLayerBottom   <- T_profiles$DepthLayerU
  if (anyNA(Val$DepthLayerBottom) | anyNA(Unc$DepthLayerBottom)) 
    stop("All section bottom layers and uncertainties must be known. Please, revise. \n")
  if (!all(diff(Val$DepthLayerBottom) > 0)) 
    stop("Section bottom layers must be ordered. Please, read the manual. \n")
  Val$MassSection        <- T_profiles$MassSection
  Unc$MassSection        <- T_profiles$MassSectionU
  if (anyNA(Val$MassSection) | anyNA(Unc$MassSection)) 
    stop("All section masses and uncertainties must be known. Please, revise. \n")

  # Radionuclides
  Val$Ra226      <- T_profiles$Ra226
  Unc$Ra226      <- T_profiles$Ra226U
  if (sum(!is.na(Val$Ra226)) == 0) 
    if (verbose) cat("   - I have not seen any Ra-226 data. You may need to use ConstantRa \n")
  Val$Po210Total    <- T_profiles$Po210Total
  Unc$Po210Total    <- T_profiles$Po210TotalU
  Val$Pb210Total    <- T_profiles$Pb210Total
  Unc$Pb210Total    <- T_profiles$Pb210TotalU
  if (is.na(Val$Pb210Total[1]) & is.na(Val$Po210Total[1])) 
    stop("I need a value of Pb-210 or Po-210 in the surface. Please, revise. \n")
  
  Val$DatePlating           <- T_profiles$DatePlating
  Val$MeasurementAlphaStart <- T_profiles$MeasurementAlphaStart
  Val$MeasurementAlphaEnd   <- T_profiles$MeasurementAlphaEnd
  
  Val$MeasurementGammaStart <- T_profiles$MeasurementGammaStart
  Val$MeasurementGammaEnd   <- T_profiles$MeasurementGammaEnd
  
  # Auto-detect proxy
  has_proxy <- any(!is.na(T_profiles$ProxyValue) & T_profiles$ProxyValue != "")
  if (has_proxy) {
    Val$ProxyValue <- T_profiles$ProxyValue
    Unc$ProxyValue <- T_profiles$ProxyU
    Val$ProxyLabel <- ProxyLabel
  } else {
    Val$ProxyValue <- NA
    Unc$ProxyValue <- NA
    Val$ProxyLabel <- NA
  }
  
  #### New fields
  # These will be populated by downstream functions
  N_vec <- rep(NA_real_, N)
  
  # Ra-226
  Val$Ra226Complete <- N_vec
  
  # Po-210
  Val$Po210Excess              <- N_vec
  Val$Po210ExcessDecay         <- N_vec
  
  # Pb-210
  Val$Pb210Excess              <- N_vec
  Val$Pb210ExcessDecay         <- N_vec
  Val$Pb210ExcessDecayComplete <- N_vec
  
  # Uncertainties (same structure)
  Unc$Ra226Complete    <- N_vec
  
  Unc$Po210Excess      <- N_vec
  Unc$Po210ExcessDecay <- N_vec
  
  Unc$Pb210Excess              <- N_vec
  Unc$Pb210ExcessDecay         <- N_vec
  Unc$Pb210ExcessDecayComplete <- N_vec
  
  #### Finalize ####
  Val$NSections <- length(Val$DepthLayerBottom)
  if (verbose) {
    cat("ReadData completed successfully\n")
    cat("   - Profiles             :", Val$NSections, "sections. \n")
    cat("   - Validation           : 1 point. \n")
  }
  if (has_proxy)
    cat("   - Proxy detected (", Val$ProxyLabel, ").\n")
  else
    cat("   - No proxy data found. \n")    

  Val <- Val[order(names(Val))]
  Unc <- Unc[order(names(Unc))]
  List <- list(Val = Val, Unc = Unc)
  invisible(List)
}
