#' @title Output merged Pb-210 dating results (CF + CFCS) in a single CSV file
#'
#' @description
#' Creates one single CSV file with per-section data in
#' and core-wide constants in the last four columns (`Name`, `Value`,
#' `Uncertainty`, `p-value`).
#'
#' @param List Dating list created by \code{ReadData} and processed by
#'   \code{Pb210CF} and/or \code{Pb210CFCS}.
#' @param verbose Logical.Print progress messages. Default is \code{TRUE}.   
#'
#' @returns Invisibly returns the final data.frame. Side effect: writes
#'   \code{"<CoreName> results.csv"} in the working directory.
#' 
#' @examples
#' \dontrun{
#' Output(List = List)}
#'
#'@seealso \code{\link{ReadData}}, \code{\link{Pb210CF}}, 
#'\code{\link{Pb210CFCS}}, \code{\link{MissingInventory}}. 
#'
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @importFrom utils write.csv
#'
#' @export

Output <- function(List, verbose = TRUE) {
  
  Val <- List$Val
  Unc <- List$Unc
  
  # Define empty values
  CoreName    <- if (!is.null(Val$Core))   Val$Core   else "Core"
  AnalystName <- if (!is.null(Val$Dating)) Val$Dating else ""
  
  # ----------------------------------------------------------------------
  # Detect models
  # ----------------------------------------------------------------------
  hasCF   <- "FluxPb210" %in% names(Val) && any(!is.na(Val$YearSectionMean))
  hasCFCS <- "YearCFCS"         %in% names(Val) && any(!is.na(Val$YearCFCS))
  hasLGIG <- "Pb210TotalCorrected" %in% names(Val) && "Po210TotalCorrected" %in% names(Val)
  
  if (!hasCF && !hasCFCS)
    stop("No dating results (neither CF nor CFCS) found in the List.\n",
         "Run Pb210CF() and/or Pb210CFCS() first.")
  
  # ----------------------------------------------------------------------
  # Internal functions
  # ----------------------------------------------------------------------
  safe_vec <- function(v, len) {
    if (is.null(v) || length(v) == 0) return(rep(NA_real_, len))
    if (length(v) >= len) return(v[1:len])
    return(c(v, rep(NA_real_, len - length(v))))
  }
  
  # Coerce NULL or zero-length scalars to NA so as.data.frame() never
  # sees a list element of length 0 alongside elements of length 1.
  safe_scalar <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA)
    x
  }
  
  format_pval <- function(p_vec) {
    if (is.null(p_vec)) return(NA)
    res <- sapply(p_vec, function(p) {
      if (is.na(p)) return(NA)
      if (p < 0.001) return(formatC(p, format = "e", digits = 2))
      return(as.character(round(p, 4)))
    })
    paste(res, collapse = "; ")
  }
  
  # ----------------------------------------------------------------------
  # Mass-depth quantities
  # ----------------------------------------------------------------------
  n_sec <- length(Val$DepthSectionMean)
  if (!"MassDepthSection" %in% names(Val)) Val$MassDepthSection <- Val$MassSection / Val$CoreSurfaceArea
  MassDepthLayerCumulative   <- Val$MassDepthLayerCumulative
  MassDepthSectionCumulative <- Val$MassDepthSectionCumulative
  
  # ----------------------------------------------------------------------
  # Activity names & Proxy detection
  # ----------------------------------------------------------------------
  if (hasLGIG) {
    total_pb_name <- "Pb210TotalCorrected"
    total_po_name <- "Po210TotalCorrected"
    ra_name       <- "Ra226"
    ex_pb_name    <- "Pb210ExcessDecayComplete"
  } else {
    # Standard workflow
    total_pb_name <- "Pb210Total"
    total_po_name <- "Po210Total"
    ra_name       <- "Ra226"
    ex_pb_name    <- "Pb210ExcessDecayComplete"
  }

  total_pb     <- if (total_pb_name %in% names(Val)) Val[[total_pb_name]] else rep(NA_real_, n_sec)
  total_po     <- if (total_po_name %in% names(Val)) Val[[total_po_name]] else rep(NA_real_, n_sec)
  unc_total_pb <- if (total_pb_name %in% names(Unc)) Unc[[total_pb_name]] else rep(NA_real_, n_sec)
  unc_total_po <- if (total_po_name %in% names(Unc)) Unc[[total_po_name]] else rep(NA_real_, n_sec)
  ra           <- if (ra_name %in% names(Val)) Val[[ra_name]] else rep(NA_real_, n_sec)
  unc_ra       <- if (ra_name %in% names(Unc)) Unc[[ra_name]] else rep(NA_real_, n_sec)
  ex_pb        <- if (ex_pb_name %in% names(Val)) Val[[ex_pb_name]] else rep(NA_real_, n_sec)
  unc_ex_pb    <- if (ex_pb_name %in% names(Unc)) Unc[[ex_pb_name]] else rep(NA_real_, n_sec)
  
  # --- Proxy Logic ---
  proxy_name_str <- "Proxy"
  proxy_val <- Val$ProxyValue
  unc_proxy <- Unc$ProxyValue
  
  # ----------------------------------------------------------------------
  # SECTION DATA
  # ----------------------------------------------------------------------
  df_main <- data.frame(
    Code                                = safe_vec(Val$SampleCode,       n_sec),
    `Layer (cm)`                        = safe_vec(Val$DepthLayerBottom, n_sec),
    `Depth (cm)`                        = safe_vec(Val$DepthSectionMean, n_sec),
    `Density (g cm-3)`                  = safe_vec(Val$DensityDry,       n_sec),
    `u(Density)`                        = safe_vec(Unc$DensityDry,       n_sec),
    `Layer accumulated mass (g cm-2)`   = MassDepthLayerCumulative[-1]/ 1e4,
    `Section accumulated mass (g cm-2)` = MassDepthSectionCumulative  / 1e4,
    `210Po (Bq kg-1)`                   = total_po,
    `u(210Po)`                          = unc_total_po,
    `210Pb (Bq kg-1)`                   = total_pb,
    `u(210Pb)`                          = unc_total_pb,
    `226Ra (Bq kg-1)`                   = ra,
    `u(226Ra)`                          = unc_ra,
    Proxy                               = proxy_val,
    `u(Proxy)`                          = unc_proxy,
    `210Pbex (Bq kg-1)`                 = safe_vec(ex_pb    [1:Val$NumberOfDatedSections], n_sec),
    `u(210Pbex)`                        = safe_vec(unc_ex_pb[1:Val$NumberOfDatedSections], n_sec),
    check.names = FALSE # does not change text style
    
  )
  
  if (hasCF) {
    df_CF <- data.frame(
      `CF bottom layer age (yr)` = safe_vec(Val$AgeLayer[-1],             n_sec),
      `u(CF bottom layer age)`   = safe_vec(Unc$AgeLayer[-1],             n_sec),
      `Top CF year (yr)`         = safe_vec(Val$YearLayer,                n_sec),
      `u(top CF year)`           = safe_vec(Unc$AgeLayer,                 n_sec),
      `Mean CF year (yr)`        = safe_vec(Val$YearSectionMean,          n_sec),
      `u(mean CF year)`          = safe_vec(Unc$AgeSectionMean,           n_sec),
      `Bottom CF year (yr)`      = safe_vec(Val$YearLayer[-1],            n_sec),
      `u(bottom CF year)`        = safe_vec(Unc$AgeLayer [-1],            n_sec),
      `MAR (g cm-2 yr-1)`        = safe_vec(Val$MassAccumulationRate,     n_sec),
      `u(MAR)`                   = safe_vec(Unc$MassAccumulationRate,     n_sec),
      `SAR (cm yr-1)`            = safe_vec(Val$SedimentAccumulationRate, n_sec),
      `u(SAR)`                   = safe_vec(Unc$SedimentAccumulationRate, n_sec),
      check.names = FALSE
    )
    cf.cols <- names(df_CF)
    df_main <- cbind(df_main, df_CF)
  }
  
  if (hasCFCS) {
    df_CFCS <- data.frame(
      `Mean CFCS age`     = safe_vec(Val$AgeCFCS,  n_sec),
      `u(mean CFCS age)`  = safe_vec(Unc$AgeCFCS,  n_sec),
      `Mean CFCS year`    = safe_vec(Val$YearCFCS, n_sec),
      `u(mean CFCS year)` = safe_vec(Unc$AgeCFCS,  n_sec),
      check.names = FALSE
    )
    df_main <- cbind(df_main, df_CFCS)
  }
  
  # If there is NO bottom missing inventory, delete CF results below Val$NumberOfDatedSections
  is_zero_or_na <- function(x) {is.null(x) || length(x) == 0 || is.na(x[1]) || x[1] == 0}
  
  if (is_zero_or_na(Val$MissingInventoryBottom) && is_zero_or_na(Val$MissingInventoryBottomBq)) {
    df_main[Val$NumberOfDatedSections:n_sec, cf.cols] <- NA
  }
  
  # Search common proxy names
  if (!is.na(Val$ProxyLabel) && nzchar(Val$ProxyLabel)) {
    names(df_main)[names(df_main) == "Proxy"] <-
      IsotopeName(Val$ProxyLabel, units = "Bq kg-1")
    
    names(df_main)[names(df_main) == "u(Proxy)"] <-
      IsotopeUncertainty(Val$ProxyLabel)
  }
  
  # ----------------------------------------------------------------------
  # Format Numbers: Handle -Inf, Round to 2 decimals (Density/MAR/SAR to 4)
  # ----------------------------------------------------------------------
  for (col in names(df_main)) {
    if (is.numeric(df_main[[col]])) {
      df_main[[col]][is.infinite(df_main[[col]])] <- NA
      
      if (grepl("Density|u(Density)|MAR|SAR", col)) {
        df_main[[col]] <- round(df_main[[col]], 4)
      } else {
        df_main[[col]] <- round(df_main[[col]], 2)
      }
    }
  }
  
  # ----------------------------------------------------------------------
  # CONSTANTS BLOCK
  # ----------------------------------------------------------------------
  const_list <- list(
    list(Name = "Core name",          Value = safe_scalar(if (!is.null(Val$Core)) Val$Core else CoreName), Uncertainty = NA, pvalue = NA),
    list(Name = "Alpha Spectrometry", Value = safe_scalar(Val$AlphaSpectrometry), Uncertainty = NA, pvalue = NA),
    list(Name = "Gamma Spectrometry", Value = safe_scalar(Val$GammaSpectrometry), Uncertainty = NA, pvalue = NA),
    list(Name = "Dating",             Value = safe_scalar(if (!is.null(Val$Dating)) Val$Dating else AnalystName), Uncertainty = NA, pvalue = NA),
    list(Name = "Report date_time",   Value = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Uncertainty = NA, pvalue = NA),
    list(Name = "Number of CF runs",  Value = safe_scalar(if (hasCF) Val$MonteCarloRunsCF else NA), Uncertainty = NA, pvalue = NA),
    list(Name = "Number of CFCS runs", Value = safe_scalar(if (hasCFCS) Val$MonteCarloRunsCFCS else NA), Uncertainty = NA, pvalue = NA),
    list(Name = "Half life Pb-210 (yr)", Value = safe_scalar(Val$HalfLifePb210), Uncertainty = safe_scalar(Unc$HalfLifePb210), pvalue = NA),
    list(Name = "Half life Po-210 (yr)", Value = safe_scalar(Val$HalfLifePo210), Uncertainty = safe_scalar(Unc$HalfLifePo210), pvalue = NA),
    list(Name = "Sampling Date",      Value = safe_scalar(as.character(Val$DateSampling)), Uncertainty = NA, pvalue = NA),
    list(Name = "Internal diameter(cm)", Value = safe_scalar(round(Val$Diameter, 2)), Uncertainty = safe_scalar(round(Unc$Diameter, 2)), pvalue = NA),
    list(Name = "Number of sections", Value = safe_scalar(Val$NSections), Uncertainty = NA, pvalue = NA),
    list(Name = "Bottom layer for dating", Value = safe_scalar(Val$LayerDating), Uncertainty = NA, pvalue = NA),
    list(Name = "Missing inventory top (Bq m-2)",    Value = safe_scalar(round(Val$MissingInventoryTop)),    Uncertainty = safe_scalar(round(Unc$MissingInventoryTop)),    pvalue = NA),
    list(Name = "Missing inventory bottom (Bq m-2)", Value = safe_scalar(round(Val$MissingInventoryBottom)), Uncertainty = safe_scalar(round(Unc$MissingInventoryBottom)), pvalue = NA),
    list(Name = "210Pb inventory (Bq m-2)",
         Value       = safe_scalar(if ("TotalInventory" %in% names(Val)) round(Val$TotalInventory) else NA),
         Uncertainty = safe_scalar(if ("TotalInventory" %in% names(Unc)) round(Unc$TotalInventory) else NA),
         pvalue = NA),
    list(Name = "210Pb Flux (Bq m-2 yr-1)",
         Value       = safe_scalar(if ("FluxPb210" %in% names(Val)) round(Val$FluxPb210) else NA),
         Uncertainty = safe_scalar(if ("FluxPb210" %in% names(Unc)) round(Unc$FluxPb210) else NA),
         pvalue = NA)
  )
  
  if (hasCFCS) {
    const_list <- c(const_list, list(
      list(Name = "CFCS Breaks (cm)", Value = paste(round(Val$BreakPoints, 2), collapse = "; "), Uncertainty = NA, pvalue = NA),
      list(Name = "CFCS MAR (g cm-2 yr-1)", Value = paste(round(Val$MassAccumulationRateMean, 4), collapse = "; "),
           Uncertainty = paste(round(Unc$MassAccumulationRateMean, 4), collapse = "; "),
           pvalue = format_pval(Val$PValueMassAccumulationRate)),
      list(Name = "CFCS SAR (cm yr-1)", Value = paste(round(Val$SedimentAccumulationRateMean, 4), collapse = "; "),
           Uncertainty = paste(round(Unc$SedimentAccumulationRateMean, 4), collapse = "; "),
           pvalue = format_pval(Val$PValueSedimentAccumulationRate))
    ))
  }
  
  const_df <- do.call(rbind, lapply(const_list, as.data.frame, stringsAsFactors = FALSE))
  colnames(const_df) <- c("Name", "Value", "Uncertainty", "p-value")
  
  # ----------------------------------------------------------------------
  # Combine into final data.frame
  # ----------------------------------------------------------------------
  max_rows <- max(nrow(df_main), nrow(const_df))
  
  if (nrow(df_main) < max_rows) {
    pad_main <- data.frame(matrix(NA, nrow = max_rows - nrow(df_main), ncol = ncol(df_main)))
    colnames(pad_main) <- colnames(df_main)
    df_main <- rbind(df_main, pad_main)
  }
  
  if (nrow(const_df) < max_rows) {
    pad_const <- data.frame(matrix(NA, nrow = max_rows - nrow(const_df), ncol = ncol(const_df)))
    colnames(pad_const) <- colnames(const_df)
    const_df <- rbind(const_df, pad_const)
  }
  
  final_df <- cbind(df_main, const_df)
  
  # ----------------------------------------------------------------------
  # 9. Write CSV
  # ----------------------------------------------------------------------
  filename <- paste0(if (!is.null(Val$Core)) Val$Core else CoreName, "_results.csv")
  write.table(final_df, file = filename, sep = ",", na = "", row.names = FALSE,
              col.names = TRUE, quote = FALSE)
  
  cat("Output file created successfully: ", filename, "\n")
  
  # ----------------------------------------------------------------------
  # 10. Write full List .OUT
  # ----------------------------------------------------------------------
  # name
  filename.out <- paste0(if (!is.null(Val$Core)) Val$Core else CoreName, ".out")
  
  # number of sections
  N <- Val$NSections
  
  # Add "u_" prefix to uncertainties to prevent column name clashes
  names(Unc) <- paste0("u_", names(Unc))
  
  # 2. Combine all variables into one flat list
  combined_list <- c(Val, Unc)
  
  # Pad every element to exactly length N with NAs
  padded_list <- lapply(combined_list, function(x) {
    if (length(x) == 0)
      return(rep(NA, N))
    c(x, rep(NA, N))[1:N]
  })
  
  # write .out
  write.table(padded_list, file = filename.out, sep = ",", na = "", row.names = FALSE,
              col.names = TRUE, quote = FALSE)
  
  invisible(final_df)
}
