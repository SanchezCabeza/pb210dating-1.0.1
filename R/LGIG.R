#' @title LGIG internal workflow: combine Po-210 (alpha) and Pb-210 (gamma) measurements
#'
#' @description
#' Operational LGIG data preparation, not for external use. 
#' Runs \code{\link{ReadData}}, \code{\link{DecayCorrection}}, and — only
#' when \eqn{^{210}}Po alpha-spectrometry data are present — combines excess
#' \eqn{^{210}}Po with excess \eqn{^{210}}Pb via \code{lm(Pb ~ Po)}, so that
#' Po-only sections are expressed on the \eqn{^{210}}Pb scale.
#'
#' When no Po data are present the function skips the combination step
#' entirely and returns the decay-corrected list ready for dating.
#' 
#' In both cases, CompleteProfiles is performed and 
#' a CSV reporting table of decay-corrected total activities
#' is written to the working directory.
#'
#' @param NameFile Name of the .CSV file with core data (see \code{\link{ReadData}})
#'   Default \code{NA} loads the built-in TehuaII example.
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#' 
#' @returns Updated dating list with:
#' \itemize{
#'   \item \code{Val$Pb210ExcessDecay}: decay-corrected excess \eqn{^{210}}Pb 
#'     (combined with \eqn{^{210}}Po where available).
#'   \item \code{Val$Pb210TotalCorrected} / \code{Unc$Pb210TotalCorrected}:
#'     decay-corrected total \eqn{^{210}}Pb for reporting.
#'   \item \code{Val$Po210TotalCorrected} / \code{Unc$Po210TotalCorrected}:
#'     decay-corrected total \eqn{^{210}}Po for reporting (when Po data exist).
#'   \item \code{Val$Pb210Source}: origin of each combined value
#'     (\code{"Pb-210"} or \code{"lm(Pb~Po)"}); \code{NA} when no Po data.
#' }
#' Side effects: PNG profile plot, optional regression + combined profile
#' plots, and a \code{"<Core> corrected.csv"} reporting file.
#'
#' @seealso \code{\link{ReadData}}, \code{\link{DecayCorrection}}, \code{\link{CompleteProfile}}. 
#'
#' @keywords internal
#' @noRd

LGIG <- function(NameFile = NA, verbose = TRUE) {
  
  if (verbose) cat("Start LGIG workflow. \n")
  
  # ------------------------------------------------------------------
  # Read data & initial plot
  # ------------------------------------------------------------------
  List <- ReadData(NameFile)
  core     <- if (!is.na(List$Val$Core)) List$Val$Core else "Core"
  hasPb    <- !all(is.na(List$Val$Pb210Total))
  hasPo    <- !all(is.na(List$Val$Po210Total))
  hasProxy <- any(!is.na(List$Val$ProxyValue))
  
  # Plot raw input profiles (total activities only — excess not yet calculated)
  PlotProfiles(List, NamePlot = paste0(core, " profiles"), Proxy = TRUE, verbose = verbose)
  
  # ------------------------------------------------------------------
  # Decay-correction
  # ------------------------------------------------------------------
  
  List <- DecayCorrection(List)
  
  Pb     <- List$Val$Pb210ExcessDecay
  PbU    <- List$Unc$Pb210ExcessDecay
  Po     <- List$Val$Po210ExcessDecay
  PoU    <- List$Unc$Po210ExcessDecay
  Ra     <- List$Val$Ra226
  RaU    <- List$Unc$Ra226
  Proxy  <- List$Val$ProxyValue
  ProxyU <- List$Unc$ProxyValue
  
  # if (all(is.na(Pb)))
  #   stop("   - Pb210ExcessDecay is all NA after DecayCorrection(). Check input file.\n")
  # 
  cat("Start combination of Po and Pb data: \n")
  
  # ------------------------------------------------------------------
  # a. No Po data — skip combination
  # ------------------------------------------------------------------
  if (!hasPo) {
    if (verbose) cat("   - No Po-210 data found. Combination step skipped.\n")
    
    List$Val$Pb210Source <- rep("Pb", List$Val$NSections)
    
    ## Reporting table (Pb only)
    Pb_corr  <- List$Val$Ra226 + Pb
    PbU_corr <- PbU
    
    report_df <- data.frame(
      SampleCode          = List$Val$SampleCode,
      LayerBottom_cm      = round(List$Val$DepthLayerBottom, 1),
      Pb210TotalCorr_Bqkg = round(Pb_corr,  2),
      u_Pb210TotalCorr    = round(PbU_corr, 2),
      Ra226_Bqkg          = round(Ra, 2),
      u_Ra226             = round(RaU, 2)
    )
    if (hasProxy) {
      report_df$Proxy   <- round(Proxy,  2)
      report_df$u_Proxy <- round(ProxyU, 2)
    }
    
    List$Val$Pb210TotalCorrected <- Pb_corr
    List$Unc$Pb210TotalCorrected <- PbU_corr
    
    utils::write.csv(report_df, paste0(core, " corrected.csv"),
              row.names = FALSE, quote = FALSE, na = "")
    
    if (verbose) {
      cat("\n   --- Decay-corrected total activities (to sampling date) ---\n")
      print(report_df, row.names = FALSE)
      cat("   -----------------------------------------------------------\n\n")
      cat("   - Reporting file written:", paste0(core, " corrected.csv"), "\n")
      cat("   - LGIG finished (Pb-210 only).\n")
      cat("   - Proceed: Equilibrium() -> Pb210CF() -> AgeModel() -> Output()\n")
    }
    List <- CompleteProfile(List, verbose = verbose)
    List$Val <- List$Val[order(names(List$Val))]
    List$Unc <- List$Unc[order(names(List$Unc))]
    return(invisible(List))
  }
  
  # ------------------------------------------------------------------
  # a. No Pb data — skip combination
  # ------------------------------------------------------------------
  if (!hasPb) {
    if (verbose) cat("   - No Pb-210 data found. Combination step skipped.\n")
    
    List$Val$Pb210Source <- rep("Po", List$Val$NSections)
    
    # Reporting table (Po only)
    Po_corr  <- List$Val$Ra226 + Po
    PoU_corr <- PoU
    
    report_df <- data.frame(
      SampleCode          = List$Val$SampleCode,
      LayerBottom_cm      = round(List$Val$DepthLayerBottom, 1),
      Po210TotalCorr_Bqkg = round(Po_corr,  2),
      u_Po210TotalCorr    = round(PoU_corr, 2),
      Ra226_Bqkg          = round(Ra, 2),
      u_Ra226             = round(RaU, 2)
    )
    if (hasProxy) {
      report_df$Proxy   <- round(Proxy,  2)
      report_df$u_Proxy <- round(ProxyU, 2)
    }
    
    List$Val$Po210TotalCorrected <- Po_corr
    List$Unc$Po210TotalCorrected <- PoU_corr
    
    # For dating, copy Po210 to Pb210
    List$Val$Pb210TotalCorrected <- Po_corr
    List$Unc$Pb210TotalCorrected <- PoU_corr
    
    utils::write.csv(report_df, paste0(core, " corrected.csv"),
              row.names = FALSE, quote = FALSE, na = "")
    
    if (verbose) {
      cat("\n   --- Decay-corrected total activities (to sampling date) ---\n")
      print(report_df, row.names = FALSE)
      cat("   -----------------------------------------------------------\n\n")
      cat("   - Reporting file written:", paste0(core, " corrected.csv"), "\n")
      cat("   - LGIG finished (Pb-210 only).\n")
      cat("   - Proceed: Equilibrium() -> Pb210CF() -> AgeModel() -> Output()\n")
    }
    
    List <- CompleteProfile(List, verbose = verbose)
    List$Val <- List$Val[order(names(List$Val))]
    List$Unc <- List$Unc[order(names(List$Unc))]
    return(invisible(List))
  }
  
  # ------------------------------------------------------------------
  # b. Po and Pb data present — combine
  # ------------------------------------------------------------------
  
  # Paired sections & regressions
  paired <- !is.na(Pb) & !is.na(Po)
  if (sum(paired) < 3) stop("   - At least 3 paired sections are required.")
  if (verbose) cat("   - Paired sections:", sum(paired), "\n")
  df_paired <- data.frame(Pb = Pb[paired], PbU = PbU[paired],
                          Po = Po[paired], PoU = PoU[paired])
  
  # lm(Pb ~ Po): Po is predictor, Pb is response — forward direction only
  fit_PbfromPo <- lm(Pb ~ Po, data = df_paired)
  
  # Forward predictionhelper  with propagated uncertainty
  .fwd <- function(fit, x) {
    p <- predict(fit,
                 newdata = data.frame(setNames(list(x), names(fit$model)[2])),
                 se.fit = TRUE)
    list(val = p$fit, unc = sqrt(summary(fit)$sigma^2 + p$se.fit^2))
  }
  
  # Combination: keep measured Pb; fill Po-only sections via lm(Pb ~ Po)
  N            <- List$Val$NSections
  combined     <- Pb
  combined_unc <- PbU
  source_type  <- rep(NA_character_, N)
  source_type[!is.na(Pb)] <- "Pb-210"
  
  idx_po_only <- is.na(Pb) & !is.na(Po)
  if (any(idx_po_only)) {
    est <- .fwd(fit_PbfromPo, Po[idx_po_only])
    combined[idx_po_only]     <- est$val
    combined_unc[idx_po_only] <- sqrt(est$unc^2 + PoU[idx_po_only]^2)
    PoU[idx_po_only]          <- combined_unc[idx_po_only]
    source_type[idx_po_only]  <- "lm(Pb~Po)"
    if (verbose)
      cat("   -", sum(idx_po_only), "Po-only section(s) filled via lm(Pb~Po).\n")
  }
  
  # Plots
  lab_pb  <- expression(""^{210}*"Pb (Bq kg"^{-1}*")")
  lab_po  <- expression(""^{210}*"Po (Bq kg"^{-1}*")")
  lab_act <- expression("Activity (Bq kg"^{-1}*")")
  
  # Regression plot: x = Po, y = Pb
  #.plot_reg <- function(core, fit, x, y, xU, yU, x_expr, y_expr) {
  .plot_reg <- function() {
    s     <- summary(fit_PbfromPo)
    x     <- df_paired$Po;  xU <- df_paired$PoU
    y     <- df_paired$Pb;  yU <- df_paired$PbU
    max_v <- max(c(x + xU, y + yU), na.rm = TRUE) * 1.1
    fn    <- paste0(core, " regression.png")
    
    .draw <- function() {
      par(pty = "s", mar = c(5, 5, 2, 2))
      plot(x, y, pch = 20, xlim = c(0, max_v), ylim = c(0, max_v),
           xlab = lab_po, ylab = lab_pb, xaxs = "i", yaxs = "i", las = 1)
      abline(0, 1, lty = 2, col = "grey70")
      segments(x - xU, y, x + xU, y, col = "grey50")
      segments(x, y - yU, x, y + yU, col = "grey50")
      abline(fit_PbfromPo, col = "red3", lwd = 2)
      legend("topleft", bty = "n", cex = 1,
             legend = c(
               paste0("Slope: ", round(coef(fit_PbfromPo)[2], 3),
                      " \u00B1 ", round(s$coefficients[2, 2], 3)),
               paste0("r = ",    round(sqrt(s$r.squared), 3),
                      " | p = ", round(s$coefficients[2, 4], 4))))
    }
    png(fn, res = 300, width = 12, height = 12, units = "cm")
    .draw(); dev.off()
    .draw()
    if (verbose) cat("   - Regression plot saved:", fn, "\n")
  }
  
  # Profile plots: "raw" (Pb vs Po) and "combined"
  .plot_profs <- function(pmode = "raw") {
    DepthL <- List$Val$DepthLayerBottom
    Depth  <- (c(0, DepthL[-length(DepthL)]) + DepthL) / 2
    fn     <- paste0(core, " ", pmode, ".png")
    
    .draw <- function() {
      par(mar = c(4, 4, 4, 1), pty = "m")
      if (pmode == "raw") {
        mx <- max(c(Pb + PbU, Po + PoU), na.rm = TRUE) * 1.1
        plot(0, 0, ylim = c(max(DepthL), 0), xlim = c(0, mx),
             axes = FALSE, ann = FALSE, type = "n", xaxs = "i", yaxs = "i")
        arrows(Pb - PbU, Depth, Pb + PbU, Depth,
               angle = 90, code = 3, length = 0.03, col = "red3")
        points(Pb, Depth, col = "red3", pch = 20)
        arrows(Po - PoU, Depth, Po + PoU, Depth,
               angle = 90, code = 3, length = 0.03, col = "steelblue")
        points(Po, Depth, col = "steelblue", pch = 20)
        legend("bottomright", bty = "n", pch = 20,
               col = c("red3", "steelblue"),
               legend = c(expression(""^{210}*"Pb"),
                          expression(""^{210}*"Po")))
      } else {
        mx <- max(combined + combined_unc, na.rm = TRUE) * 1.1
        col_map <- c("Pb-210"     = "red3",
                     "lm(Pb~Po)" = "steelblue")
        col_pts <- col_map[source_type]
        col_pts[is.na(col_pts)] <- "grey"
        plot(combined, Depth, ylim = c(max(DepthL), 0), xlim = c(0, mx),
             axes = FALSE, ann = FALSE, col = col_pts, pch = 20,
             xaxs = "i", yaxs = "i")
        arrows(combined - combined_unc, Depth,
               combined + combined_unc, Depth,
               angle = 90, code = 3, length = 0.03, col = col_pts)
        used <- unique(source_type[!is.na(source_type)])
        legend("bottomright", bty = "n", pch = 20, cex = 1,
               col = col_map[used], legend = used)
      }
      axis(2, las = 1); axis(3); box()
      mtext("Depth (cm)", 2, 2.5)
      mtext(lab_act, 3, 2.5)
    }
    png(fn, res = 300, width = 10, height = 15, units = "cm")
    .draw(); dev.off()
    .draw()
    if (verbose) cat("   - Profile plot saved:", fn, "\n")
  }
  
  .plot_reg()
  .plot_profs("raw")
  .plot_profs("combined")
  
  # ------------------------------------------------------------------
  # Reporting table
  # ------------------------------------------------------------------
  # Po expressed on Pb scale via lm(Pb ~ Po) for all measured Po sections
  Po_on_Pb_scale     <- rep(NA_real_, N)
  Po_on_Pb_scale_unc <- rep(NA_real_, N)
  idx_po <- !is.na(Po)
  which(idx_po)
  if (any(idx_po)) {
    est <- .fwd(fit_PbfromPo, Po[idx_po])
    Po_on_Pb_scale    [idx_po] <- est$val
    Po_on_Pb_scale_unc[idx_po] <- sqrt(est$unc^2 + PoU[idx_po]^2)
  }
  
  Pb_corr  <- List$Val$Ra226 + Pb
  PbU_corr <- PbU
  Po_corr  <- List$Val$Ra226 + Po_on_Pb_scale
  PoU_corr <- Po_on_Pb_scale_unc
  
  report_df <- data.frame(
    SampleCode          = List$Val$SampleCode,
    LayerBottom_cm      = round(List$Val$DepthLayerBottom, 1),
    Pb210TotalCorr_Bqkg = round(Pb_corr,  2),
    u_Pb210TotalCorr    = round(PbU_corr, 2),
    Po210TotalCorr_Bqkg = round(Po_corr,  2),
    u_Po210TotalCorr    = round(PoU_corr, 2),
    Ra226_Bqkg          = round(Ra, 2),
    u_Ra226             = round(RaU, 2)
  )
  
  if (hasProxy) {
    report_df$Proxy   <- round(List$Val$ProxyValue, 2)
    report_df$u_Proxy <- round(List$Unc$ProxyValue, 2)
  }
  
  utils::write.csv(report_df, paste0(core, " corrected.csv"),
            row.names = FALSE, quote = FALSE, na = "")
  if (verbose) cat("   - Reporting file written:", paste0(core, " corrected.csv"), "\n")
  
  # ------------------------------------------------------------------
  # Update List and finalise
  # ------------------------------------------------------------------
  List$Val$Pb210ExcessDecay    <- combined
  List$Unc$Pb210ExcessDecay    <- combined_unc
  List$Val$Pb210Source         <- source_type
  List$Val$Pb210TotalCorrected <- Pb_corr
  List$Unc$Pb210TotalCorrected <- PbU_corr
  List$Val$Po210TotalCorrected <- Po_corr
  List$Unc$Po210TotalCorrected <- PoU_corr
  
  List     <- CompleteProfile(List, verbose = verbose)
  List$Val <- List$Val[order(names(List$Val))]
  List$Unc <- List$Unc[order(names(List$Unc))]
  
  cat("   - LGIG finished. Pb210ExcessDecay now holds the combined values.\n")
  cat("   - Proceed: Equilibrium() -> Pb210CF() -> AgeModel() -> Output()\n")
  invisible(List)
}
