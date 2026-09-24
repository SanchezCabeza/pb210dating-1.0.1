#' @title Plots accumulation rates of sediment cores dated with the CF model
#'
#' @description
#' This function produces two independent PNG plots for the \strong{MassAccumulationRate} (MAR) 
#' and the \strong{SedimentAccumulationRate} (SAR) of cores previously dated with 
#' the Constant Flux (CF) model only. 
#' It uses the per-section rates calculated by \code{\link{Pb210CF}} (not the CFCS means).
#'
#' @param List Dating object, created by \code{\link{ReadData}} and dated with \code{\link{Pb210CF}}.
#' @param YearMin Minimum year to be plotted. Default is \code{1900}.
#' @param YearMax Maximum year to be plotted (auto-calculated if NA).
#' @param MARMin Minimum MAR to be plotted. Default is \code{0}.
#' @param MARMax Maximum MAR to be plotted (auto-calculated if NA).
#' @param SARMin Minimum SAR to be plotted. Default is \code{0}.
#' @param SARMax Maximum SAR to be plotted (auto-calculated if NA).
#' @param verbose Logical. If TRUE, prints progress logs (default TRUE).
#' @param ... Extra arguments passed to plotting functions.
#'
#' @returns Two PNG files (`Core MAR.png` and `Core SAR.png`) + combined plot in the RStudio pane.
#' 
#' @examples 
#' \dontrun{
#' After \code{\link{Pb210CF}} has been run: 
#' PlotMARSAR(List = List)}
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

PlotMARSAR <- function(List = List, YearMin = 1900, YearMax = NA, 
                       MARMin = 0, MARMax = NA, 
                       SARMin = 0, SARMax = NA, 
                       verbose = TRUE, ...){    
  
  if (verbose) cat("Plot CF variable MAR and SAR. \n")
  
  # Safety check: must have CF per-section rates
  if (!"MassAccumulationRate" %in% names(List$Val) || 
      !"SedimentAccumulationRate" %in% names(List$Val)) {
    stop("PlotMARSAR is designed exclusively for the CF model. The required variables MassAccumulationRate and SedimentAccumulationRate are missing. Run Pb210CF() first (or use PlotProfiles() / AgeModel() for CFCS).")
  }
  
  # Plot names
  Core <- List$Val$Core
  Name1 <- paste0(Core, "_MAR.png")
  Name2 <- paste0(Core, "_SAR.png")
  Name3 <- paste0(Core, "_MAR+SAR.png")
  
  # Read data using revised vocabulary (CF per-section rates)
  YearSection  <- List$Val$YearSectionMean
  YearSectionU <- List$Unc$YearSectionMean
  MARSection   <- List$Val$MassAccumulationRate
  SARSection   <- List$Val$SedimentAccumulationRate
  MARSectionU  <- List$Unc$MassAccumulationRate
  SARSectionU  <- List$Unc$SedimentAccumulationRate
  
  # Singularities
  MARSectionU[is.infinite(MARSection)] <- Inf
  SARSectionU[is.infinite(SARSection)] <- Inf
  
  # Limit to finite values
  Select <- is.finite(YearSection) & !is.na(MARSection) & !is.na(SARSection)
  if (sum(Select) == 0) stop("No finite accumulation rates found to plot.")
  
  YearSection   <- YearSection[Select]
  YearSectionU  <- YearSectionU[Select]
  MARSection    <- MARSection[Select]
  MARSectionU   <- MARSectionU[Select]
  SARSection    <- SARSection[Select]
  SARSectionU   <- SARSectionU[Select]
  
  # Calculate limits if not provided
  TimeMax <- max(YearSection[is.finite(YearSection)])
  TimeMin <- min(YearSection[is.finite(YearSection)])
  if (is.na(YearMax)) YearMax <- ceiling(TimeMax / 10) * 10
  if (is.na(YearMin)) YearMin <- floor(TimeMin / 10) * 10
  
  MarMax2 <- max(MARSection[is.finite(MARSection)] + MARSectionU[is.finite(MARSectionU)], na.rm = TRUE)
  if (is.na(MARMax)) MARMax <- ceiling(MarMax2 * 10) / 10
  if (is.na(MARMin)) MARMin <- 0
  
  SarMax2 <- max(SARSection[is.finite(SARSection)] + SARSectionU[is.finite(SARSectionU)], na.rm = TRUE)
  if (is.na(SARMax)) SARMax <- ceiling(SarMax2 * 10) / 10
  if (is.na(SARMin)) SARMin <- 0
  
  # Infinite points
  XInfinite <- YearSection[is.infinite(MARSection)]
  MARInfinite <- rep(MARMin, length(XInfinite))
  SARInfinite <- rep(SARMin, length(XInfinite))
  
  # Suppress arrow warnings
  oldw <- getOption("warn")
  options(warn = -1)
  
  # Plot MAR
  png(Name1, res = 300, width = 15, height = 10, units = "cm")
  par(mar = c(4,5,1,1))
  plot(YearSection, MARSection, xlim = c(YearMin, YearMax), ylim = c(MARMin, MARMax),
       col = "blue", pch = 20, xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
  arrows(x0 = YearSection, x1 = YearSection, y0 = MARSection - MARSectionU, y1 = MARSection + MARSectionU,
         angle = 90, code = 3, length = 0.03, col = "blue")
  arrows(y0 = MARSection, y1 = MARSection, x0 = YearSection - YearSectionU, x1 = YearSection + YearSectionU,
         angle = 90, code = 3, length = 0.03, col = "blue")
  axis(1); axis(2, las = 1, cex = 1); box()
  mtext("Year", 1, 2.5, cex = 1.25)
  mtext(expression(paste("MAR (g cm"^"-2", " yr"^"-1", ")")), 2, 2.5, cex = 1.25, padj = -0.5)
  points(XInfinite, MARInfinite, col = "red", pch = 16)
  dev.off()
  
  # Plot SAR
  png(Name2, res = 300, width = 15, height = 10, units = "cm")
  par(mar = c(4,5,1,1))
  plot(YearSection, SARSection, xlim = c(YearMin, YearMax), ylim = c(SARMin, SARMax),
       col = "blue", pch = 20, xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
  arrows(x0 = YearSection, x1 = YearSection, y0 = SARSection - SARSectionU, y1 = SARSection + SARSectionU,
         angle = 90, code = 3, length = 0.03, col = "blue")
  arrows(y0 = SARSection, y1 = SARSection, x0 = YearSection - YearSectionU, x1 = YearSection + YearSectionU,
         angle = 90, code = 3, length = 0.03, col = "blue")
  axis(1); axis(2, las = 1, cex = 1); box()
  mtext("Year", 1, 2.5, cex = 1.25)
  mtext(expression(paste("SAR (cm yr"^"-1", ")")), 2, 2.5, cex = 1.25, padj = -0.5)
  points(XInfinite, SARInfinite, col = "red", pch = 16)
  dev.off()
  
  # Combined MAR + SAR PNG
  png(Name3, res = 300, width = 15, height = 10, units = "cm")
  par(mar = c(4,5,1,5))
  plot(YearSection, MARSection, xlim = c(YearMin, YearMax), ylim = c(MARMin, MARMax),
       col = "blue", pch = 20, xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
  arrows(x0 = YearSection, x1 = YearSection, y0 = MARSection - MARSectionU, y1 = MARSection + MARSectionU,
         angle = 90, code = 3, length = 0.03, col = "blue")
  arrows(y0 = MARSection, y1 = MARSection, x0 = YearSection - YearSectionU, x1 = YearSection + YearSectionU,
         angle = 90, code = 3, length = 0.03, col = "blue")
  axis(1); axis(2, las = 1, cex = 1); box()
  mtext("Year", 1, 2.5, cex = 1.25)
  mtext(expression(paste("MAR (g cm"^"-2", " yr"^"-1", ")")), 2, 2.5, cex = 1.25, padj = -0.5)
  points(XInfinite, MARInfinite, col = "red", pch = 16)
  
  par(new = TRUE)
  plot(YearSection, SARSection, xlim = c(YearMin, YearMax), ylim = c(SARMin, SARMax),
       col = "black", pch = 20, xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
  arrows(x0 = YearSection, x1 = YearSection, y0 = SARSection - SARSectionU, y1 = SARSection + SARSectionU,
         angle = 90, code = 3, length = 0.03, col = "black")
  arrows(y0 = SARSection, y1 = SARSection, x0 = YearSection - YearSectionU, x1 = YearSection + YearSectionU,
         angle = 90, code = 3, length = 0.03, col = "black")
  axis(4, las = 1, cex = 1)
  mtext(expression(paste("SAR (cm yr"^"-1", ")")), 4, 2.5, cex = 1.25, padj = 0.5)
  points(XInfinite, SARInfinite, col = "red", pch = 16)
  legend("topleft", bty = "n", col = c("blue", "black"), pch = 20, legend = c("MAR", "SAR"))
  dev.off()
  
  
  # Combined RStudio pane
  par(mar = c(4,5,1,5))
  plot(YearSection, MARSection, xlim = c(YearMin, YearMax), ylim = c(MARMin, MARMax),
       col = "blue", pch = 20, xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
  arrows(x0 = YearSection, x1 = YearSection, y0 = MARSection - MARSectionU, y1 = MARSection + MARSectionU,
         angle = 90, code = 3, length = 0.03, col = "blue")
  arrows(y0 = MARSection, y1 = MARSection, x0 = YearSection - YearSectionU, x1 = YearSection + YearSectionU,
         angle = 90, code = 3, length = 0.03, col = "blue")
  axis(1); axis(2, las = 1, cex = 1); box()
  mtext("Year", 1, 2.5, cex = 1.25)
  mtext(expression(paste("MAR (g cm"^"-2", " yr"^"-1", ")")), 2, 2.5, cex = 1.25, padj = -0.5)
  points(XInfinite, MARInfinite, col = "red", pch = 16)
  
  par(new = TRUE)
  plot(YearSection, SARSection, xlim = c(YearMin, YearMax), ylim = c(SARMin, SARMax),
       col = "black", pch = 20, xaxs = "i", yaxs = "i", axes = FALSE, ann = FALSE)
  arrows(x0 = YearSection, x1 = YearSection, y0 = SARSection - SARSectionU, y1 = SARSection + SARSectionU,
         angle = 90, code = 3, length = 0.03, col = "black")
  arrows(y0 = SARSection, y1 = SARSection, x0 = YearSection - YearSectionU, x1 = YearSection + YearSectionU,
         angle = 90, code = 3, length = 0.03, col = "black")
  axis(4, las = 1, cex = 1)
  mtext(expression(paste("SAR (cm yr"^"-1", ")")), 4, 2.5, cex = 1.25, padj = 0.5)
  points(XInfinite, SARInfinite, col = "red", pch = 16)
  legend("topleft", bty = "n", col = c("blue", "black"), pch = 20, legend = c("MAR", "SAR"))
  
  options(warn = oldw)
}
