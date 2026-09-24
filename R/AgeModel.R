#' @title Pb-210 age model and validation
#'
#' @description
#' This function plots the age model after \eqn{^{210}}Pb dating and, if provided, 
#' the validation information (e.g., independent chronological markers). 
#' It also provides the statistical significance of the validation. 
#'
#' @importFrom graphics barplot polygon text
#' 
#' @param List Dating object, created by \code{\link{ReadData}} and transformed 
#' with the package functions.
#' @param PlotName Name of the plot (without extension) to be produced. 
#' The extension will automatically be set to \code{.png}. Default is \code{"Age Model"}.
#' @param DepthMin Minimum depth (cm) of the age model plot. Default is \code{NA}.
#' @param DepthMax Maximum depth (cm) of the age model plot. Default is \code{NA}.
#' @param TimeMin Minimum time (decimal age) of the age model plot. Default is \code{NA}.
#' @param TimeMax Maximum time (decimal age) of the age model plot. Default is \code{NA}.
#' @param SavePNG Logical. If TRUE (default), saves a high-resolution PNG file.
#' @param verbose Logical. If TRUE, prints progress logs to the console.
#' @param ... Extra arguments passed to plotting functions (e.g., \code{\link[graphics]{par}}).
#' 
#' @details
#' The age model is a representation of age (year CE) versus depth. 
#' For each point, both uncertainty bars are plotted. 
#' Minimum and maximum values for plotting are automatically calculated if not 
#' explicitly provided. 
#'
#' For validation, two approaches are used. 
#' For visual validation, 1-\eqn{\sigma} (68\%, red) and 2-\eqn{\sigma} (95\%, pink) 
#' ellipses of each validation point are plotted. 
#' The validation ellipses are calculated assuming that the 
#' variables (year and depth of the event) follow a uniform (square) distribution. 
#'
#' The validation ellipses include the user-provided label. 
#' As in most cases these will correspond to radionuclides, numbers are superscript,
#' so we recommend using the format "123Element" (e.g., 137Cs). 
#'
#' For statistical validation, a significance (p-value) of the Mahalanobis distance between the validation point and each age-model curve, is also provided. 
#' 
#' @returns Generates a \code{.PNG} plot of the age model with validation ellipses
#' and p-values. 
#' The plot is also reproduced in the R graphics device (e.g., RStudio plots pane).
#'
#' @examples
#' \dontrun{
#' a <- ReadData()
#' a <- Equilibrium(a, interactive = FALSE)
#' a.cf <- Pb210CF(a)
#' AgeModel(a.cf)}
#'
#' @seealso \code{\link{Pb210CF}}, \code{\link{Pb210CFCS}}
#'
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

AgeModel <- function(List = List, PlotName = "Age_Model", 
                     DepthMin = NA, DepthMax = NA, TimeMin = 1900, TimeMax = NA, 
                     SavePNG = TRUE, verbose = TRUE, ...){
  
  Val <- List$Val
  Unc <- List$Unc
  
  PlotName1 <- paste0(PlotName, ".png")
  if (verbose) cat("Age model plot:\n")
  
  #### Detect models ####
  hasCF   <- "YearSectionMean" %in% names(Val) && any(!is.na(Val$YearSectionMean))
  hasCFCS <- "YearCFCS"        %in% names(Val) && any(!is.na(Val$YearCFCS))
  
  if (verbose) {
    if (hasCF)   cat("   - CF age model detected.\n")
    if (hasCFCS) cat("   - CFCS age model detected.\n")
  }
  
  #### Get LayerDating ###
  LayerDating <- if ("LayerDating" %in% names(Val) && !is.na(Val$LayerDating)) {
    Val$LayerDating
  } else {
    tail(Val$DepthLayerBottom, 1)
  }
  if (verbose) cat("   - Equilibrium layer (LayerDating):", LayerDating, "cm\n")
  
  # Get Validation
  NValidation <- length(Val$ValidationDepthMean)
  DepthRef    <- Val$ValidationDepthMean
  DepthLower  <- Val$ValidationDepthLower
  DepthUpper  <- Val$ValidationDepthUpper
  YearRef     <- Val$ValidationYearMean
  YearLower   <- Val$ValidationYearLower
  YearUpper   <- Val$ValidationYearUpper
  
  # Label isotope expressions for validation text
  LabelExpr <- IsotopeLegend(Val$ValidationLabel)
  
  ## Extract CF data, ONLY up to LayerDating
  if (hasCF) {
    # number of dated sections
    n_cf      <- length(Val$YearSectionMean)
    depth_all <- Val$DepthSectionMean[1:n_cf]
    # select sections
    selCF      <- !is.na(Val$YearSectionMean) & is.finite(Val$YearSectionMean) & depth_all <= (LayerDating + 1e-6)
    YearCF     <- Val$YearSectionMean[selCF]
    DepthCF    <- depth_all[selCF]
    YearCFU    <- Unc$YearSectionMean[selCF]
    DepthCFU   <- Unc$DepthSectionMean[selCF]
    # Layer extent (top and bottom of each section) for vertical bars
    layer_bot  <- Val$DepthLayerBottom
    layer_top  <- c(0, layer_bot[-length(layer_bot)])
    DepthCFBot <- layer_bot[selCF]
    DepthCFTop <- layer_top[selCF]
    
    if (verbose) cat("   - CF points plotted:", length(YearCF), "\n")
  } else {
    YearCF <- DepthCF <- YearCFU <- DepthCFU <- DepthCFBot <- DepthCFTop <- numeric(0)
  }
  
  ## Extract CFCS data, ONLY up to LayerDating
  if (hasCFCS) {
    # number of dated sections
    n_cfcs    <- length(Val$YearCFCS)
    depth_all <- Val$DepthSectionMean[1:n_cfcs]
    # select sections
    selCFCS      <- !is.na(Val$YearCFCS) & is.finite(Val$YearCFCS) & depth_all <= (LayerDating + 1e-6)
    YearCFCS     <- Val$YearCFCS[selCFCS]
    DepthCFCS    <- Val$DepthSectionMean[selCFCS]
    YearCFCSU    <- Unc$YearCFCS[selCFCS]
    DepthCFCSU   <- Unc$DepthSectionMean[selCFCS]
    # Layer extent (top and bottom of each section) for vertical bars
    layer_bot    <- Val$DepthLayerBottom
    layer_top    <- c(0, layer_bot[-length(layer_bot)])
    DepthCFCSBot <- layer_bot[selCFCS]
    DepthCFCSTop <- layer_top[selCFCS]
    if (verbose) cat("   - CFCS points plotted:", length(YearCFCS), "\n")
  } else {
    YearCFCS <- DepthCFCS <- YearCFCSU <- DepthCFCSU <- DepthCFCSBot <- DepthCFCSTop <- numeric(0)
  }
  
  #### Plot limits ####
  all_depths <- c(DepthCF, DepthCFCS)
  all_years  <- c(YearCF, YearCFCS)
  
  if (is.na(DepthMax)) DepthMax <- ceiling(max(all_depths, na.rm = TRUE) / 5) * 5
  if (is.na(DepthMin)) DepthMin <- 0
  YearMax <- max(all_years, na.rm = TRUE)
  YearMin <- min(all_years, na.rm = TRUE)
  if (is.na(TimeMax)) TimeMax <- ceiling(YearMax / 10) * 10
  if (is.na(TimeMin)) TimeMin <- floor(YearMin / 10) * 10
  
  #### Statistical validation ####
  # p-values
  CF.p   <- NA
  CFCS.p <- NA
  
  if (NValidation > 0) {
    ## Validation uncertainties (equivalent 1σ of uniform distributions)
    SyVal <- (YearUpper[1]  - YearLower[1])  / (2 * sqrt(3))
    SdVal <- (DepthLower[1] - DepthUpper[1]) / (2 * sqrt(3))
    
    #### Internal function ####
    .curve.test <- function(Year, Depth, YearU, DepthU) {
      M2.min <- Inf
      for(i in 1:(length(Year)-1)) {
        dy <- Year [i+1] - Year [i]
        dd <- Depth[i+1] - Depth[i]
        ## Average uncertainties over the segment
        Sy <- sqrt(((YearU [i] + YearU [i+1])/2)^2 + SyVal^2)
        Sd <- sqrt(((DepthU[i] + DepthU[i+1])/2)^2 + SdVal^2)
        ## Work in normalized coordinates
        x1 <- Year[i]     / Sy
        y1 <- Depth[i]    / Sd
        x2 <- Year[i+1]   / Sy
        y2 <- Depth[i+1]  / Sd
        xv <- YearRef[1]  / Sy
        yv <- DepthRef[1] / Sd
        vx <- x2 - x1
        vy <- y2 - y1
        t <- ((xv-x1)*vx + (yv-y1)*vy)/(vx^2 + vy^2)
        t <- max(0,min(1,t))
        xc <- x1 + t*vx
        yc <- y1 + t*vy
        # statistcs
        M2 <- (xv-xc)^2 + (yv-yc)^2
        if(M2 < M2.min) M2.min <- M2
      }
      # probability
      p <- 1-stats::pchisq(M2.min,df=2)
      list(M2 = M2.min, p = p)
    }
    
    #### CF ####
    if(hasCF && length(YearCF)>1){
      A <- .curve.test(YearCF, DepthCF, YearCFU, DepthCFU)
      CF.p  <- A$p
      CF.M2 <- A$M2
    }
    
    #### CFCS ####
    if(hasCFCS && length(YearCFCS)>1){
      A <- .curve.test(YearCFCS, DepthCFCS, YearCFCSU, DepthCFCSU)
      CFCS.p  <- A$p
      CFCS.M2 <- A$M2
    }
  }
  
  #### Drawing function ####
  .draw_age_model <- function() {
    plot(TimeMin, DepthMin, xlim = c(TimeMin, TimeMax), ylim = c(DepthMax, DepthMin),
         col = "white", pch = 20, xaxs = "i", yaxs = "i",
         axes = FALSE, ann = FALSE)
    
    # Validation ellipses
    if (NValidation > 0) {
      AreaPlot   <- (TimeMax - TimeMin) * (DepthMax - DepthMin)
      AreaSquare <- (YearUpper - YearLower) * (DepthLower - DepthUpper)
      Intensity  <- pmax(0, 1 - AreaSquare / AreaPlot)
      
      for (i in seq_len(NValidation)) {
        EllipseFull(x0 = YearRef[i], y0 = DepthRef[i],
                    a = abs(YearUpper[i] - YearLower[i]),
                    b = abs(DepthLower[i] - DepthUpper[i]),
                    col = "red", intensity = Intensity[i] / 4)
        EllipseFull(x0 = YearRef[i], y0 = DepthRef[i],
                    a = 0.5 * abs(YearUpper[i] - YearLower[i]),
                    b = 0.5 * abs(DepthLower[i] - DepthUpper[i]),
                    col = "red", intensity = Intensity[i])
      }
      points(YearRef, DepthRef, col = "black", pch = 20)
      
      for (i in seq_along(LabelExpr)) {
        text(x = YearRef[i], y = DepthRef[i], 
             labels = LabelExpr[i], pos = 2, cex = 0.75)
      }
    }
    
    # CF model (blue) - now correctly limited
    if (hasCF && length(YearCF) > 0) {
      points(YearCF, DepthCF, col = "blue", pch = 20, cex = 1.2)
      # age uncertainties
      arrows(YearCF, DepthCFTop, YearCF, DepthCFBot,
             angle=90, code=3, length=0.03, col="blue")
      arrows(YearCF - YearCFU, DepthCF, YearCF + YearCFU, DepthCF,
             angle=90, code=3, length=0.03, col="blue")
    }
    
    # CFCS model (black)
    if (hasCFCS && length(YearCFCS) > 0) {
      points(YearCFCS, DepthCFCS, col = "black", pch = 20)
      arrows(YearCFCS, DepthCFCSTop, YearCFCS, DepthCFCSBot,
             angle=90, code=3, length=0.03, col="black")
      arrows(YearCFCS - YearCFCSU, DepthCFCS, YearCFCS + YearCFCSU, DepthCFCS,
             angle=90, code=3, length=0.03, col="black")
    }
    
    # axes
    axis(2, las = 1, cex = 1)
    axis(3)
    box()
    mtext("Depth (cm)", side = 2, line = 2.5, cex = 1.25)
    mtext("Year", side = 3, line = 2.5, cex = 1.25)
    
    # Legend
    legend.text <- character(0)
    legend.col  <- character(0)
    if(hasCF){
      legend.text <- c(legend.text, paste0("CF (", ifelse(CF.p<0.05,"p<0.05","p>0.05"), ")"))
      legend.col <- c(legend.col,"blue")
    }
    
    if(hasCFCS){
      legend.text <- c(legend.text, paste0("CFCS (", ifelse(CFCS.p<0.05,"p<0.05","p>0.05"), ")"))
      legend.col <- c(legend.col,"black")
    }
    
    legend("topleft", legend=legend.text, col=legend.col, pch=20, bty="n")
  }
  
  # Generate PNG + screen plot
  oldw <- getOption("warn")
  options(warn = -1)
  
  if (SavePNG) {
    png(PlotName1, res = 300, width = 15, height = 15, units = "cm")
    par(mar = c(1, 4, 4, 1))
    .draw_age_model()
    dev.off()
  }
  
  # plots pane
  par(mar = c(1, 4, 4, 1))
  .draw_age_model()
  
  options(warn = oldw)
  
  # return statistics
  return(c(CF.p = CF.p, CF.M2 = CF.M2, CFCS.p = CFCS.p, CFCS.M2 = CFCS.M2))
  invisible(NULL)
}
