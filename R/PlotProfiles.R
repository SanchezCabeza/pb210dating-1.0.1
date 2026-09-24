#' @title Plot of sediment profiles: flexible selection of Pb-210, Po-210, Ra-226 and proxy
#' 
#' @description
#' Plots sediment core profiles of any combination of total/excess/complete \eqn{^{210}Pb}, 
#' total/excess/complete \eqn{^{210}Po}, total/complete \eqn{^{226}Ra} 
#' and one proxy (e.g. \eqn{^{137}Cs}). 
#' The user can select which variables to plot via the \code{PlotActivities} 
#' argument. If left \code{NULL}, the function checks the input data 
#' and plots total values. 
#' 
#' @param List Dating list created by \code{ReadData} and processed by other package functions.
#' @param PlotActivities Character vector specifying which activities to plot on the 
#'   main (left/top) axis. Possible values:  
#'   \code{"Pb210Total"}, \code{"Pb210Excess"}, 
#'   \code{"Pb210ExcessDecay}, \code{"Pb210ExcessDecayComplete} (this is the one used for dating), 
#'   \code{"Po210Total"}, \code{"Po210Excess}, \code{"Po210ExcessDecay}, 
#'   \code{"Po210ExcessDecayComplete}, 
#'   \code{"Ra226}, \code{"Ra226Complete}, \code{"Proxy"}. 
#'   The argument is \strong{case-insensitive} (e.g. "pb210excess", "PB210EXCESS", or "Pb210Excess" all work). 
#'   Default \code{NULL} automatically identifies and plots 
#'   \code{"Pb210Total"}, \code{"Po210Total"} and \code{"Ra226}.
#' @param DepthMin Minimum depth to plot (cm). If NA, estimated automatically.
#' @param DepthMax Maximum depth to plot (cm). If NA, estimated automatically.
#' @param PbMax Axis limits when \code{PlotActivities = NULL}). 
#' If not NULL, axis limit is computed automatically from the selected radionuclides.
#' @param RaMax Same as above. 
#' @param ProxyMax Same as above. 
#' @param Proxy Logical. Needed when \code{PlotActivities = NULL}).
#' @param NamePlot Base name for the PNG file (default "Profiles").
#' @param SavePNG Logical. If TRUE (default), saves a high-resolution PNG file.
#' @param verbose Logical. If TRUE, prints progress logs to the console (default \code{TRUE}).
#' @param ... Further arguments passed to graphical parameters.
#' 
#' @returns Invisibly returns NULL. Side effects: optionally creates a .png file and displays a plot.
#' 
#' @examples 
#' \dontrun{
#' List <- ReadData() # Reads demonstration core data (TehuaII.rda)
#' PlotProfiles(List = List, SavePNG = TRUE, PlotActivities = "Po210Total")}
#' 
#' @seealso \code{\link{ReadData}}, 
#' \code{\link{CompleteProfile}}, \code{\link{ConstantRa}}, 
#' \code{\link{DecayCorrection}}
#'
#' @references
#' Sanchez-Cabeza, J. A. & Ruiz-Fernandez, A. C. (2012). 
#' ^{210}Pb sediment radiochronology: an integrated formulation and 
#' classification of dating models. 
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' @export

PlotProfiles <- function(List, 
                         DepthMin = NA, DepthMax = NA,
                         PbMax    = NA, RaMax    = NA, 
                         ProxyMax = NA, Proxy    = FALSE,  
                         NamePlot = "Profiles", SavePNG  = TRUE, 
                         PlotActivities = NULL,   
                         verbose = TRUE,
                         ...) {
  
  if (verbose) cat("Start plotting profiles. \n")
  
  # ------------------------------------------------------------------
  # Read core geometry
  # ------------------------------------------------------------------
  DepthL <- List$Val$DepthLayerBottom
  Depth  <- RollMean2V(c(0, DepthL))
  
  # ------------------------------------------------------------------
  # Initialize data availability flags
  # ------------------------------------------------------------------
  hasPb <- hasPo <- hasRa <- hasProxy <- FALSE
  
  # Check what data is actually present in the List
  if (!is.null(List$Val$ProxyValue)) hasProxy <- any(!is.na(List$Val$ProxyValue))
  if (!is.null(List$Val$Pb210Total)) hasPb    <- any(!is.na(List$Val$Pb210Total))
  if (!is.null(List$Val$Po210Total)) hasPo    <- any(!is.na(List$Val$Po210Total))
  if (!is.null(List$Val$Ra226))      hasRa    <- any(!is.na(List$Val$Ra226))
  
  # ------------------------------------------------------------------
  # Determine what can be plot
  # ------------------------------------------------------------------
  key_map <- c(
    "pb210total"               = "Pb210Total",
    "pb210excess"              = "Pb210Excess",
    "pb210excessdecay"         = "Pb210ExcessDecay",
    "pb210excessdecaycomplete" = "Pb210ExcessDecayComplete",
    "po210excess"              = "Po210Excess",
    "po210total"               = "Po210Total",
    "po210excessdecay"         = "Po210ExcessDecay",
    "po210excessdecaycomplete" = "Po210ExcessDecayComplete",
    "ra226"                    = "Ra226",
    "ra226complete"            = "Ra226Complete"
    # Proxy removed from key_map to prevent double-plotting
  )
  
  if (!is.null(PlotActivities)) {
    user_input <- tolower(trimws(as.character(PlotActivities)))
    radionuclides <- key_map[user_input]
    radionuclides <- radionuclides[!is.na(radionuclides)]
    proxy_requested <- any(user_input == "proxy")
  } else {
    default_activities <- character(0)
    if (hasPb) default_activities <- c(default_activities, "Pb210Total")
    if (hasPo) default_activities <- c(default_activities, "Po210Total")
    if (hasRa) default_activities <- c(default_activities, "Ra226")
    
    radionuclides   <- default_activities
    proxy_requested <- Proxy || (hasProxy && length(default_activities) == 0) 
  }
  
  # ------------------------------------------------------------------
  # Explicitly legend labels 
  # ------------------------------------------------------------------
  act_info <- list(
    Pb210Total               = list(col = "#d7191c", lab = IsotopeLegend("210Pb", "Total")),
    Pb210Excess              = list(col = "#d7191c", lab = IsotopeLegend("210Pb", "Excess")),
    Pb210ExcessDecay         = list(col = "#d7191c", lab = IsotopeLegend("210Pb", "Excess-Decay")),
    Pb210ExcessDecayComplete = list(col = "#d7191c", lab = IsotopeLegend("210Pb", "Excess-Decay-Complete")),
    Po210Total               = list(col = "#fdae61", lab = IsotopeLegend("210Po", "Total")),
    Po210Excess              = list(col = "#fdae61", lab = IsotopeLegend("210Po", "Excess")),
    Po210ExcessDecay         = list(col = "#fdae61", lab = IsotopeLegend("210Po", "Excess-Decay")),
    Ra226                    = list(col = "#2c7bb6", lab = IsotopeLegend("226Ra")),
    Ra226Complete            = list(col = "#2c7bb6", lab = IsotopeLegend("226Ra", "Complete"))
  )
  
  # ------------------------------------------------------------------
  # Collect radionuclide data
  # ------------------------------------------------------------------
  plot_data <- list()
  max_act   <- 0
  
  for (vname in radionuclides) {
    if (vname %in% names(List$Val) && !all(is.na(List$Val[[vname]]))) {
      vv <- List$Val[[vname]]
      uu <- List$Unc[[vname]]
      if (is.null(uu)) uu <- rep(0, length(vv))
      
      plot_data[[vname]] <- list(
        val = vv,
        unc = uu,
        col = act_info[[vname]]$col,
        lab = act_info[[vname]]$lab
      )
      max_act <- max(max_act, max(vv + uu, na.rm = TRUE), na.rm = TRUE)
    } else {
      warning("   - Requested activity ", vname, " not found or all NA - skipped.")
    }
  }
  
  # ------------------------------------------------------------------
  # Proxy handling
  # ------------------------------------------------------------------
  has_proxy <- proxy_requested && 
    !is.null(List$Val$ProxyLabel) && 
    !all(is.na(List$Val$ProxyValue))
  
  proxy_col <- "#008837"   # defined unconditionally; only used when has_proxy = TRUE
  
  if (has_proxy) {
    X        <- List$Val$ProxyValue
    XU       <- List$Unc$ProxyValue
    if (is.null(XU)) XU <- rep(0, length(X))
    # ProxyLabel non-NULL is already guaranteed by the has_proxy condition above
    rawLabel <- List$Val$ProxyLabel
    XLabel   <- IsotopeLabel(rawLabel)
  }
  
  # ------------------------------------------------------------------
  # Plot boundaries setup
  # ------------------------------------------------------------------
  if (is.na(DepthMax)) DepthMax <- ceiling(max(Depth, na.rm = TRUE))
  if (is.na(DepthMin)) DepthMin <- 0
  
  if (length(plot_data) > 0) {
    if (is.na(PbMax) && is.na(RaMax)) {
      ActivityMax <- ceiling(max_act)
    } else {
      ActivityMax <- max(PbMax, RaMax, na.rm = TRUE)
    }
  } else {
    ActivityMax <- 1
  }
  
  if (has_proxy && is.na(ProxyMax)) {
    ProxyMax <- ceiling(max(X + XU, na.rm = TRUE))
  }
  
  draw_profile_plot <- function() {
    par(mar = c(4, 4, 4, 1))
    
    # Using type = "n" prevents artifacts from col = "white"
    plot(0, 0, ylim = c(DepthMax, DepthMin), xlim = c(0, ActivityMax),
         type = "n", xaxs = "i", yaxs = "i",
         axes = FALSE, ann = FALSE)
    
    axis(2, las = 1, cex = 1)
    axis(3)
    box()
    
    mtext("Depth (cm)", side = 2, line = 2.5)
    mtext(expression(paste("Activity (Bq kg"^{-1}, ")")), side = 3, line = 2.5)
    
    for (pd in plot_data) {
      arrows(x0 = pd$val - pd$unc, x1 = pd$val + pd$unc,
             y0 = Depth, y1 = Depth,
             angle = 90, code = 3, length = 0.03, col = pd$col)
      points(pd$val, Depth, col = "black", bg = pd$col, pch = 21)
    }
    
    if (has_proxy) {
      rad_labs <- lapply(plot_data, `[[`, "lab")
      all_labs <- c(rad_labs, list(XLabel))
      all_cols <- c(unlist(lapply(plot_data, `[[`, "col")), proxy_col)
      
      legend("bottomright", bty = "n",
             col = "black", pt.bg = all_cols, pch = 21, cex = 0.7,
             legend = do.call(c, all_labs))
      
      par(new = TRUE)
      
      # Using type = "n" prevents white ghost points from rendering over the grid
      plot(X, Depth, ylim = c(DepthMax, DepthMin), xlim = c(0, ProxyMax),
           type = "n", xaxs = "i", yaxs = "i",
           axes = FALSE, ann = FALSE)
      
      axis(1)
      mtext(XLabel, side = 1, line = 2.5)
      
      arrows(x0 = X - XU, x1 = X + XU, y0 = Depth, y1 = Depth,
             angle = 90, code = 3, length = 0.03, col = proxy_col)
      points(X, Depth, col = "black", bg = proxy_col, pch = 21)
      
    } else if (length(plot_data) > 0) {
      rad_labs <- lapply(plot_data, `[[`, "lab")
      rad_cols <- sapply(plot_data, `[[`, "col")
      legend("bottomright", bty = "n",
             col = "black", pt.bg = rad_cols, pch = 21, cex = 0.7,
             legend = do.call(c, rad_labs))
    }
  }
  
  # ------------------------------------------------------------------
  # Render Plots (PNG File & Screen)
  # ------------------------------------------------------------------
  if (SavePNG) {
    Name <- paste0(NamePlot, ".png")
    png(Name, res = 300, width = 10, height = 15, units = "cm")
    draw_profile_plot()
    dev.off()
    if (verbose) cat("   - PNG saved:", Name, "\n")
  }
  
  draw_profile_plot()
  
  invisible(NULL)
}
