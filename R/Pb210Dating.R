#' @title Complete Pb-210 dating workflow
#'
#' @description
#' Function wrapper to run the complete workflow for a sediment core suitable for \eqn{^{210}}Pb dating. 
#' This function is suitable for calculations with a well-known core, 
#' in which the equilibrium depth is well constrained. 
#' If run without careful core exploration, 
#' the risk of bad results is high and it is discouraged.
#'
#' @param FileName Character. CSV file containing the sediment-core data.
#' Default = NA uses the package example.
#' @param LayerDating Numeric. Bottom depth (cm) of the dating interval.
#' If NA (default), the equilibrium layer is selected automatically
#' using Equilibrium().
#' @param verbose Logical.Print progress messages. Default is \code{TRUE}.
#'
#' @returns
#' Dating list containing all calculated variables from the CF and CFCS
#' models.
#' 
#' @examples 
#' \dontrun{
#' Pb210Dating(FileName = List, LayerDating = 17)}
#'
#' @export

Pb210Dating <- function(FileName = NA, 
                        LayerDating = NA,
                        verbose = TRUE)
{
  if (verbose)
    cat("Starting complete Pb-210 dating workflow.\n")
  
  ## Read data
  List <- ReadData(NameFile = FileName, verbose = verbose)
  
  ## Plot original profiles
  PlotProfiles(List, PlotActivities = c("Pb210Total", "Po210Total", "Ra226", "Proxy"),
    NamePlot = paste0(List$Val$Core, "_profiles"), verbose = verbose)
  
  ## Estimate supported Ra if needed
  if (is.null(List$Val$Ra226) || all(is.na(List$Val$Ra226))) 
    List <- ConstantRa(List, mode = "auto", verbose = verbose)
  
  ## Decay correction
  List <- DecayCorrection(List, verbose = verbose)
  
  ## Complete missing profile
  List <- CompleteProfile(List, verbose = verbose)
  
  ## Select dating interval
  if (is.na(LayerDating)) {
    List <- Equilibrium(List, Interactive = FALSE, verbose = verbose)
  } else {
    List <- EquilibriumSet(List, LayerDating = LayerDating, verbose = verbose)
  }
  
  ## CFCS model
  List <- Pb210CFCS(List, verbose = verbose)
  
  ## CF model
  List <- Pb210CF(List, verbose = verbose)
  
  ## Age model
  AgeModel(List, PlotName = paste0(List$Val$Core, "_AgeModel"), verbose = verbose)
  
  ## Output
  Output(List)
  
  if (verbose)
    cat("Pb-210 dating workflow completed.\n")
  
  invisible(List)
}
