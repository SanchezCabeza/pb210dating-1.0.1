#' @title Pb-210 dating of sediment cores
#'
#' @description
#' 
#' The \code{pb210dating} package provides a comprehensive toolset to date sediment
#' cores using the classical \eqn{^{210}}Pb models: the Constant Flux (CF)
#' and the Constant Flux Constant Sedimentation (CFCS) models.
#'
#' It includes functions for data import, decay correction, supported activity
#' (\eqn{^{226}}Ra) calculation, equilibrium estimation, Monte Carlo uncertainty propagation,
#' visualization, and export of results.
#'
#' @importFrom grDevices png dev.off colorRampPalette
#' @importFrom graphics par plot points arrows axis box mtext legend abline segments
#' @importFrom stats lm approx approxfun predict sd rnorm runif na.exclude coef setNames
#' @importFrom utils read.csv tail write.table
#' @importFrom lubridate decimal_date
#' 
#' @details
#' This package follows the nomenclature and formulations of
#' Sanchez-Cabeza and Ruiz-Fernandez (2012) and the Monte Carlo approach
#' described in Sanchez-Cabeza et al. (2014).
#'
#' \strong{Input}: 
#' The package uses a single CSV file with a strict 21-column structure
#' (see \code{\link{ReadData}} for details). The file contains:
#' \itemize{
#'   \item Columns 1–18: Core profile data (SampleCode, depths, masses, activities, dates, etc.).
#'   \item Columns 19–21: Metadata block (Core name, analysts, sampling date,
#'         diameter, half-lives, validation markers, etc.).
#' }
#' 
#' \strong{The dating object}:
#' All functions work by progressively updating a \strong{List} object containing
#' two sub-lists:
#' \itemize{
#'   \item \code{List$Val}: values (activities, ages, accumulation rates, etc.)
#'   \item \code{List$Unc}: corresponding uncertainties
#' }
#'
#' \strong{Recommended workflows}: 
#'
#' \strong{CF Model}:
#' \enumerate{
#'   \item \code{\link{ReadData}}
#'   \item \code{\link{DecayCorrection}} (if needed)
#'   \item \code{\link{CompleteProfile}} (optional)
#'   \item \code{\link{Equilibrium}} or \code{\link{EquilibriumSet}}
#'   \item \code{\link{Pb210CF}}
#'   \item \code{\link{AgeModel}} and \code{\link{PlotMARSAR}}
#'   \item \code{\link{Output}}
#' }
#'
#' \strong{CFCS Model}:
#' \enumerate{
#'   \item \code{\link{ReadData}}
#'   \item \code{\link{DecayCorrection}} (if needed)
#'   \item \code{\link{Pb210CFCS}} (with \code{Breaks})
#'   \item \code{\link{AgeModel}}
#'   \item \code{\link{Output}}
#' }
#'
#' The function \code{\link{Pepare4Dating}} performs the standard preprocessing 
#' workflow required before applying the \eqn{^{210}}Pb dating models.
#' It is a convenience wrapper that sequentially executes \code{\link{ReadData}}, 
#' \code{\link{DecayCorrection}}, and \code{\link{CompleteProfile}}. 
#' 
#' The function \code{\link{Pb210Dating}} is a convenience  wrapper to 
#' run the complete workflow for a sediment core suitable for 210Pb dating. 
#' This function could be used with a well-known core, 
#' in which the equilibrium depth is well constrained. 
#'
#' @author
#' \itemize{
#'   \item Joan-Albert Sanchez-Cabeza (\email{j.a.sanchez@@cmarl.unam.mx})
#'   \item Ana Carolina Ruiz-Fernández
#'   \item David Moriña Soler
#' }
#'
#' @references
#' Sanchez-Cabeza, J. A., & Ruiz-Fernandez, A. C. (2012).
#' \eqn{^{210}}Pb sediment radiochronology: an integrated formulation and
#' classification of dating models.
#' \emph{Geochimica et Cosmochimica Acta}, 82, 183-200.
#'
#' Sanchez-Cabeza, J. A., Ruiz-Fernández, A. C., Ontiveros-Cuadras, J. F., 
#' Bernal, L. H. P., & Olid, C. (2014). 
#' Monte Carlo uncertainty calculation of \eqn{^{210}}Pb chronologies and accumulation 
#' rates of sediments and peat bogs. \emph{Quaternary Geochronology}, 23, 80-93.
#'
#' @seealso
#' \code{\link{ReadData}},
#' \code{\link{Pb210CF}},
#' \code{\link{Pb210CFCS}},
#' \code{\link{AgeModel}}
#'
#' @examples
#' \dontrun{
#'   # Basic CF workflow
#'   List <- ReadData("my_core.csv")
#'   List <- DecayCorrection(List)
#'   List <- Equilibrium(List)
#'   List <- Pb210CF(List)
#'   AgeModel(List)
#'   Output(List)}
#'
#' @keywords internal
"_PACKAGE"
