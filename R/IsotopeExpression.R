# =============================================================================
# IsotopeExpression.R
# Radionuclide label expressions for plotting in the pb210dating package.
# =============================================================================

# -----------------------------------------------------------------------------
# .ParseIsotopeString  (internal, not exported)
# -----------------------------------------------------------------------------
# Parse a plain-text isotope string into its mass and element components.
#
# Accepted formats:
#   "137Cs"      mass number immediately followed by element symbol
#   "Cs137"      element symbol followed by mass number
#   "Cs-137"     element symbol, dash, mass number
#   "cs137"      case-insensitive variants of any of the above
#   "239,240Pu"  comma-separated mass numbers (combined isotopes, one element)
#   "239+240Pu"  plus-separated mass numbers
#   "239-240Pu"  hyphen-separated mass numbers
#
# Returns a list with:
#   $mass    character — mass number string, normalised to comma separator
#   $element character — element symbol, first letter capitalised
# -----------------------------------------------------------------------------
.ParseIsotopeString <- function(x) {
  
  x <- trimws(x)
  
  known_elements <- c(
    "H","He","Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
    "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
    "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
    "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
    "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
    "Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
  )
  
  # Pattern 1: mass-first  e.g. "137Cs", "239,240Pu", "239+240Pu"
  m1 <- regmatches(x, regexec("^([0-9][0-9,+\\-]*)([A-Za-z]{1,3})$", x))[[1]]
  
  # Pattern 2: element-first  e.g. "Cs137", "Cs-137", "Pb-210"
  m2 <- regmatches(x, regexec("^([A-Za-z]{1,3})[\\-]?([0-9][0-9,+\\-]*)$", x))[[1]]
  
  if (length(m1) == 3) {
    mass_raw <- m1[2]; elem_raw <- m1[3]
  } else if (length(m2) == 3) {
    elem_raw <- m2[2]; mass_raw <- m2[3]
  } else {
    stop(sprintf(
      "Cannot parse isotope string '%s'.\n  Accepted formats: '137Cs', 'Cs137', 'Cs-137', '239,240Pu', '239+240Pu'.",
      x))
  }
  
  elem <- paste0(toupper(substr(elem_raw, 1, 1)),
                 tolower(substr(elem_raw, 2, nchar(elem_raw))))
  
  if (!elem %in% known_elements)
    warning(sprintf("'%s' is not a recognised element symbol.", elem))
  
  # Normalise separators (+, -) to comma for superscript rendering
  mass_clean <- gsub("[+\\-]", ",", mass_raw)
  mass_clean <- gsub(",+$", "", mass_clean)
  
  list(mass = mass_clean, element = elem)
}

# -----------------------------------------------------------------------------
# IsotopeName  (internal, not exported)
# -----------------------------------------------------------------------------
# Such as 210Pbex (Bq kg-1)
# -----------------------------------------------------------------------------
IsotopeName <- function(x, suffix = NULL, units = NULL) {
  p <- .ParseIsotopeString(x)
  lbl <- paste0(p$mass, p$element)
  if (!is.null(suffix) && nzchar(suffix))
    lbl <- paste0(lbl, suffix)
  if (!is.null(units) && nzchar(units))
    lbl <- paste0(lbl, " (", units, ")")
  lbl
}

# -----------------------------------------------------------------------------
# IsotopeUncertainty  (internal, not exported)
# -----------------------------------------------------------------------------
# Such as u(210Pbex)
# -----------------------------------------------------------------------------

IsotopeUncertainty <- function(x, suffix = NULL) {
  p <- .ParseIsotopeString(x)
  lbl <- paste0(p$mass, p$element)
  if (!is.null(suffix) && nzchar(suffix))
    lbl <- paste0(lbl, suffix)
  paste0("u(", lbl, ")")
}

# -----------------------------------------------------------------------------
# .BuildExprString  (internal, not exported)
# -----------------------------------------------------------------------------
# Assemble a single plotmath string from parsed components.
# -----------------------------------------------------------------------------

.BuildExprString <- function(mass, elem, suffix, units) {
  
  s <- sprintf('{}^"%s"*"%s"', mass, elem)
  
  if (!is.null(suffix) && nzchar(suffix)) {
    if (startsWith(suffix, "_")) {
      sc <- sub("^_\\{?", "", suffix); sc <- sub("\\}$", "", sc)
      s  <- paste0(s, sprintf('["%s"]', sc))
    } else if (startsWith(suffix, "^")) {
      sc <- sub("^\\^\\{?", "", suffix); sc <- sub("\\}$", "", sc)
      s  <- paste0(s, sprintf('^"%s"', sc))
    } else {
      s <- paste0(s, sprintf('*"(%s)"', suffix))
    }
  }
  
  if (!is.null(units) && nzchar(units))
    s <- paste0(s, "~", units)
  
  s
}


#' Isotope label expression for plot axes and titles
#'
#' Generates a \code{\link[grDevices]{plotmath}} expression that renders a
#' radionuclide symbol with the mass number as a superscript immediately
#' preceding the element symbol (e.g. \eqn{^{137}}Cs, \eqn{^{210}}Pb), in
#' accordance with IUPAC typographic conventions.  The returned object is
#' accepted directly by the \code{xlab}, \code{ylab}, and \code{main}
#' arguments of all base-R plotting functions, as well as by
#' \code{\link[graphics]{axis}}, \code{\link[graphics]{mtext}}, and
#' \code{\link[graphics]{text}}.
#'
#' @param x A single character string specifying the isotope in any of the
#'   plain-text notations commonly used in the literature.  Accepted formats
#'   are illustrated below; matching is case-insensitive for the element
#'   symbol. Alternatively, a two-element character vector
#'   \code{c(mass, element)} may be supplied, e.g. \code{c("239,240", "Pu")}.
#'
#'   \describe{
#'     \item{\code{"137Cs"}}{Mass number immediately followed by element symbol.}
#'     \item{\code{"Cs137"}, \code{"Cs-137"}}{Element symbol before mass number,
#'       optionally separated by a hyphen.}
#'     \item{\code{"239,240Pu"}}{Comma-separated mass numbers for combined
#'       isotopes sharing one element symbol (e.g. \eqn{^{239,240}}Pu).}
#'     \item{\code{"239+240Pu"}, \code{"239-240Pu"}}{Plus- or hyphen-separated
#'       mass numbers; both are normalised to a comma in the rendered label.}
#'   }
#'
#' @param suffix An optional character string appended immediately after the
#'   element symbol.  Three formats are recognised:
#'   \describe{
#'     \item{\code{"_{ex}"}}{Renders as a plotmath subscript:
#'       \eqn{^{210}}Pb\eqn{_{ex}}.}
#'     \item{\code{"^{*}"}}{Renders as a plotmath superscript after the symbol: 
#'       \eqn{^{210}}Pb\eqn{"^{*}"}.}
#'     \item{any other string}{Rendered in parentheses:
#'       \eqn{^{210}}Pb(total).}
#'   }
#'   Default \code{NULL} (no suffix).
#'
#' @param units An optional plotmath string appended after the isotope label,
#'   separated by a thin space (\code{~}).  Standard plotmath syntax applies:
#'   use \code{~} for spaces and \code{^{-1}} for negative exponents, e.g.
#'   \code{"Bq~kg^{-1}"}. Default \code{NULL} (no units).
#'
#' @return A parsed \code{\link[base]{expression}} object of length 1, suitable
#'   for direct use wherever R's plotmath rendering is needed.
#'
#' @seealso \code{\link{IsotopeLegend}} for legend entries;
#'   \code{\link[grDevices]{plotmath}} for the full plotmath syntax reference.
#'
#' @examples
#' # Axis label with subscript suffix and units
#' plot(1:10, runif(10),
#'      ylab = IsotopeLabel("Pb-210", suffix = "_{ex}", units = "Bq~kg^{-1}"))
#'
#' @export
IsotopeLabel <- function(x, suffix = NULL, units = NULL) {
  
  if (length(x) == 2 && !grepl("[A-Za-z]", x[1])) {
    parsed <- list(mass = x[1], element = x[2])
  } else {
    parsed <- .ParseIsotopeString(x)
  }
  
  parse(text = .BuildExprString(parsed$mass, parsed$element, suffix, units))
}


#' Isotope label expressions for plot legends
#'
#' Generates a single \code{\link[base]{expression}} object of length
#' \eqn{n} containing one typeset radionuclide label per entry, suitable for
#' direct use as the \code{legend} argument of
#' \code{\link[graphics]{legend}}.
#'
#' @details
#' \code{\link[graphics]{legend}} requires its \code{legend} argument to be
#' either a plain character vector or a \emph{single} \code{expression} object
#' whose elements correspond one-to-one to legend entries.  This function
#' ensures that requirement by calling \code{\link[base]{parse}} on the full
#' character vector of plotmath strings at once, which returns an expression of
#' the correct length.
#'
#' Each element is formatted by the same rules as \code{\link{IsotopeLabel}};
#' the \code{suffix} and \code{units} arguments are recycled across all entries
#' if a vector shorter than \code{isotopes} is supplied.
#'
#' @param isotopes A character vector of isotope strings in any format accepted
#'   by \code{\link{IsotopeLabel}}.
#' @param suffix A character vector of suffix strings (recycled).  See the
#'   \code{suffix} argument of \code{\link{IsotopeLabel}} for recognised
#'   formats.  Default \code{NULL} (no suffix for any entry).
#' @param units A character vector of plotmath unit strings (recycled).
#'   Default \code{NULL} (no units for any entry).
#'
#' @return A parsed \code{\link[base]{expression}} object of length
#'   \code{length(isotopes)}, ready to be passed as \code{legend = } in a call
#'   to \code{\link[graphics]{legend}}.
#'
#' @seealso \code{\link{IsotopeLabel}} for single labels on axes and titles.
#'
#' @examples
#' plot(1:10, runif(10), type = "n",
#'      xlab = "depth (cm)", ylab = "activity (Bq/kg)")
#' points(1:10, runif(10), col = "steelblue", pch = 16)
#' points(1:10, runif(10), col = "firebrick", pch = 17)
#' points(1:10, runif(10), col = "darkgreen", pch = 15)
#' legend("topright",
#'        legend = IsotopeLegend(c("137Cs", "241Am", "239,240Pu")),
#'        col    = c("steelblue", "firebrick", "darkgreen"),
#'        pch    = c(16, 17, 15))
#'
#' # With per-entry suffixes and shared units
#' legend("bottomleft",
#'        legend = IsotopeLegend(c("210Pb", "210Pb", "Ra-226"),
#'                                suffix = c("_{ex}", "_{total}", ""),
#'                                units  = "Bq~kg^{-1}"),
#'        col = c("black", "grey50", "red"),
#'        lty = 1)
#'
#' @export
IsotopeLegend <- function(isotopes, suffix = NULL, units = NULL) {
  
  n <- length(isotopes)
  if (is.null(suffix)) suffix <- rep("", n)
  if (is.null(units))  units  <- rep("", n)
  suffix <- rep_len(suffix, n)
  units  <- rep_len(units,  n)
  
  expr_strings <- vapply(seq_len(n), function(i) {
    parsed <- .ParseIsotopeString(isotopes[i])
    .BuildExprString(parsed$mass, parsed$element,
                       if (nzchar(suffix[i])) suffix[i] else NULL,
                       if (nzchar(units[i]))  units[i]  else NULL)
  }, character(1))
  
  # parse() on a length-n character vector returns expression of length n,
  # which is the type required by legend().
  parse(text = expr_strings)
}
