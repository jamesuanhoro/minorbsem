methods::setClass("CmdStanMCMC")
methods::setClass("CmdStanFit")
methods::setClass("R6")
methods::setClassUnion(
  "cmdstan_classes", c("CmdStanMCMC", "CmdStanFit", "R6")
)

#' A class for setting up priors.
#'
#' @slot lkj_shape (positive real) The shape parameter of the LKJ-prior on the
#' interfactor correlation matrix in confirmatory factor models.
#' @slot ml_par (real) The location parameter of the normal prior on loadings.
#' @slot sl_par (positive real)
#' The scale parameter of the normal prior on loadings.
#' @slot rs_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the residual standard deviations.
#' @slot rc_par (positive real) The shape parameter of the Beta(rc_par, rc_par)
#' prior on the residual error correlations.
#' @slot sc_par (positive real)
#' The scale parameter of the normal prior on coefficients.
#' @slot rm_par (positive real) The scale parameter of the normal prior
#' on the tau / CRMR parameter.
#'
#' @name mbsempriors-class
#' @rdname mbsempriors-class
#' @export
methods::setClass(
  "mbsempriors",
  methods::representation(
    lkj_shape = "numeric",
    ml_par = "numeric",
    sl_par = "numeric",
    rs_par = "numeric",
    rc_par = "numeric",
    sc_par = "numeric",
    rm_par = "numeric"
  ),
  prototype = list(
    lkj_shape = 2.0,
    ml_par = 0.0,
    sl_par = 1.0,
    rs_par = 1.0,
    rc_par = 2.0,
    sc_par = 0.5,
    rm_par = 0.15
  )
)

#' Set priors in package
#'
#' @description Modify default priors in package.
#' @param lkj_shape (positive real) The shape parameter of the LKJ-prior on the
#' interfactor correlation matrix in confirmatory factor models.
#' @param ml_par (real) The location parameter of the normal prior on loadings.
#' @param sl_par (positive real)
#' The scale parameter of the normal prior on loadings.
#' @param rs_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the residual standard deviations.
#' @param rc_par (positive real) The shape parameter of the Beta(rc_par, rc_par)
#' prior on the residual error correlations.
#' @param sc_par (positive real)
#' The scale parameter of the normal prior on coefficients.
#' @param rm_par (positive real) The scale parameter of the normal prior
#' on the tau / CRMR parameter.
#' @returns An object of \code{\link{mbsempriors-class}}
#' @examples
#' # Change LKJ shape parameter only
#' custom_priors <- new_mbsempriors(lkj_shape = 1.0)
#' \dontrun{
#' minorbsem("F1 =~ x1 + x2 + x3
#'            F2 =~ x4 + x5 + x6
#'            F3 =~ x7 + x8 + x9", HS,
#'   priors = custom_priors
#' )
#' }
#' @export
new_mbsempriors <- function(
    lkj_shape = 2.0,
    ml_par = 0.0,
    sl_par = 1.0,
    rs_par = 1.0,
    rc_par = 2.0,
    sc_par = 0.5,
    rm_par = 0.15) {
  mb_priors_object <- methods::new("mbsempriors")
  mb_priors_object <- methods::initialize(
    mb_priors_object,
    lkj_shape = lkj_shape,
    ml_par = ml_par,
    sl_par = sl_par,
    rs_par = rs_par,
    rc_par = rc_par,
    sc_par = sc_par,
    rm_par = rm_par
  )
  return(mb_priors_object)
}

methods::setClass(
  "mbsem",
  methods::representation(
    major_parameters = "data.frame",
    minor_factor_matrix = "data.frame",
    data_list = "list",
    priors = "mbsempriors",
    stan_fit = "cmdstan_classes",
    version = "character"
  )
)

new_mbsem <- function() {
  methods::new("mbsem")
}
