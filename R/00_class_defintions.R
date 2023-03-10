methods::setClass("CmdStanMCMC")
methods::setClass("CmdStanFit")
methods::setClass("R6")
methods::setClassUnion("stan_fit_classes", c("CmdStanMCMC", "CmdStanFit", "R6"))

#' A class for setting up priors.
#'
#' @slot lkj_shape (positive real) The shape parameter of the LKJ-prior on the
#' interfactor correlation matrix in confirmatory factor models.
#' @slot sl_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the hyper-parameter of the standard
#' deviation of loadings.
#' @slot rs_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the residual standard deviations.
#' @slot rc_par (positive real) The shape parameter of the Beta(rc_par, rc_par)
#' prior on the residual error correlations.
#' @slot sc_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the hyper-parameter of the standard
#' deviations of coefficients; SD(coefs) vary by outcome.
#' @slot fc_par (positive real) The shape parameter of the Beta(rc_par, rc_par)
#' prior on the inter-factor correlations in latent regression models.
#' @slot mln_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the
#' log-scale meta-regression intercept.
#' Relevant only for meta-analysis.
#' @slot mlb_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the
#' log-scale meta-regression coefficients.
#' Relevant only for meta-analysis.
#'
#' @name mbsempriors-class
#' @rdname mbsempriors-class
#' @export
methods::setClass(
  "mbsempriors",
  methods::representation(
    lkj_shape = "numeric",
    sl_par = "numeric",
    rs_par = "numeric",
    rc_par = "numeric",
    sc_par = "numeric",
    fc_par = "numeric",
    mln_par = "numeric",
    mlb_par = "numeric"
  ),
  prototype = list(
    lkj_shape = 2.0,
    sl_par = 1.0,
    rs_par = 2.5,
    rc_par = 2.0,
    sc_par = 1.0,
    fc_par = 2.0,
    mln_par = 2.5,
    mlb_par = 1.0
  )
)

#' Set priors in package
#'
#' @description Modify default priors in package.
#' @param lkj_shape (positive real) The shape parameter of the LKJ-prior on the
#' interfactor correlation matrix in confirmatory factor models.
#' @param sl_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the hyper-parameter of the standard
#' deviation of loadings.
#' @param rs_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the residual standard deviations.
#' @param rc_par (positive real) The shape parameter of the Beta(rc_par, rc_par)
#' prior on the residual error correlations.
#' @param sc_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the hyper-parameter of the standard
#' deviations of coefficients; SD(coefs) vary by outcome.
#' @param fc_par (positive real) The shape parameter of the Beta(rc_par, rc_par)
#' prior on the inter-factor correlations in latent regression models.
#' @param mln_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the
#' log-scale meta-regression intercept.
#' Relevant only for meta-analysis.
#' @param mlb_par (positive real) The scale parameter of the
#' Student-t(df = 3, loc = 0) prior on the
#' log-scale meta-regression coefficients.
#' Relevant only for meta-analysis.
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
    sl_par = 1.0,
    rs_par = 2.5,
    rc_par = 2.0,
    sc_par = 1.0,
    fc_par = 2.0,
    mln_par = 2.5,
    mlb_par = 1.0) {
  mb_priors_object <- methods::new("mbsempriors")
  mb_priors_object <- methods::initialize(
    mb_priors_object,
    lkj_shape = lkj_shape,
    sl_par = sl_par,
    rs_par = rs_par,
    rc_par = rc_par,
    sc_par = sc_par,
    fc_par = fc_par,
    mln_par = mln_par,
    mlb_par = mlb_par
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
    stan_fit = "stan_fit_classes",
    version = "character"
  )
)

new_mbsem <- function() {
  methods::new("mbsem")
}
