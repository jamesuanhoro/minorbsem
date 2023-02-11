#' Stan data helper function
#' @description A function that creates data list object passed to Stan
#' @param lavaan_object lavaan fit object of corresponding model
#' @param method (character) One of "normal", "lasso"
#' @param lkj_shape (positive real) The shape parameter of the LKJ-prior on the
#' interfactor correlation matrix.
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
#' @returns Data list object used in fitting Stan model
#' @keywords internal
create_data_list <- function(
    lavaan_object = NULL,
    method = "normal",
    lkj_shape = 2.0,
    sl_par = 1.0,
    rs_par = 2.5,
    rc_par = 2.0,
    sc_par = 1.0) {
  data_list <- list()

  # Retrieve parameter structure from lavaan
  param_structure <- lavaan::lavInspect(lavaan_object)

  # Set method
  data_list$method <- method_hash(method)

  # Has data?
  data_list$Y <- lavaan_object@Data@X[[1]]
  data_list$has_data <- ifelse(is.null(data_list$Y), 0, 1)

  # Set up priors
  if (any(c(lkj_shape, sl_par, rs_par, rc_par, sc_par) <= 0)) {
    stop("lkj_shape, sl_par, rs_par, rc_par, sc_par must all exceed 0")
  }
  # Shape parameter for LKJ of interfactor corr
  data_list$shape_phi_c <- lkj_shape
  data_list$sl_par <- sl_par # sigma loading parameter
  data_list$rs_par <- rs_par # residual sd parameter
  data_list$rc_par <- rc_par # residual corr parameter
  data_list$sc_par <- sc_par # sigma coefficients parameter
  # TODO: either make this an argument or set std.lv = TRUE
  data_list$shape_beta <- 2.0

  # Sample cov
  data_list$S <- lavaan_object@SampleStats@cov[[1]]
  # Number of items
  data_list$Ni <- nrow(data_list$S)
  # Sample size
  data_list$Np <- lavaan_object@SampleStats@nobs[[1]]

  # Loading pattern, 0s and 1s
  data_list$loading_pattern <- (param_structure$lambda > 0) * 1
  # Number of factors
  data_list$Nf <- ncol(data_list$loading_pattern)

  # Is this an SEM or a CFA?
  psi_mat <- param_structure$psi
  sum_off_diag_psi <- sum(psi_mat[lower.tri(psi_mat)])
  if (is.null(param_structure$beta)) {
    # This is a CFA
    data_list$sem_indicator <- 0
    # Set to 0 for uncorrelated factors, 1 for correlated
    data_list$corr_fac <- ifelse(sum_off_diag_psi == 0, 0, 1)
  } else {
    # This is an SEM
    data_list$sem_indicator <- 1
    # Factor correlation matrix
    data_list$F_corr_mat <- matrix(byrow = TRUE, ncol = 2, nrow = 0)
    if (sum_off_diag_psi > 0) {
      # Get which elements in theta are non-zero
      data_list$F_corr_mat <- which(psi_mat != 0, arr.ind = TRUE)
      # Eliminate diagonal elements
      data_list$F_corr_mat <- data_list$F_corr_mat[
        data_list$F_corr_mat[, 1] != data_list$F_corr_mat[, 2],
      ]
      # Eliminate duplicate rows
      data_list$F_corr_mat <- unique(t(apply(data_list$F_corr_mat, 1, sort)))
    }
    data_list$Nf_corr <- nrow(data_list$F_corr_mat)
    data_list$coef_pattern <- (param_structure$beta > 0) * 1
  }

  # Marker variables per factor
  # Each factor should have one unique indicator or stop!
  data_list$markers <- array(dim = data_list$Nf)
  unique_indicators <- which(rowSums(data_list$loading_pattern) == 1)
  unique_indicators <- as.matrix(data_list$loading_pattern[unique_indicators, ])
  if (any(colSums(unique_indicators) == 0)) {
    notice <- paste0(
      "Each factor must have at least one indicator unique to it.", "\n",
      "This is to ensure the sign/direction of the factor ",
      "does not flip across iterations.", "\n",
      "Also, note that the package can only fit standard ",
      "CFAs and SEM models,", "\n",
      "no higher-order factors, MIMIC, multilevel SEM, ",
      "path analysis models, ..."
    )
    stop(notice)
  }
  unique_indicators <- rowSums(data_list$loading_pattern) == 1
  for (j in seq_len(ncol(data_list$loading_pattern))) {
    data_list$markers[j] <- which(
      data_list$loading_pattern[, j] == 1 & unique_indicators
    )[1]
  }

  # Check for correlated error terms
  # Number of correlated errors
  data_list$error_mat <- matrix(byrow = TRUE, ncol = 2, nrow = 0)
  theta_mat <- param_structure$theta
  sum_off_diag_theta <- sum(theta_mat[lower.tri(theta_mat)])
  if (sum_off_diag_theta > 0) {
    # Get which elements in theta are non-zero
    data_list$error_mat <- which(theta_mat != 0, arr.ind = TRUE)
    # Eliminate diagonal elements
    data_list$error_mat <- data_list$error_mat[
      data_list$error_mat[, 1] != data_list$error_mat[, 2],
    ]
    # Eliminate duplicate rows
    data_list$error_mat <- unique(t(apply(data_list$error_mat, 1, sort)))
  }
  data_list$Nce <- nrow(data_list$error_mat)

  return(data_list)
}
