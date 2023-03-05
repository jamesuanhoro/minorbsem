#' Stan data helper function
#' @description A function that creates data list object passed to Stan
#' @param lavaan_object lavaan fit object of corresponding model
#' @param method (character) One of "normal", "lasso",
#' "logistic", "GDP", "WB", "WB-cond", "none"
#' @param simple_struc (LOGICAL) Only relevant for CFAs.
#' If TRUE: assume simple structure;
#' If FALSE: estimate all cross-loadings using generalized
#' double Pareto priors.
#' @param priors An object of \code{\link{mbsempriors-class}}.
#' See \code{\link{new_mbsempriors}} for more information.
#' @returns Data list object used in fitting Stan model
#' @keywords internal
create_data_list <- function(
    lavaan_object = NULL,
    method = "normal",
    simple_struc = TRUE,
    priors = NULL) {
  data_list <- list()

  # Retrieve parameter structure from lavaan
  param_structure <- lavaan::lavInspect(lavaan_object)

  # Set method
  data_list$method <- method_hash(method)

  # Set simple structure to 0 by default, change within CFA section
  data_list$complex_struc <- 0

  # Has data?
  data_list$Y <- lavaan_object@Data@X[[1]]
  data_list$has_data <- ifelse(is.null(data_list$Y), 0, 1)

  methods::validObject(priors) # validate priors
  # Shape parameter for LKJ of interfactor corr
  data_list$shape_phi_c <- priors@lkj_shape
  data_list$sl_par <- priors@sl_par # sigma loading parameter
  data_list$rs_par <- priors@rs_par # residual sd parameter
  data_list$rc_par <- priors@rc_par # residual corr parameter
  data_list$sc_par <- priors@sc_par # sigma coefficients parameter
  data_list$fc_par <- priors@fc_par # factor correlation parameter

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
    data_list$complex_struc <- as.integer(ifelse(isFALSE(simple_struc), 1, 0))
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
  for (j in seq_len(ncol(data_list$loading_pattern))) {
    data_list$markers[j] <- which(
      data_list$loading_pattern[, j] == 1
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
