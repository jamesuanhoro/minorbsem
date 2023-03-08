#' Stan data helper function
#' @description A function that creates data list object passed to Stan
#' @param lavaan_object lavaan fit object of corresponding model
#' @inheritParams minorbsem
#' @returns Data list object used in fitting Stan model
#' @keywords internal
create_data_list_meta <- function(
    lavaan_object = NULL,
    method = "normal",
    simple_struc = TRUE,
    priors = NULL) {
  data_list <- list()

  # Get number of groups
  data_list$Ng <- lavaan::lavInspect(lavaan_object, "ngroups")

  if (data_list$Ng < 2) {
    stop(paste0(
      "Only one group found. ",
      "Meta-analysis requires multiple samples."
    ))
  }

  # Retrieve parameter structure from lavaan
  param_structure <- lavaan::lavInspect(lavaan_object)[[1]]

  # Set method
  data_list$method <- method_hash(method)

  # Set simple structure to 0 by default, change within CFA section
  data_list$complex_struc <- 0

  methods::validObject(priors) # validate priors
  # Shape parameter for LKJ of interfactor corr
  data_list$shape_phi_c <- priors@lkj_shape
  data_list$sl_par <- priors@sl_par # sigma loading parameter
  data_list$rs_par <- priors@rs_par # residual sd parameter
  data_list$rc_par <- priors@rc_par # residual corr parameter
  data_list$sc_par <- priors@sc_par # sigma coefficients parameter
  data_list$fc_par <- priors@fc_par # factor correlation parameter
  data_list$mln_par <- priors@mln_par # meta-reg intercept
  data_list$mlb_par <- priors@mlb_par # meta-reg coefficient

  # Sample cov
  data_list$S <- lapply(lavaan::lavInspect(
    lavaan_object, "SampStat"
  ), "[[", "cov")
  # Number of items
  data_list$Ni <- nrow(data_list$S[[1]])
  # Sample size
  data_list$Np <- lavaan::lavInspect(lavaan_object, "nobs")

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
    stop("Only CFAs are implemented right now.")
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

  # For now, no moderators
  data_list$p <- 0
  data_list$X <- matrix(nrow = data_list$Ng, ncol = data_list$p)

  # Handle missing data
  data_list <- prepare_missing_data_list(data_list)

  # 1 for meta-analysis
  data_list$meta <- 1

  return(data_list)
}
