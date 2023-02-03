create_data_list <- function(
    lavaan_object = NULL,
    lkj_shape = 2,
    sl_par = 1,
    rs_par = 2.5,
    rc_par = 2.0) {
  data_list <- list()

  # Retrieve parameter structure from lavaan
  param_structure <- lavaan::lavInspect(lavaan_object)

  # Has data?
  data_list$Y <- lavaan_object@Data@X[[1]]
  data_list$has_data <- ifelse(is.null(data_list$Y), 0, 1)

  # Set up priors
  if (any(c(lkj_shape, sl_par, rs_par, rc_par) <= 0)) {
    stop("lkj_shape, sl_par, rs_par, rc_par must all exceed 0")
  }
  data_list$shape_phi_c <- lkj_shape # Shape parameter for LKJ of interfactor corr
  data_list$sl_par <- sl_par # sigma loading parameter
  data_list$rs_par <- rs_par # residual sd parameter
  data_list$rc_par <- rc_par # residual corr parameter

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
  Psi <- param_structure$psi
  sum_off_diag_psi <- sum(Psi[lower.tri(Psi)])
  if (is.null(param_structure$beta)) {
    # This is a CFA
    # Set to 0 for uncorrelated factors, 1 for correlated
    data_list$corr_fac <- ifelse(sum_off_diag_psi == 0, 0, 1)
  } else {
    # This is an SEM (on hold right now)
    # data_list$Nf_corr <- 0
    # (data_list$F_corr_mat <- matrix(byrow = TRUE, ncol = 2, nrow = data_list$Nf_corr))
  }

  # Marker variables per factor
  # Each factor should have one unique indicator or stop!
  data_list$markers <- array(dim = data_list$Nf)
  unique_indicators <- which(rowSums(data_list$loading_pattern) == 1)
  unique_indicators <- data_list$loading_pattern[unique_indicators, ]
  if (any(colSums(unique_indicators) == 0)) {
    notice <- paste0(
      "Each factor must have at least one indicator unique to it.", "\n",
      "This is to ensure the sign/direction of the factor does not flip across iterations.", "\n",
      "Also, note that the package can only fit standard CFAs and SEM models,", "\n",
      "no higher-order factors, MIMIC, multilevel SEM, path analysis models, ..."
    )
    stop(notice)
  }
  unique_indicators <- rowSums(data_list$loading_pattern) == 1
  for (j in 1:ncol(data_list$loading_pattern)) {
    data_list$markers[j] <- which(
      data_list$loading_pattern[, j] == 1 & unique_indicators
    )[1]
  }

  # Check for correlated error terms
  # Number of correlated errors
  data_list$error_mat <- matrix(byrow = TRUE, ncol = 2, nrow = 0)
  Theta <- param_structure$theta
  sum_off_diag_theta <- sum(Theta[lower.tri(Theta)])
  if (sum_off_diag_theta > 0) {
    # Get which elements in theta are non-zero
    data_list$error_mat <- which(Theta != 0, arr.ind = TRUE)
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
