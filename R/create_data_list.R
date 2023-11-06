#' Stan data helper function
#' @description A function that creates data list object passed to Stan
#' @param lavaan_object lavaan fit object of corresponding model
#' @inheritParams minorbsem
#' @returns Data list object used in fitting Stan model
#' @keywords internal
create_data_list <- function(
    lavaan_object = NULL,
    method = "normal",
    simple_struc = TRUE,
    priors = NULL,
    compute_ll = FALSE,
    model = NULL,
    orthogonal = NULL) {
  data_list <- list()

  # Retrieve parameter structure from lavaan
  param_structure <- lavaan::lavInspect(lavaan_object)

  # Set method
  data_list$method <- method_hash(method)

  # Set simple structure to 0 by default, change within CFA section
  data_list$complex_struc <- 0

  methods::validObject(priors) # validate priors
  # Shape parameter for LKJ of interfactor corr
  data_list$shape_phi_c <- priors@lkj_shape
  data_list$sl_par <- priors@sl_par # sigma loading parameter
  data_list$rm_par <- priors@rm_par # sigma(tau) parameter
  data_list$fs_par <- priors@fs_par # sigma factor sd parameter
  data_list$rs_par <- priors@rs_par # residual sd parameter
  data_list$rc_par <- priors@rc_par # residual corr parameter
  data_list$sc_par <- priors@sc_par # sigma coefficients parameter
  data_list$fc_par <- priors@fc_par # factor correlation parameter

  # Sample cov
  data_list$S <- lavaan_object@SampleStats@cov[[1]]
  rownames(data_list$S) <- rownames(param_structure$lambda)
  colnames(data_list$S) <- rownames(param_structure$lambda)
  # Number of items
  data_list$Ni <- nrow(data_list$S)
  # Sample size
  data_list$Np <- lavaan_object@SampleStats@nobs[[1]]

  # Has data?
  data_list$Y <- lavaan_object@Data@X[[1]]
  data_list$has_data <- 1
  if (is.null(data_list$Y)) {
    data_list$has_data <- 0
    data_list$Y <- matrix(nrow = 0, ncol = data_list$Ni)
  } else {
    data_list$Y <- t(t(data_list$Y) - colMeans(data_list$Y))
  }
  # return ll?
  # yes, when requested, data available and method != 90
  # TODO: add a warning
  data_list$ret_ll <- isTRUE(compute_ll) * data_list$has_data *
    (data_list$method != 90)

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

    param_structure <- lavaan::lavInspect(lavaan::cfa(
      model,
      sample.cov = data_list$S, sample.nobs = data_list$Np,
      std.lv = TRUE,
      do.fit = FALSE,
      likelihood = "wishart",
      do.fit = FALSE,
      ceq.simple = TRUE,
      orthogonal = orthogonal
    ))
    partab <- lavaan::lavaanify(
      model,
      ceq.simple = TRUE, std.lv = TRUE
    )

    # Loading pattern
    data_list$loading_pattern <- param_structure$lambda
    # Number of factors
    data_list$Nf <- ncol(data_list$loading_pattern)
    # location(loading) parameter
    data_list$load_est <- matrix(
      priors@ml_par, data_list$Ni, data_list$Nf
    )
    # sigma(loading) parameter
    data_list$load_se <- matrix(
      priors@sl_par, data_list$Ni, data_list$Nf
    )
    # get fixed loadings
    data_list$loading_fixed <- matrix(-999, data_list$Ni, data_list$Nf)
    fix_load <- partab[partab$op == "=~" & partab$free == 0, ]
    fix_col_ids <- unname(sapply(
      fix_load$lhs, function(x) which(x == colnames(param_structure$lambda))
    ))
    fix_row_ids <- unname(sapply(
      fix_load$rhs, function(x) which(x == rownames(param_structure$lambda))
    ))
    for (i in seq_len(length(fix_col_ids))) {
      data_list$loading_fixed[
        fix_row_ids[i], fix_col_ids[i]
      ] <- fix_load$ustart[i]
    }

    # res-var
    theta_var_diag <- diag(param_structure$theta)
    data_list$res_var_pattern <- theta_var_diag
    theta_zeroes <- theta_var_diag != 0
    if (sum(theta_zeroes) > 0) {
      theta_var_diag[theta_zeroes] <-
        theta_var_diag[theta_zeroes] - min(theta_var_diag[theta_zeroes]) + 1
      data_list$res_var_pattern <- theta_var_diag
    }
    data_list$res_var_fixed <- array(999, data_list$Ni)
    ind_names <- rownames(param_structure$theta)
    fix_rv <- partab[
      partab$op == "~~" & partab$free == 0 &
        partab$lhs %in% ind_names & partab$lhs == partab$rhs,
    ]
    fix_ind_ids <- unname(sapply(
      fix_rv$lhs, function(x) which(x == ind_names)
    ))
    for (i in seq_len(length(fix_ind_ids))) {
      data_list$res_var_fixed[fix_ind_ids[i]] <- fix_rv$ustart[i]
    }

    # Set to 0 for uncorrelated factors, 1 for correlated
    data_list$corr_mask <- diag(data_list$Nf)
    data_list$corr_mask[lower.tri(data_list$corr_mask)] <-
      (psi_mat[lower.tri(psi_mat)] != 0) + 0
    data_list$corr_mask[upper.tri(data_list$corr_mask)] <-
      t(data_list$corr_mask)[upper.tri(data_list$corr_mask)]
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
      data_list$loading_pattern[, j] != 0
    )[1]
  }

  # Check for correlated error terms
  # Number of correlated errors
  data_list$error_mat <- matrix(ncol = 2, nrow = 0)
  theta_corr_mat <- param_structure$theta
  diag(theta_corr_mat) <- 0
  theta_zeroes <- theta_corr_mat != 0
  if (sum(theta_zeroes) > 0) {
    theta_corr_mat[theta_zeroes] <-
      theta_corr_mat[theta_zeroes] - min(theta_corr_mat[theta_zeroes]) + 1
    data_list$error_mat <- which(theta_zeroes, arr.ind = TRUE)
    data_list$error_mat <- data_list$error_mat[
      data_list$error_mat[, 1] > data_list$error_mat[, 2], ,
      drop = FALSE
    ]
  }
  data_list$error_pattern <- theta_corr_mat
  data_list$Nce <- nrow(data_list$error_mat)

  return(data_list)
}
