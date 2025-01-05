#' Stan data helper function for path analysis
#' @description A function that creates data list object passed to Stan
#' @param lavaan_object lavaan fit object of corresponding model
#' @param partab lavaanify result of corresponding model
#' @param old_names (Optional) Variable name order of original
#' correlation matrix, used to reorder acov_mat.
#' @inheritParams minorbsem
#' @returns Data list object used in fitting Stan model
#' @keywords internal
create_data_list_pa <- function(
    lavaan_object = NULL,
    method = "normal",
    correlation = correlation,
    priors = NULL,
    compute_ll = FALSE,
    partab = NULL,
    centered = TRUE,
    acov_mat = NULL,
    old_names = NULL) {
  data_list <- list()

  # Retrieve parameter structure from lavaan
  param_structure <- lavaan::lavInspect(lavaan_object)

  # Set method
  data_list$method <- method_hash(method)

  methods::validObject(priors) # validate priors
  # Shape parameter for LKJ of interfactor corr
  data_list$shape_phi_c <- priors@lkj_shape
  data_list$sl_par <- priors@sl_par # sigma loading parameter
  data_list$rm_par <- priors@rm_par # sigma(tau) parameter
  data_list$rs_par <- priors@rs_par # residual sd parameter
  data_list$rc_par <- priors@rc_par # residual corr parameter
  data_list$sc_par <- priors@sc_par # sigma coefficients parameter

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
  # yes, when requested, data available and method != 90 and NOT Corr-Struc-Ana
  # TODO: add a warning
  data_list$ret_ll <- isTRUE(compute_ll) * data_list$has_data *
    (data_list$method != 90) * !isTRUE(correlation)

  # Loading pattern, 0s and 1s
  data_list$loading_pattern <- (param_structure$lambda > 0) * 1

  # Assume CFA by default
  data_list$pa_indicator <- 1
  data_list$sem_indicator <- 0

  # Loading pattern
  data_list$loading_pattern <- param_structure$lambda
  # Number of factors
  data_list$Nf <- ncol(data_list$loading_pattern)
  # location(loading) parameter
  data_list$load_est <- matrix(
    priors@ml_par, data_list$Ni, data_list$Nf,
    dimnames = dimnames(data_list$loading_pattern)
  )
  # sigma(loading) parameter
  data_list$load_se <- matrix(
    priors@sl_par, data_list$Ni, data_list$Nf,
    dimnames = dimnames(data_list$loading_pattern)
  )
  # get fixed loadings
  data_list$loading_fixed <- matrix(
    -999, data_list$Ni, data_list$Nf,
    dimnames = dimnames(data_list$loading_pattern)
  )
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

  all_zero_loadings <- all(data_list$loading_pattern == 0)
  all_zero_loadings_fix <- all(data_list$loading_fixed == -999)
  if (isFALSE(all_zero_loadings) || isFALSE(all_zero_loadings_fix)) {
    err_msg <- paste0(
      "This model has latent variables, ",
      "this function is only for path analysis"
    )
    stop(err_msg)
  }

  # res-var
  theta_var_diag <- diag(param_structure$psi)
  data_list$res_var_pattern <- theta_var_diag
  theta_zeroes <- theta_var_diag != 0
  if (sum(theta_zeroes) > 0) {
    theta_var_diag[theta_zeroes] <-
      theta_var_diag[theta_zeroes] - min(theta_var_diag[theta_zeroes]) + 1
    data_list$res_var_pattern <- theta_var_diag
  }
  data_list$res_var_fixed <- array(999, data_list$Ni)
  ind_names <- rownames(param_structure$psi)
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

  # Set up coefficient table
  data_list$coef_pattern <- matrix(0, data_list$Ni, data_list$Ni)
  data_list$coef_est <- matrix(0, data_list$Ni, data_list$Ni)
  data_list$coef_se <- matrix(priors@sc_par, data_list$Ni, data_list$Ni)
  data_list$coef_fixed <- matrix(-999, data_list$Ni, data_list$Ni)
  if (!is.null(param_structure$beta)) {
    # This is an SEM
    data_list$sem_indicator <- 1
    # Factor coefficient matrix
    data_list$coef_pattern <- param_structure$beta
    dimnames(data_list$coef_est) <- dimnames(data_list$coef_pattern)
    dimnames(data_list$coef_se) <- dimnames(data_list$coef_pattern)
    dimnames(data_list$coef_fixed) <- dimnames(data_list$coef_pattern)
    beta_zeroes <- data_list$coef_pattern != 0
    data_list$coef_pattern[beta_zeroes] <-
      data_list$coef_pattern[beta_zeroes] -
      min(data_list$coef_pattern[beta_zeroes]) + 1
    # get fixed coefs
    fix_coef <- partab[partab$op == "~" & partab$free == 0, ]
    fix_col_ids <- unname(sapply(
      fix_coef$lhs, function(x) which(x == colnames(param_structure$lambda))
    ))
    fix_row_ids <- unname(sapply(
      fix_coef$rhs, function(x) which(x == colnames(param_structure$lambda))
    ))
    for (i in seq_len(length(fix_col_ids))) {
      data_list$coef_fixed[
        fix_row_ids[i], fix_col_ids[i]
      ] <- fix_coef$ustart[i]
    }
  } else {
    err_msg <- "Model does not specify any coefficients. Not a path analysis"
    stop(err_msg)
  }

  # Check for correlated error terms
  # Number of correlated errors
  data_list$error_mat <- matrix(ncol = 2, nrow = 0)
  theta_corr_mat <- param_structure$psi
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

  # Check for conditional independence
  cond_ind_list <- dagitty::impliedConditionalIndependencies(
    dagitty::lavaanToGraph(partab)
  )
  ind_names <- rownames(param_structure$psi)
  data_list$cond_ind_mat <- matrix(
    0,
    nrow = nrow(param_structure$psi), ncol = ncol(param_structure$psi),
    dimnames = dimnames(param_structure$psi)
  )
  if (length(cond_ind_list) > 0) {
    for (i in seq_along(cond_ind_list)) {
      ci <- cond_ind_list[[i]]
      x_ind <- which(ind_names == ci$X)
      y_ind <- which(ind_names == ci$Y)
      data_list$cond_ind_mat[x_ind, y_ind] <- 1
      data_list$cond_ind_mat[y_ind, x_ind] <- 1
    }
  }
  # don't for any correlations estimated by lavaan
  data_list$cond_ind_mat <-
    data_list$cond_ind_mat * (data_list$error_pattern == 0)

  data_list$correlation <- 0
  if (isTRUE(correlation)) {
    data_list$centered <- ifelse(isFALSE(centered), 0, 1)
    data_list$correlation <- 1
    r_mat <- stats::cov2cor(data_list$S)
    r_vec <- r_mat[lower.tri(data_list$S, diag = FALSE)]
    data_list$r_obs_vec <- .g_map(r_vec)
    tmp_acov <- acov_mat
    if (!is.null(tmp_acov)) {
      if (!isTRUE(all.equal(old_names, rownames(param_structure$lambda)))) {
        tmp_acov <- .fix_acov(
          acov_mat, old_names, rownames(param_structure$lambda)
        )
      }
    }
    data_list$r_obs_vec_cov <- .get_avar_mat(
      r_mat, data_list$Np, tmp_acov
    )
  }

  return(data_list)
}
