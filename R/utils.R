# Helper functions in package

#' Stan data helper function
#' @description A function that creates data list object passed to Stan
#' @param lavaan_object lavaan fit object of corresponding model
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
#' @export
create_data_list <- function(
    lavaan_object = NULL,
    lkj_shape = 2.0,
    sl_par = 1.0,
    rs_par = 2.5,
    rc_par = 2.0,
    sc_par = 1.0) {
  data_list <- list()

  # Retrieve parameter structure from lavaan
  param_structure <- lavaan::lavInspect(lavaan_object)

  # Has data?
  data_list$Y <- lavaan_object@Data@X[[1]]
  data_list$has_data <- ifelse(is.null(data_list$Y), 0, 1)

  # Set up priors
  if (any(c(lkj_shape, sl_par, rs_par, rc_par, sc_par) <= 0)) {
    stop("lkj_shape, sl_par, rs_par, rc_par, sc_par must all exceed 0")
  }
  data_list$shape_phi_c <- lkj_shape # Shape parameter for LKJ of interfactor corr
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
  Psi <- param_structure$psi
  sum_off_diag_psi <- sum(Psi[lower.tri(Psi)])
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
      data_list$F_corr_mat <- which(Psi != 0, arr.ind = TRUE)
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

#' Clean up Stan fit helper function
#' @description A function that cleans up model returned by Stan
#' @param stan_fit Stan fit
#' @param data_list Data list object passed to Stan
#' @returns A list of major structural parameters,
#' minor factors standardized residual covariance matrix, and original Stan fit
#' @keywords internal
#' @export
clean_up_stan_fit <- function(
    stan_fit,
    data_list) {
  result <- list()

  indicator_labels <- rownames(data_list$loading_pattern)
  factor_labels <- colnames(data_list$loading_pattern)

  rms_result <- posterior::summarise_draws(
    posterior::as_draws(stan_fit, variable = "rms_src")
  )
  rms_result$group <- "RMSE(resids)"
  rms_result$op <- ""
  rms_result$from <- ""
  rms_result$to <- ""

  load_idxs <- paste0("Load_mat[", apply(which(
    data_list$loading_pattern == 1,
    arr.ind = TRUE
  ), 1, paste0, collapse = ","), "]")
  load_result <- posterior::summarise_draws(
    posterior::as_draws(stan_fit, variable = load_idxs)
  )
  load_result$group <- "Factor loadings"
  load_result$op <- "=~"
  load_result$from <- factor_labels[as.integer(
    gsub("Load_mat\\[\\d+,|\\]", "", load_result$variable)
  )]
  load_result$to <- indicator_labels[as.integer(
    gsub("Load_mat\\[|,\\d+\\]", "", load_result$variable)
  )]

  rv_result <- posterior::summarise_draws(
    posterior::as_draws(stan_fit, variable = "res_var")
  )
  rv_result$group <- "Residual variances"
  rv_result$op <- "~~"
  rv_result$from <- indicator_labels[as.integer(
    gsub("res_var\\[|\\]", "", rv_result$variable)
  )]
  rv_result$to <- indicator_labels[as.integer(
    gsub("res_var\\[|\\]", "", rv_result$variable)
  )]

  rc_result <- matrix(nrow = 0, ncol = nrow(rv_result))
  if (data_list$Nce > 0) {
    rc_result <- posterior::summarise_draws(
      posterior::as_draws(stan_fit, variable = "res_cor")
    )
    rc_result$group <- "Error correlations"
    rc_result$op <- "~~"
    rc_result$from <- indicator_labels[data_list$error_mat[, 1]]
    rc_result$to <- indicator_labels[data_list$error_mat[, 2]]
  }

  coef_result <- matrix(nrow = 0, ncol = nrow(rv_result))
  phi_result <- matrix(nrow = 0, ncol = nrow(rv_result))
  pv_result <- matrix(nrow = 0, ncol = nrow(rv_result))
  if (data_list$sem_indicator == 0) {
    phi_result <- posterior::summarise_draws(
      posterior::as_draws(stan_fit, variable = "phi_mat")
    )
    phi_duplicates <- duplicated(phi_result[, c("mean", "median", "sd")])
    phi_result <- phi_result[!phi_duplicates, ]
    phi_result$group <- "Inter-factor correlations"
    phi_result$op <- "~~"
    phi_result$from <- factor_labels[as.integer(
      gsub("phi_mat\\[\\d+,|\\]", "", phi_result$variable)
    )]
    phi_result$to <- factor_labels[as.integer(
      gsub("phi_mat\\[|,\\d+\\]", "", phi_result$variable)
    )]
    phi_result <- phi_result[phi_result$from != phi_result$to, ]
  } else if (data_list$sem_indicator == 1) {
    # Get interfactor correlations
    if (data_list$Nf_corr > 0) {
      phi_result <- posterior::summarise_draws(
        posterior::as_draws(stan_fit, variable = "phi_cor")
      )
      phi_result$group <- "Inter-factor correlations"
      phi_result$op <- "~~"
      phi_result$from <- factor_labels[data_list$F_corr_mat[, 1]]
      phi_result$to <- factor_labels[data_list$F_corr_mat[, 2]]
    }

    # Get factor variances
    pv_result <- posterior::summarise_draws(
      posterior::as_draws(stan_fit, variable = "phi_var")
    )
    pv_result$group <- "Factor variances"
    pv_result$op <- "~~"
    pv_result$from <- factor_labels[as.integer(
      gsub("phi_var\\[|\\]", "", pv_result$variable)
    )]
    pv_result$to <- factor_labels[as.integer(
      gsub("phi_var\\[|\\]", "", pv_result$variable)
    )]

    # Get factor coefficients
    coef_idxs <- paste0("Coef_mat[", apply(which(
      data_list$coef_pattern == 1,
      arr.ind = TRUE
    ), 1, paste0, collapse = ","), "]")
    coef_result <- posterior::summarise_draws(
      posterior::as_draws(stan_fit, variable = coef_idxs)
    )
    coef_result$group <- "Latent regression coefficients"
    coef_result$op <- "~"
    coef_result$from <- factor_labels[as.integer(
      gsub("Coef_mat\\[|,\\d+\\]", "", coef_result$variable)
    )]
    coef_result$to <- factor_labels[as.integer(
      gsub("Coef_mat\\[\\d+,|\\]", "", coef_result$variable)
    )]
  }

  result$major_parameters <- rbind(
    rms_result, coef_result, load_result, phi_result, pv_result,
    rv_result, rc_result
  )

  result$major_parameters$ess_bulk <- round(result$major_parameters$ess_bulk, 1)
  result$major_parameters$ess_tail <- round(result$major_parameters$ess_tail, 1)
  result$major_parameters <- result$major_parameters[
    ,
    c(
      "group", "from", "op", "to",
      "mean", "median",
      "sd", "mad", "q5", "q95",
      "rhat", "ess_bulk", "ess_tail"
    )
  ]

  result$minor_factor_matrix <- posterior::summarise_draws(
    posterior::as_draws(stan_fit, variable = "Resid")
  )

  result$data_list <- data_list
  result$stan_fit <- stan_fit

  return(result)
}

#' Add row header helper function
#' @description A function that adds row headers to kable object
#' @param kbl_object Kable object
#' @param table_to_print Table used in Kable object
#' @param search_term Label to search for in table and row header text
#' @param extra Extra text to print in row header
#' @returns If search_term present in term, modified Kable object,
#' otherwise: original object is returned
#' @keywords internal
#' @export
add_row_header <- function(kbl_object, table_to_print, search_term, extra = NULL) {
  if (any(table_to_print[1] == search_term)) {
    kbl_object <- kbl_object |>
      kableExtra::pack_rows(
        paste0(search_term, " ", extra),
        which(table_to_print[1] == search_term)[1],
        rev(which(table_to_print[1] == search_term))[1]
      )
  }
  return(kbl_object)
}

# load_result <- rename_post_df_columns(
#   load_result, factor_labels, indicator_labels, "lo:", "=~",
#   TRUE, "Load_mat\\[\\d+,|\\]", "Load_mat\\[|,\\d+\\]")

#' Rename columns of posterior data.frame prior by parameter type
#'
#' @param df Posterior data.frame
#' @param labels_1 Labels for first variable
#' @param labels_2 Labels for second variable
#' @param begin_name String to begin columns names with
#' @param operation Operation, one of =~, ~, ~~
#' @param search_replace TRUE/FALSE for type of column naming
#' @param search_term_1 Terms to use in replacement for first variable
#' @param search_term_2 Terms to use in replacement for second variable
#' @returns A renamed posterior data.frame
#' @keywords internal
#' @export
rename_post_df_columns <- function(
    df,
    labels_1,
    labels_2,
    begin_name,
    operation,
    search_replace,
    search_term_1,
    search_term_2) {
  col_names <- colnames(df)
  len_vars <- length(col_names) - 3
  col_names <- col_names[1:len_vars]
  part_1 <- vector("character")
  part_2 <- vector("character")
  if (search_replace == TRUE) {
    part_1 <- labels_1[as.integer(
      gsub(search_term_1, "", col_names)
    )]
    part_2 <- labels_2[as.integer(
      gsub(search_term_2, "", col_names)
    )]
  } else {
    part_1 <- labels_1[search_term_1]
    part_2 <- labels_2[search_term_2]
  }
  col_names <- paste0(begin_name, part_1, operation, part_2)
  colnames(df)[1:len_vars] <- col_names

  return(df)
}

#' Pre-plotting Stan fit helper function
#' @description A function that cleans up model results to prepare for plotting
#' @param stan_fit Stan fit
#' @param data_list Data list object passed to Stan
#' @returns A data.frame that is ready for plotting
#' @keywords internal
#' @export
prepare_stan_plot <- function(
    stan_fit,
    data_list) {
  result <- list()

  indicator_labels <- rownames(data_list$loading_pattern)
  factor_labels <- colnames(data_list$loading_pattern)

  rms_result <- posterior::as_draws_df(
    posterior::as_draws(stan_fit, variable = "rms_src")
  )
  colnames(rms_result)[1] <- "rm:rms_src"

  load_idxs <- paste0("Load_mat[", apply(which(
    data_list$loading_pattern == 1,
    arr.ind = TRUE
  ), 1, paste0, collapse = ","), "]")
  load_result <- posterior::as_draws_df(
    posterior::as_draws(stan_fit, variable = load_idxs)
  )
  load_result <- rename_post_df_columns(
    load_result, factor_labels, indicator_labels, "lo:", "=~",
    TRUE, "Load_mat\\[\\d+,|\\]", "Load_mat\\[|,\\d+\\]"
  )

  rv_result <- posterior::as_draws_df(
    posterior::as_draws(stan_fit, variable = "res_var")
  )
  rv_result <- rename_post_df_columns(
    rv_result, indicator_labels, indicator_labels, "ev:", "~~",
    TRUE, "res_var\\[|\\]", "res_var\\[|\\]"
  )

  rc_result <- matrix(nrow = nrow(rv_result), ncol = 0)
  if (data_list$Nce > 0) {
    rc_result <- posterior::as_draws_df(
      posterior::as_draws(stan_fit, variable = "res_cor")
    )
    rc_result <- rename_post_df_columns(
      rc_result, indicator_labels, indicator_labels, "rc:", "~~",
      FALSE, data_list$error_mat[, 1], data_list$error_mat[, 2]
    )
  }

  coef_result <- matrix(nrow = nrow(rv_result), ncol = 0)
  phi_result <- matrix(nrow = nrow(rv_result), ncol = 0)
  pv_result <- matrix(nrow = nrow(rv_result), ncol = 0)
  if (data_list$sem_indicator == 0) {
    phi_result <- posterior::as_draws_df(
      posterior::as_draws(stan_fit, variable = "phi_mat")
    )
    phi_duplicates <- duplicated(as.list(phi_result))
    phi_result <- phi_result[!phi_duplicates]
    phi_result <- rename_post_df_columns(
      phi_result, factor_labels, factor_labels, "fc:", "~~",
      TRUE, "phi_mat\\[\\d+,|\\]", "phi_mat\\[|,\\d+\\]"
    )
  } else if (data_list$sem_indicator == 1) {
    # Get interfactor correlations
    if (data_list$Nf_corr > 0) {
      phi_result <- posterior::as_draws_df(
        posterior::as_draws(stan_fit, variable = "phi_cor")
      )
      phi_result <- rename_post_df_columns(
        phi_result, factor_labels, factor_labels, "fc:", "~~",
        FALSE, data_list$F_corr_mat[, 1], data_list$F_corr_mat[, 2]
      )
    }

    # Get factor variances
    pv_result <- posterior::as_draws_df(
      posterior::as_draws(stan_fit, variable = "phi_var")
    )
    pv_result <- rename_post_df_columns(
      pv_result, factor_labels, factor_labels, "fv:", "~~",
      TRUE, "phi_var\\[|\\]", "phi_var\\[|\\]"
    )

    # Get factor coefficients
    coef_idxs <- paste0("Coef_mat[", apply(which(
      data_list$coef_pattern == 1,
      arr.ind = TRUE
    ), 1, paste0, collapse = ","), "]")
    coef_result <- posterior::as_draws_df(
      posterior::as_draws(stan_fit, variable = coef_idxs)
    )
    coef_result <- rename_post_df_columns(
      coef_result, factor_labels, factor_labels, "co:", "~",
      TRUE, "Coef_mat\\[|,\\d+\\]", "Coef_mat\\[\\d+,|\\]"
    )
  }

  resid_result <- posterior::as_draws_df(
    posterior::as_draws(stan_fit, variable = "Resid")
  )
  resid_result <- rename_post_df_columns(
    resid_result, indicator_labels, indicator_labels, "re:", "~~",
    TRUE, "Resid\\[|,\\d+\\]", "Resid\\[\\d+,|\\]"
  )

  result <- cbind(
    rms_result, coef_result, load_result, phi_result, pv_result,
    rv_result, rc_result, resid_result
  )

  # Drop duplicated columns
  result <- result[!duplicated(as.list(result))]
  # Drop unchanging columns
  result <- result[apply(result, 2, var) > 0]

  return(result)
}
