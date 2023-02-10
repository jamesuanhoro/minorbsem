#' Pre-plotting Stan fit helper function
#' @description A function that cleans up model results to prepare for plotting
#' @param object (mbsem) An object of class mbsem
#' returned by minorbsem.
#' @returns A data.frame that is ready for plotting
#' @keywords internal
prepare_stan_plot_data <- function(object) {
  result <- list()

  data_list <- object@data_list
  stan_fit <- object@stan_fit

  indicator_labels <- rownames(data_list$loading_pattern)
  factor_labels <- colnames(data_list$loading_pattern)

  rms_result <- posterior::as_draws_df(
    posterior::as_draws(stan_fit, variable = "rms_src")
  )
  colnames(rms_result)[1] <- "rm: rms_src"

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
  result <- result[apply(result, 2, stats::var) > 0]

  varying <- colnames(result)[
    which(regexpr("\\:", colnames(result)) > 0)
  ]

  result_long <- stats::reshape(
    result,
    varying = list(varying),
    times = varying,
    idvar = paste0(".", c("chain", "iteration", "draw")),
    timevar = "parameter", v.names = "value",
    direction = "long"
  )
  rownames(result_long) <- NULL

  result_long$param_class <- unlist(lapply(strsplit(
    result_long$parameter, ":"
  ), "[[", 1))

  result_long$parameter <- factor(
    result_long$parameter,
    unique(result_long$parameter)
  )
  result_long$param_class <- factor(
    result_long$param_class,
    c("rm", "co", "lo", "fc", "fv", "ev", "rc", "re")
  )

  return(result_long)
}
