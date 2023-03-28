#' Pre-plotting Stan fit helper function
#' @description A function that cleans up model results to prepare for plotting
#' @param object (mbsem) An object of \code{\link{mbsem-class}}
#' returned by \code{\link{minorbsem}}.
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

  rmsea_result <- matrix(nrow = nrow(rms_result), ncol = 0)
  if (data_list$meta == 1) {
    rmsea_result <- posterior::as_draws_df(
      posterior::as_draws(stan_fit, variable = "rmsea_mn")
    )
    colnames(rmsea_result)[1] <- "rm: rmsea"
  }

  load_idxs <- paste0("Load_mat[", apply(which(
    data_list$loading_pattern >= ifelse(data_list$complex_struc == 1, -999, 1),
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
  rsq_result <- matrix(nrow = nrow(rv_result), ncol = 0)
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
    rsq_result <- posterior::as_draws_df(
      posterior::as_draws(stan_fit, variable = "r_square")
    )
    rsq_result <- rename_post_df_columns(
      rsq_result, factor_labels, factor_labels, "rsq:", "~~",
      TRUE, "r_square\\[|\\]", "r_square\\[|\\]"
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
    rmsea_result, rms_result, coef_result, load_result, phi_result,
    rsq_result, rv_result, rc_result, resid_result
  )

  # ID first set of chain/iter/draw
  keep_cols <- c(
    which(colnames(result) == ".chain")[1],
    which(colnames(result) == ".iteration")[1],
    which(colnames(result) == ".draw")[1]
  )

  # Drop duplicated columns except for chain/iter/draw
  dup_list <- duplicated(as.list(result))
  # Drop unchanging columns
  no_change <- apply(result, 2, stats::var) < .Machine$double.eps
  dup_change_list <- dup_list | no_change
  dup_change_list[keep_cols] <- FALSE
  result <- result[!dup_change_list]

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
    c("rm", "co", "lo", "fc", "rsq", "ev", "rc", "re")
  )

  return(result_long)
}
