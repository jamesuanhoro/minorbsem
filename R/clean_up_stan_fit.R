#' Clean up Stan fit helper function
#' @description A function that cleans up model returned by Stan
#' @param stan_fit Stan fit
#' @param data_list Data list object passed to Stan
#' @returns An object of \code{\link{mbsem-class}}
#' @keywords internal
clean_up_stan_fit <- function(
    stan_fit,
    data_list) {
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

  major_parameters <- as.data.frame(rbind(
    rms_result, coef_result, load_result, phi_result, pv_result,
    rv_result, rc_result
  ))

  major_parameters$ess_bulk <- round(major_parameters$ess_bulk, 1)
  major_parameters$ess_tail <- round(major_parameters$ess_tail, 1)
  major_parameters <- major_parameters[
    ,
    c(
      "group", "from", "op", "to",
      "mean", "median",
      "sd", "mad", "q5", "q95",
      "rhat", "ess_bulk", "ess_tail"
    )
  ]

  minor_factor_matrix <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws(stan_fit, variable = "Resid")
  ))

  mbsem_result <- new_mbsem()
  mbsem_result <- methods::initialize(
    mbsem_result,
    major_parameters = major_parameters,
    minor_factor_matrix = minor_factor_matrix,
    data_list = data_list,
    stan_fit = stan_fit,
    version = minorbsem_version()
  )

  return(mbsem_result)
}
