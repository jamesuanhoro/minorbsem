# Helper functions in package

#' Package on attach message
#' @returns R version minimum as string
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  version <- minorbsem_version()
  packageStartupMessage(" ")
  packageStartupMessage(strrep("#", 79))
  packageStartupMessage("This is ", paste(pkgname, version))
  packageStartupMessage(
    "\nAll users of R (or SEM) are invited to report bugs, ",
    "submit functions or ideas\nfor functions. ",
    "An efficient way to do this is to open an issue on GitHub\n",
    "https://github.com/jamesuanhoro/minorbsem/issues/."
  )
  packageStartupMessage(strrep("#", 79))
}

#' Get package version function
#' @returns Package version as string
#' @keywords internal
minorbsem_version <- function() {
  version <- read.dcf(
    system.file("DESCRIPTION", package = "minorbsem"),
    fields = "Version"
  )[1]
  return(version)
}

#' Get R minimum function
#' @returns R version minimum as string
#' @keywords internal
minimum_r_version <- function() {
  r_version_info <- read.dcf(
    system.file("DESCRIPTION", package = "minorbsem"),
    fields = "Depends"
  )[1]
  r_version <- trimws(gsub("R \\(>=|\\)", "", r_version_info))
  return(r_version)
}

#' Run generic CFA and SEM function
#' @returns NULL
#' @export
init_minorbsem <- function() {
  model_syntax <- "
  F1 =~ x1 + x2 + x3\n F2 =~ x4 + x5 + x6\n F3 =~ x7 + x8 + x9"
  minorbsem(
    model_syntax, minorbsem::HS,
    warmup = 200, sampling = 200, chains = 1, show = FALSE,
    method = "none", refresh = 0, show_messages = FALSE
  )
  model_syntax <- "
  ind60 =~ x1 + x2 + x3\n dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8\n dem60 ~ ind60\n dem65 ~ ind60 + dem60"
  minorbsem(
    model_syntax, minorbsem::PD,
    warmup = 200, sampling = 200, chains = 1, show = FALSE,
    method = "none", refresh = 0, show_messages = FALSE
  )
  return(NULL)
}

#' Check user input function
#' @description A function that checks user input for adequacy
#' and fails on inadequate input.
#' @param type string for type of check
#' @param object_1 Object to check
#' @param object_2 Object to check
#' @param object_3 Object to check
#' @param object_4 Object to check
#' @returns NULL
#' @keywords internal
user_input_check <- function(
    type,
    object_1 = NULL,
    object_2 = NULL,
    object_3 = NULL,
    object_4 = NULL) {
  if (type == "model" && is.null(object_1)) {
    stop("Model cannot be null")
  }

  if (type == "priors") {
    if (!inherits(object_1, "mbsempriors")) {
      stop("See ?new_mbsempriors for how to set up priors.")
    }
  }

  if (grepl("method", type)) {
    accepted_methods <- method_hash()
    if (!tolower(object_1) %in% tolower(accepted_methods)) {
      err_msg <- paste0(
        "method must be one of the following: ",
        paste0("\"", accepted_methods, "\"", collapse = ", ")
      )
      stop(err_msg)
    }
  }

  if (type == "data") {
    if (is.null(object_1) && (is.null(object_2) || is.null(object_3))) {
      stop(paste0(
        "User must provide either:\n\t",
        "(i) a dataset or\n\t",
        "(ii) sample covariance and sample size"
      ))
    }
  }

  return(NULL)
}

#' target fitter function
#' @description A function that takes user input and fits the
#' Stan model.
#' @param data_list Data list object passed to Stan
#' @inheritParams minorbsem
#' @returns Fitted Stan model
#' @keywords internal
target_fitter <- function(
    data_list,
    seed,
    warmup,
    sampling,
    refresh,
    adapt_delta,
    max_treedepth,
    chains,
    ncores,
    show_messages) {
  init_resid <- function() {
    list(
      rms_src_p = array(.025, (data_list$method != 100) * 1),
      resids = array(
        0,
        (data_list$method < 90) * (data_list$Ni^2 - data_list$Ni) / 2
      ),
      loadings_complex = array(
        0,
        data_list$complex_struc * sum(
          data_list$loading_pattern == 0 & data_list$loading_fixed == -999
        )
      )
    )
  }

  if (data_list$sem_indicator == 0) {
    mod_resid <- instantiate::stan_package_model(
      name = "cfa_resid", package = "minorbsem"
    )
  } else if (data_list$sem_indicator == 1) {
    mod_resid <- instantiate::stan_package_model(
      name = "sem_resid", package = "minorbsem"
    )
  }

  message("Fitting Stan model ...")

  stan_fit <- mod_resid$sample(
    data = data_list,
    seed = seed,
    iter_warmup = warmup,
    iter_sampling = sampling,
    refresh = refresh,
    init = init_resid,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    chains = chains,
    parallel_chains = ncores,
    show_messages = show_messages
  )

  return(stan_fit)
}

#' Select random method and random case function
#'
#' @returns randomly selected method and case
#' @keywords internal
random_method_selection <- function() {
  accepted_methods <- method_hash()

  method <- sample(accepted_methods, 1)

  case_fun <- sample(1:2, 1)
  if (case_fun == 1) {
    method <- tolower(method)
  } else if (case_fun == 2) {
    method <- toupper(method)
  }

  print(method)
  return(method)
}

#' converter function
#' @description A function that swaps type from string to integer
#' and vice-versa
#' @param search_term Integer or string or NULL
#' @param str_list Character vector of accepted values
#' @returns If search_term is integer, returns string and vice-versa
#' @keywords internal
converter_helper <- function(search_term, str_list) {
  if (is.null(search_term)) {
    converted_value <- names(str_list)
  } else if (is.integer(search_term) || is.numeric(search_term)) {
    search_term <- as.integer(search_term)
    converted_value <- names(str_list)[which(search_term == str_list)]
  } else if (is.character(search_term)) {
    idx <- which(tolower(search_term) == tolower(names(str_list)))
    converted_value <- as.integer(str_list[idx])
  }

  return(converted_value)
}

#' Method hash function
#' @description A function that swaps method from string to integer
#' and vice-versa
#' @inheritParams converter_helper
#' @returns If search_term is integer, returns string and vice-versa
#' @keywords internal
method_hash <- function(search_term = NULL) {
  # Reserving 90+ for methods with no residuals
  list_methods <- c(
    "normal" = 1,
    "lasso" = 2,
    "logistic" = 3,
    "GDP" = 4,
    "WB" = 90,
    "WB-cond" = 91,
    "WW" = 92,
    "none" = 100
  )
  converted_value <- converter_helper(search_term, list_methods)
  return(converted_value)
}

#' Posterior summary helper function
#' @description A function that slightly modifies the default
#' summary function in posterior package
#' @param stan_fit Fitted Stan object
#' @param variable Variable(s) to search for in Stan
#' @param interval Confidence interval to select
#' @param major If TRUE, add some preamble for printing the major
#' parameters table.
#' @returns Summary of posterior draws
#' @keywords internal
mbsem_post_sum <- function(stan_fit, variable, interval = .9, major = FALSE) {
  draws <- posterior::as_draws(stan_fit$draws(variable))

  lo_lim <- (1.0 - interval) / 2.0
  up_lim <- 1.0 - lo_lim # nolint
  sum_stats <- posterior::summarise_draws(
    draws, "mean", "median", "sd", "mad",
    ~ quantile(.x, probs = c(lo_lim, up_lim), na.rm = TRUE)
  )
  convergence_metrics <- posterior::summarise_draws(
    draws, posterior::default_convergence_measures()
  )

  result <- as.data.frame(merge(
    sum_stats, convergence_metrics,
    by = "variable"
  ))

  var_order <- sum_stats$variable
  # Without this, the row order within param type changes after merging
  result <- result[match(var_order, result$variable), ]

  if (isTRUE(major)) {
    result <- cbind(
      variable = result$variable,
      group = "", from = "", op = "", to = "",
      result[, -1]
    )
  }

  return(result)
}

#' Create parameter list for plotting
#' @param data_list Data list object passed to Stan
#' @returns Name of varying model parameters
#' @keywords internal
get_param_plot_list <- function(data_list) {
  param_structure <- data_list$loading_pattern
  fac_names <- colnames(param_structure)
  ind_names <- rownames(param_structure)

  rms_params <- array(dim = 0)
  if (data_list$method < 90) {
    rms_params <- c("RMSE" = "rms_src")
  }

  coef_idxs <- matrix(nrow = 0, ncol = 2)
  rsq_idxs <- array(dim = 0)
  if (data_list$sem_indicator == 0) {
    phi_idxs <- which(
      lower.tri(data_list$corr_mask) & data_list$corr_mask == 1,
      arr.ind = TRUE
    )
    load_idxs <- which(
      data_list$loading_pattern >=
        ifelse(data_list$complex_struc == 1, -999, 1) &
        data_list$loading_fixed == -999,
      arr.ind = TRUE
    )
    rv_idxs <- which(data_list$res_var_pattern != 0)
  } else if (data_list$sem_indicator == 1) {
    phi_idxs <- data_list$F_corr_mat
    load_idxs <- which(data_list$loading_pattern >= 1, arr.ind = TRUE)
    rv_idxs <- seq_len(data_list$Ni)
    coef_idxs <- which(
      data_list$coef_pattern == 1,
      arr.ind = TRUE
    )
    rsq_idxs <- which(rowSums(data_list$coef_pattern) >= 1)
  }

  coef_params <- array(dim = 0)
  if (nrow(coef_idxs) > 0) {
    coef_params <- paste0("Coef_mat[", apply(
      coef_idxs, 1, paste0,
      collapse = ","
    ), "]")
    names(coef_params) <- apply(coef_idxs, 1, function(x) {
      paste0(fac_names[x[2]], "~", fac_names[x[1]])
    })
  }

  rsq_params <- array(dim = 0)
  if (length(rsq_idxs) > 0) {
    rsq_params <- paste0("r_square[", rsq_idxs, "]")
    names(rsq_params) <- paste0("rsq:", fac_names[rsq_idxs])
  }

  phi_params <- array(dim = 0)
  if (nrow(phi_idxs) > 0 && data_list$sem_indicator == 0) {
    phi_params <- paste0("phi_mat[", apply(
      phi_idxs, 1, paste0,
      collapse = ","
    ), "]")
    names(phi_params) <- apply(phi_idxs, 1, function(x) {
      paste0(fac_names[x[1]], "~~", fac_names[x[2]])
    })
  } else if (nrow(phi_idxs) > 0 && data_list$sem_indicator == 1) {
    phi_params <- paste0("phi_cor[", seq_len(data_list$Nf_corr), "]")
    names(phi_params) <- apply(phi_idxs, 1, function(x) {
      paste0(fac_names[x[1]], "~~", fac_names[x[2]])
    })
  }

  load_params <- array(dim = 0)
  if (nrow(load_idxs) > 0) {
    load_params <- paste0(
      "Load_mat[", apply(load_idxs, 1, paste0, collapse = ","), "]"
    )
    names(load_params) <- apply(load_idxs, 1, function(x) {
      paste0(fac_names[x[2]], "=~", ind_names[x[1]])
    })
  }

  rv_params <- array(dim = 0)
  if (length(rv_idxs) > 0) {
    rv_params <- paste0("res_var[", rv_idxs, "]")
    names(rv_params) <- paste0(ind_names[rv_idxs], "~~", ind_names[rv_idxs])
  }

  rc_params <- array(dim = 0)
  if (data_list$Nce > 0) {
    rc_idxs <- data_list$error_mat
    rc_params <- paste0("res_cor[", seq_len(data_list$Nce), "]")
    names(rc_params) <- apply(rc_idxs, 1, function(x) {
      paste0(ind_names[x[1]], "~~", ind_names[x[2]])
    })
  }

  params <- c(
    rms_params, coef_params, rsq_params, phi_params, load_params,
    rv_params, rc_params
  )

  return(params)
}

#' Modify major parameters table helper function
#' @description A function that adds user friendly descriptions to the
#' major parameters table
#' @param major_parameters Major paramters table
#' @param idxs Relevant rows indexes
#' @param group Parameter group
#' @param op lavaan style operator
#' @param from Variable from
#' @param to Variable to
#' @returns Updated major parameters table
#' @keywords internal
modify_major_params <- function(
    major_parameters,
    idxs,
    group = "",
    op = "",
    from = "",
    to = "") {
  result <- major_parameters

  if (length(idxs) > 0) {
    result[idxs, ]$group <- group
    result[idxs, ]$op <- op
    result[idxs, ]$from <- from
    result[idxs, ]$to <- to
  }

  return(result)
}

#' Create major parameters helper function
#' @description A function that creates the table of major parameters
#' @param stan_fit Fitted Stan object
#' @param data_list Data list object passed to Stan
#' @param interval Confidence interval to select
#' @returns Summary of posterior draws
#' @keywords internal
create_major_params <- function(stan_fit, data_list, interval = .9) {
  indicator_labels <- rownames(data_list$loading_pattern)
  factor_labels <- colnames(data_list$loading_pattern)

  params <- c("ppp", "rms_src")
  from_list <- c("PPP", "RMSE")

  load_idxs <- paste0("Load_mat[", apply(which(
    data_list$loading_pattern >= ifelse(data_list$complex_struc == 1, -999, 1),
    arr.ind = TRUE
  ), 1, paste0, collapse = ","), "]")
  params <- c(params, load_idxs)

  if (data_list$Nce > 0) {
    params <- c(params, "res_cor")
  }

  if (data_list$sem_indicator == 0) {
    params <- c(params, "res_var")

    params <- c(params, "phi_mat")
  } else if (data_list$sem_indicator == 1) {
    # Get interfactor correlations
    if (data_list$Nf_corr > 0) {
      params <- c(params, "phi_cor")
    }

    # Get R-square
    params <- c(params, "r_square")

    # Get factor coefficients
    coef_idxs <- paste0("Coef_mat[", apply(which(
      data_list$coef_pattern == 1,
      arr.ind = TRUE
    ), 1, paste0, collapse = ","), "]")
    params <- c(params, coef_idxs)
  }

  major_parameters <- mbsem_post_sum(
    stan_fit = stan_fit, variable = params, interval = interval, major = TRUE
  )

  # Dump duplicates CFA inter-factor correlations here
  phi_rows_idx <- grepl("phi\\_mat", major_parameters$variable)
  if (sum(phi_rows_idx) > 0) {
    phi_rows <- major_parameters[phi_rows_idx, , drop = FALSE]
    phi_rows$ind_1 <- as.integer(gsub(
      "phi\\_mat\\[|,\\d+\\]", "", phi_rows$variable
    ))
    phi_rows$ind_2 <- as.integer(gsub(
      "phi\\_mat\\[\\d+,|\\]", "", phi_rows$variable
    ))
    dump_rows <- phi_rows$ind_1 >= phi_rows$ind_2
    phi_rows_idx[phi_rows_idx == TRUE] <- dump_rows
    major_parameters <- major_parameters[!phi_rows_idx, , drop = FALSE]
  }

  mid_cols <- which(
    colnames(major_parameters) %in% c("median", "sd", "mad") |
      grepl("\\%", colnames(major_parameters))
  )
  major_parameters[
    major_parameters$variable == "ppp",
    mid_cols
  ] <- NA_real_

  major_parameters <- modify_major_params(
    major_parameters,
    which(major_parameters$variable %in% c("ppp", "rmsea_mn")),
    group = "Goodness of fit",
    from = from_list[1]
  )
  major_parameters <- modify_major_params(
    major_parameters,
    which(major_parameters$variable == "rms_src"),
    group = "Goodness of fit",
    from = from_list[2]
  )

  idxs <- which(grepl("Load\\_mat", major_parameters$variable))
  major_parameters <- modify_major_params(
    major_parameters, idxs,
    group = "Factor loadings", op = "=~",
    from = factor_labels[as.integer(
      gsub("Load_mat\\[\\d+,|\\]", "", major_parameters[idxs, ]$variable)
    )],
    to = indicator_labels[as.integer(
      gsub("Load_mat\\[|,\\d+\\]", "", major_parameters[idxs, ]$variable)
    )]
  )

  idxs <- which(grepl("res\\_cor", major_parameters$variable))
  major_parameters <- modify_major_params(
    major_parameters, idxs,
    group = "Error correlations", op = "~~",
    from = indicator_labels[data_list$error_mat[, 1]],
    to = indicator_labels[data_list$error_mat[, 2]]
  )

  idxs <- which(grepl("res\\_var", major_parameters$variable))
  major_parameters <- modify_major_params(
    major_parameters, idxs,
    group = "Residual variances", op = "~~",
    from = indicator_labels[as.integer(
      gsub("res_var\\[|\\]", "", major_parameters[idxs, ]$variable)
    )],
    to = indicator_labels[as.integer(
      gsub("res_var\\[|\\]", "", major_parameters[idxs, ]$variable)
    )]
  )

  idxs <- which(grepl("phi\\_mat", major_parameters$variable))
  major_parameters <- modify_major_params(
    major_parameters, idxs,
    group = "Inter-factor correlations", op = "~~",
    from = factor_labels[as.integer(
      gsub("phi_mat\\[\\d+,|\\]", "", major_parameters[idxs, ]$variable)
    )],
    to = factor_labels[as.integer(
      gsub("phi_mat\\[|,\\d+\\]", "", major_parameters[idxs, ]$variable)
    )]
  )
  major_parameters <- major_parameters[
    !(major_parameters$group == "Inter-factor correlations" &
      major_parameters$from == major_parameters$to),
  ]

  idxs <- which(grepl("phi\\_cor", major_parameters$variable))
  major_parameters <- modify_major_params(
    major_parameters, idxs,
    group = "Inter-factor correlations", op = "~~",
    from = factor_labels[data_list$F_corr_mat[, 1]],
    to = factor_labels[data_list$F_corr_mat[, 2]]
  )

  idxs <- which(grepl("r\\_square", major_parameters$variable))
  major_parameters <- modify_major_params(
    major_parameters, idxs,
    group = "R square", op = "~~",
    from = factor_labels[as.integer(
      gsub("r_square\\[|\\]", "", major_parameters[idxs, ]$variable)
    )],
    to = factor_labels[as.integer(
      gsub("r_square\\[|\\]", "", major_parameters[idxs, ]$variable)
    )]
  )

  idxs <- which(grepl("Coef\\_mat", major_parameters$variable))
  major_parameters <- modify_major_params(
    major_parameters, idxs,
    group = "Latent regression coefficients", op = "~",
    from = factor_labels[as.integer(
      gsub("Coef_mat\\[|,\\d+\\]", "", major_parameters[idxs, ]$variable)
    )],
    to = factor_labels[as.integer(
      gsub("Coef_mat\\[\\d+,|\\]", "", major_parameters[idxs, ]$variable)
    )]
  )

  major_parameters$ess_bulk <- round(major_parameters$ess_bulk, 1)
  major_parameters$ess_tail <- round(major_parameters$ess_tail, 1)

  target <- c(
    "Goodness of fit",
    "Dispersion between and within clusters",
    "Latent regression coefficients", "R square",
    "Factor loadings", "Inter-factor correlations",
    "Residual variances", "Error correlations"
  )

  new_order <- unlist(sapply(target, function(x) {
    which(major_parameters$group == x)
  }))
  major_parameters <- major_parameters[new_order, -1]
  new_order <- unlist(sapply(c("PPP", "RMSEA", "RMSE"), function(x) {
    which(major_parameters$from[1:2] == x)
  }))
  major_parameters[1:2, ] <- major_parameters[new_order, ]

  return(major_parameters)
}

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
  col_names <- paste0(begin_name, " ", part_1, operation, part_2)
  colnames(df)[1:len_vars] <- col_names

  return(df)
}
