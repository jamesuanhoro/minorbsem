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
    "All users of R (or SEM) are invited to submit functions ",
    "or ideas for functions."
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

  if (regexpr("method", type) > 0) {
    accepted_methods <- method_hash()
    if (type == "method-meta") {
      accepted_methods <- accepted_methods[
        which(regexpr("WB", accepted_methods) < 0)
      ]
    }
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

  if (type == "data-meta") {
    if (
      (is.null(object_1) || is.null(object_2)) &&
        (is.null(object_3) || is.null(object_4))
    ) {
      stop(paste0(
        "User must provide either:\n\t",
        "(i) a dataset and group variable or\n\t",
        "(ii) sample covariance and sample size"
      ))
    }
  }

  if (type == "type-meta") {
    accepted_types <- type_hash()
    err_msg <- paste0(
      "type must be one of the following: ",
      paste0(
        "\"", accepted_types, "\"", collapse = ", "
      )
    )
    if (is.null(object_1)) stop(err_msg)
    if (!tolower(object_1) %in% accepted_types) stop(err_msg)
  }

  return(NULL)
}

#' Fill in missing values in sample covs function
#' @inheritParams meta_mbcfa
#' @returns list of filled in sample covs
#' @keywords internal
fill_in_missing_covs <- function(sample_cov) {
  sample_cov <- lapply(seq_along(sample_cov), function(i) {
    s_mat <- sample_cov[[i]]
    diag(s_mat)[is.na(diag(s_mat))] <- 1
    s_mat[is.na(s_mat)] <- 999
    s_mat
  })
  return(sample_cov)
}

#' Create missing data matrices function
#' @inheritParams create_data_list_meta
#' @returns list of missing data matrices
#' @keywords internal
prepare_missing_data_list <- function(data_list = NULL) {
  data_list$S <- fill_in_missing_covs(data_list$S)

  # Indicators of whether a variable is present by study
  data_list$valid_var <- do.call(cbind, lapply(data_list$S, function(s_mat) {
    apply(s_mat, 1, function(s) 0 + !(sum((s == 999)) == data_list$Ni - 1))
  }))

  # Count of missing covariance/correlation elements
  data_list$Nmiss <- sum(unlist(lapply(
    seq_len(length(data_list$S)), function(i) {
      idx <- which(data_list$valid_var[, i] == 1)
      s_mat <- data_list$S[[i]][idx, idx]
      sum(s_mat == 999) / 2
    }
  )))

  # Indicator of whether a valid variable has any missing covariance
  # elements by study
  data_list$miss_ind <- matrix(unlist(lapply(
    seq_len(length(data_list$S)), function(i) {
      s_mat <- data_list$S[[i]]
      idx <- which(data_list$valid_var[, i] == 1)
      res <- rep(0, data_list$Ni)
      res[idx] <- apply(s_mat[idx, idx], 2, function(s) any(s == 999) + 0)
      res
    }
  )), data_list$Ni)

  # Number of items that have missing covariance elements
  data_list$Nitem_miss <- sum(data_list$miss_ind)

  return(data_list)
}

#' Select random method and random case function
#'
#' @param meta IF TRUE, assume meta-analysis
#' @returns randomly selected method and case
#' @keywords internal
random_method_selection <- function(meta = FALSE) {
  accepted_methods <- method_hash()
  if (isTRUE(meta)) {
    accepted_methods <- accepted_methods[
      which(regexpr("WB", accepted_methods) < 0)
    ]
  }

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
    "none" = 100
  )
  converted_value <- converter_helper(search_term, list_methods)
  return(converted_value)
}

#' Type-meta hash function
#' @description A function that swaps type from string to integer
#' and vice-versa
#' @param elaborate (LOGICAL) If TRUE, print full names, otherwise
#' print abbreviations
#' @inheritParams converter_helper
#' @returns If search_term is integer, returns string and vice-versa
#' @keywords internal
type_hash <- function(search_term = NULL, elaborate = FALSE) {
  list_types <- c(
    "fe" = 1,
    "re" = 2
    # Add "dep" = 3 when ready
  )
  if (isTRUE(elaborate)) {
    # Only use when feeding integers
    list_types <- c(
      "Fixed-effects" = 1,
      "Random-effects" = 2
      # Add "Dependent-samples" = 3 when ready
    )
  }
  converted_value <- converter_helper(search_term, list_types)
  return(converted_value)
}

#' Multivariate normal density function
#'
#' @param x_mat the data matrix
#' @param mu mean vector
#' @param s_mat covariance matrix
#' @returns sum of casewise log-likelihood
#' @keywords internal
mb_ldmvn <- function(x_mat, mu, s_mat) {
  k <- nrow(x_mat)
  rooti <- backsolve(chol(s_mat), diag(k))
  quads <- colSums((crossprod(rooti, x_mat - mu))^2)
  return(-(k / 2) * log(2 * pi) + sum(log(diag(rooti))) - .5 * quads)
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
add_row_header <- function(
    kbl_object,
    table_to_print,
    search_term,
    extra = NULL) {
  if (any(table_to_print[1] == search_term)) {
    kbl_object <- kableExtra::pack_rows(
      kbl_object,
      paste0(search_term, " ", extra),
      which(table_to_print[1] == search_term)[1],
      rev(which(table_to_print[1] == search_term))[1]
    )
  }
  return(kbl_object)
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
  draws <- posterior::subset_draws(
    posterior::as_draws(stan_fit),
    variable = variable
  )

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
  if (data_list$meta == 1) {
    params[1] <- "rmsea_mn"
    from_list[1] <- "RMSEA"
  }

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
  duplicates <- duplicated(
    major_parameters[, c("mean", "median", "sd", "mad", "rhat", "ess_bulk")]
  )
  major_parameters <- major_parameters[
    !(duplicates & regexpr("phi\\_mat", major_parameters$variable) > 0),
  ]

  major_parameters[
    major_parameters$variable == "ppp", c("sd", "mad")
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

  idxs <- which(regexpr("Load\\_mat", major_parameters$variable) > 0)
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

  idxs <- which(regexpr("res\\_cor", major_parameters$variable) > 0)
  major_parameters <- modify_major_params(
    major_parameters, idxs,
    group = "Error correlations", op = "~~",
    from = indicator_labels[data_list$error_mat[, 1]],
    to = indicator_labels[data_list$error_mat[, 2]]
  )

  idxs <- which(regexpr("res\\_var", major_parameters$variable) > 0)
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

  idxs <- which(regexpr("phi\\_mat", major_parameters$variable) > 0)
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

  idxs <- which(regexpr("phi\\_cor", major_parameters$variable) > 0)
  major_parameters <- modify_major_params(
    major_parameters, idxs,
    group = "Inter-factor correlations", op = "~~",
    from = factor_labels[data_list$F_corr_mat[, 1]],
    to = factor_labels[data_list$F_corr_mat[, 2]]
  )

  idxs <- which(regexpr("r\\_square", major_parameters$variable) > 0)
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

  idxs <- which(regexpr("Coef\\_mat", major_parameters$variable) > 0)
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
    "Goodness of fit", "Latent regression coefficients", "R square",
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

#' Plotting param_type validation
#'
#' @param param_type param_type for plotting
#' @returns NULL
#' @keywords internal
validate_param_type <- function(param_type) {
  if (any(
    !param_type %in% c("all", "rm", "lo", "ev", "rc", "fc", "rsq", "co", "re")
  ) ||
    is.null(param_type)) {
    stop(paste0(
      "All param_type options must be in ",
      "c(\"all\", \"rm\", \"lo\", \"ev\", \"rc\", \"fc\", ",
      "\"rsq\", \"co\", \"re\")"
    ))
  }
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

#' Include minor factor residuals in model implied covariance matrix,
#' helper function
#'
#' @param omega_mat model implied covariance matrix to be updated
#' @param params vector of posterior samples from a single iteration
#' @param data_list Data list object passed to Stan
#' @returns A single model-impled covariance matrix
#' @keywords internal
include_residuals <- function(omega_mat, params, data_list) {
  if (data_list$method >= 90) {
    # there is no residual to include for this method
    return(omega_mat)
  }
  tv <- diag(omega_mat)
  n_re <- data_list$Ni * (data_list$Ni - 1) / 2
  re <- params[paste0("resids[", 1:n_re, "]")]
  rm <- params["rms_src_p[1]"]
  pos <- 0
  for (i in 2:data_list$Ni) {
    for (j in 1:(i - 1)) {
      pos <- pos + 1
      omega_mat[i, j] <- omega_mat[i, j] + re[pos] * rm * sqrt(tv[i] * tv[j])
      omega_mat[j, i] <- omega_mat[i, j]
    }
  }
  return(omega_mat)
}

#' Create model implied covariance matrix from CFA, helper function
#'
#' @param params vector of posterior samples from a single iteration
#' @param data_list Data list object passed to Stan
#' @param include_residuals (LOGICAL)
#' TRUE: Include minor factor residual covariances
#' in model-implied covariance matrix;
#' FALSE: Exclude them. If TRUE, different
#' models fit to the data will hardly be distinguishable.
#' @param all_lo loading indexes
#' @param all_ev error variance indexes
#' @param all_ph factor correlation indexes
#' @returns A single model-impled covariance matrix
#' @keywords internal
create_single_cfa_vcov_row <- function(
    params,
    data_list,
    include_residuals,
    all_lo,
    all_ev,
    all_ph) {
  lo_mat <- matrix(params[all_lo], nrow = data_list$Ni, ncol = data_list$Nf)
  ev <- params[all_ev]

  ph_mat <- diag(data_list$Nf)
  if (data_list$corr_fac == 1) {
    ph_mat <- matrix(params[all_ph], nrow = data_list$Nf, ncol = data_list$Nf)
  }

  lpl_mat <- lo_mat %*% ph_mat %*% t(lo_mat)

  lpe_mat <- matrix(0, nrow = data_list$Ni, ncol = data_list$Nce)
  if (ncol(lpe_mat) > 0) {
    rc <- params[paste0("res_cor[", 1:data_list$Nce, "]")]
    for (i in 1:data_list$Nce) {
      lpe_mat[data_list$error_mat[i, 1], i] <- sqrt(
        abs(rc[i]) * ev[data_list$error_mat[i, 1]]
      )
      lpe_mat[data_list$error_mat[i, 2], i] <- sign(rc[i]) * sqrt(
        abs(rc[i]) * ev[data_list$error_mat[i, 2]]
      )
    }
  }
  lpe_sq_mat <- tcrossprod(lpe_mat)

  d_ast <- ev - diag(lpe_sq_mat)

  omega_mat <- lpl_mat + lpe_sq_mat + diag(d_ast)

  if (include_residuals == TRUE) {
    omega_mat <- include_residuals(omega_mat, params, data_list)
  }

  return(omega_mat)
}

#' Create model implied covariance matrix from SEM, helper function
#'
#' @param params vector of posterior samples from a single iteration
#' @param data_list Data list object passed to Stan
#' @param include_residuals (LOGICAL)
#' TRUE: Include minor factor residual covariances
#' in model-implied covariance matrix; FALSE: Exclude them. If TRUE, different
#' models fit to the data will hardly be distinguishable.
#' @param all_lo loading indexes
#' @param all_ev error variance indexes
#' @param all_ph factor correlation indexes
#' @param all_co latent coefficient indexes
#' @param all_fv factor variance indexes
#' @returns A single model-impled covariance matrix
#' @keywords internal
create_single_sem_vcov_row <- function(
    params,
    data_list,
    include_residuals,
    all_lo,
    all_ev,
    all_ph,
    all_co,
    all_fv) {
  lo_mat <- matrix(params[all_lo], nrow = data_list$Ni, ncol = data_list$Nf)
  co_mat <- matrix(params[all_co], nrow = data_list$Nf, ncol = data_list$Nf)
  ev <- params[all_ev]
  fv <- params[all_fv]

  lombi_mat <- lo_mat %*% solve(diag(data_list$Nf) - co_mat)

  fpe_mat <- matrix(0, nrow = data_list$Nf, ncol = data_list$Nf_corr)
  if (ncol(fpe_mat) > 0) {
    fc <- params[all_ph]
    for (i in 1:data_list$Nf_corr) {
      fpe_mat[data_list$F_corr_mat[i, 1], i] <- sqrt(
        abs(fc[i]) * ev[data_list$F_corr_mat[i, 1]]
      )
      fpe_mat[data_list$F_corr_mat[i, 2], i] <- sign(fc[i]) * sqrt(
        abs(fc[i]) * ev[data_list$F_corr_mat[i, 2]]
      )
    }
  }
  fpe_sq_mat <- tcrossprod(fpe_mat)
  fvr <- fv - diag(fpe_sq_mat)

  lpe_mat <- matrix(0, nrow = data_list$Ni, ncol = data_list$Nce)
  if (ncol(lpe_mat) > 0) {
    rc <- params[paste0("res_cor[", 1:data_list$Nce, "]")]
    for (i in 1:data_list$Nce) {
      lpe_mat[data_list$error_mat[i, 1], i] <- sqrt(
        abs(rc[i]) * ev[data_list$error_mat[i, 1]]
      )
      lpe_mat[data_list$error_mat[i, 2], i] <- sign(rc[i]) * sqrt(
        abs(rc[i]) * ev[data_list$error_mat[i, 2]]
      )
    }
  }
  lpe_sq_mat <- tcrossprod(lpe_mat)

  d_ast <- ev - diag(lpe_sq_mat)

  omega_mat <- lombi_mat %*% (
    fpe_sq_mat + diag(fvr)
  ) %*% t(lombi_mat) + lpe_sq_mat + diag(d_ast)

  if (include_residuals == TRUE) {
    omega_mat <- include_residuals(omega_mat, params, data_list)
  }

  return(omega_mat)
}

#' Create model implied cov matrix or log-likelihood helper function
#'
#' @param mat Matrix of posterior samples
#' @param data_list Data list object passed to Stan
#' @param include_residuals (LOGICAL)
#' TRUE: Include minor factor residual covariances
#' in model-implied covariance matrix; FALSE: Exclude them. If TRUE, different
#' models fit to the data will hardly be distinguishable.
#' @param return_ll (LOGICAL)
#' TRUE: Return matrix of log-likelihood
#' FALSE: Return matrix of covariances
#' @returns Returns a matrix, which one depends on
#' \code{return_ll} argument.
#' @keywords internal
create_mi_vcov_ll <- function(
    mat,
    data_list,
    include_residuals,
    return_ll = FALSE) {
  all_ev <- paste0("res_var[", 1:data_list$Ni, "]")
  all_lo <- paste0("Load_mat[", apply(which(
    data_list$loading_pattern != 2,
    arr.ind = TRUE
  ), 1, paste0, collapse = ","), "]")
  all_ph <- NULL

  returned_mat <- matrix()
  if (data_list$method == 91) {
    sigma_cols <- colnames(mat)
    sigma_cols <- sigma_cols[which(regexpr("Sigma\\[", sigma_cols) > 0)]
    returned_mat <- mat[, sigma_cols]
    if (isTRUE(return_ll)) {
      mu <- rep(0, data_list$Ni)
      y_dat_t <- t(data_list$Y) - colMeans(data_list$Y)
      returned_mat <- apply(returned_mat, 1, function(m) {
        m_vcov <- matrix(m, nrow = data_list$Ni, ncol = data_list$Ni)
        mb_ldmvn(y_dat_t, mu, m_vcov)
      })
    } else {
      returned_mat <- t(returned_mat)
    }
  } else if (data_list$sem_indicator == 0) {
    if (data_list$corr_fac == 1) {
      all_ph <- paste0("phi_mat[", apply(which(
        diag(data_list$Nf) != 2,
        arr.ind = TRUE
      ), 1, paste0, collapse = ","), "]")
    }

    if (isFALSE(return_ll)) {
      returned_mat <- apply(
        mat, 1, create_single_cfa_vcov_row,
        data_list = data_list, include_residuals = include_residuals,
        all_lo = all_lo, all_ev = all_ev, all_ph = all_ph
      )
    } else if (isTRUE(return_ll)) {
      mu <- rep(0, data_list$Ni)
      y_dat_t <- t(data_list$Y) - colMeans(data_list$Y)
      returned_mat <- apply(mat, 1, function(m) {
        m_vcov <- create_single_cfa_vcov_row(
          m,
          data_list = data_list, include_residuals = include_residuals,
          all_lo = all_lo, all_ev = all_ev, all_ph = all_ph
        )
        mb_ldmvn(y_dat_t, mu, m_vcov)
      })
    }
  } else if (data_list$sem_indicator == 1) {
    # Use _u Coefs and Loads as these are unstandardized
    all_lo <- paste0("Load_mat_u[", apply(which(
      data_list$loading_pattern != 2,
      arr.ind = TRUE
    ), 1, paste0, collapse = ","), "]")
    all_co <- paste0("Coef_mat_u[", apply(which(
      data_list$coef_pattern != 2,
      arr.ind = TRUE
    ), 1, paste0, collapse = ","), "]")
    all_fv <- paste0("phi_var[", 1:data_list$Nf, "]")
    if (data_list$Nf_corr > 0) {
      all_ph <- paste0("phi_cor[", 1:data_list$Nf_corr, "]")
    }

    if (isFALSE(return_ll)) {
      returned_mat <- apply(
        mat, 1, create_single_sem_vcov_row,
        data_list = data_list, include_residuals = include_residuals,
        all_lo = all_lo, all_ev = all_ev, all_ph = all_ph,
        all_co = all_co, all_fv = all_fv
      )
    } else if (isTRUE(return_ll)) {
      mu <- rep(0, data_list$Ni)
      y_dat_t <- t(data_list$Y) - colMeans(data_list$Y)
      returned_mat <- apply(mat, 1, function(m) {
        m_vcov <- create_single_sem_vcov_row(
          m,
          data_list = data_list, include_residuals = include_residuals,
          all_lo = all_lo, all_ev = all_ev, all_ph = all_ph,
          all_co = all_co, all_fv = all_fv
        )
        mb_ldmvn(y_dat_t, mu, m_vcov)
      })
    }
  }

  return(returned_mat)
}
