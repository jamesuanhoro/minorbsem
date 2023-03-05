#' Pretty print model results
#'
#' @description Nice printing of model results,
#' optionally produces HTML document
#' @param object (mbsem) An object of \code{\link{mbsem-class}}
#' returned by \code{\link{minorbsem}}.
#' @param interval (real in (0, 1)) Credible interval to return.
#' @param digits (positive integer) Number of decimal places to print in table
#' @param simple (Logical) TRUE to produce table with less information
#' about parameters;
#' FALSE: produces table with more information
#' @param save_html (string) Optional file name to save table as HTML
#' @returns NULL
#' @examples
#' \dontrun{
#' fit <- minorbsem("F1 =~ x1 + x2 + x3
#'                   F2 =~ x4 + x5 + x6
#'                   F3 =~ x7 + x8 + x9", HS)
#' pretty_print_summary(fit)
#' pretty_print_summary(fit, simple = FALSE)
#' }
#' @export
pretty_print_summary <- function(
    object, interval = .9, digits = 3, simple = TRUE,
    save_html = NULL) {
  stopifnot(inherits(object, "mbsem"))

  if (interval <= 0 || interval >= 1 || !is.numeric(interval)) {
    stop("Interval must be a number between 0 and 1")
  }

  table_to_print <- create_major_params(
    stan_fit = object@stan_fit,
    data_list = object@data_list,
    interval = interval
  )
  method_str <- method_hash(object@data_list$method)
  n_obs <- object@data_list$Np

  if (isTRUE(simple)) {
    if (object@data_list$method == 4) {
      table_to_print[2, "mean"] <- table_to_print[2, "median"]
      table_to_print[2, "sd"] <- table_to_print[2, "mad"]
    }
    # Use index because of quantile name changes
    table_to_print <- table_to_print[, c(1:5, 7, 9:10, 11:12)]
  }

  if (object@data_list$method == 100) {
    # Remove RMSE index
    table_to_print[2, 5:ncol(table_to_print)] <- NA_real_
  }

  if (object@data_list$method %in% 90:99) {
    # Remove RMSE index
    table_to_print[2, 2] <- "RMSEA"
  }

  result <- kableExtra::kbl(
    table_to_print[, -1],
    row.names = FALSE,
    booktabs = TRUE,
    caption = paste0(
      "Parameter estimates (method = ", method_str,
      ", sample size = ", n_obs, ")"
    ),
    digits = digits
  )

  if (simple && object@data_list$method == 4) {
    result <- kableExtra::footnote(
      result,
      general = paste0(
        "Mean and SD RMSE are median and mad respectively ",
        "because RMSE is usually extremely right-skewed."
      )
    )
  }

  result <- kableExtra::kable_paper(result)
  result <- kableExtra::kable_styling(result, full_width = FALSE)

  result <- add_row_header(
    result, table_to_print,
    "Goodness of fit", ""
  )
  result <- add_row_header(
    result, table_to_print,
    "Latent regression coefficients", "(outcome ~ predictor)"
  )
  result <- add_row_header(
    result, table_to_print,
    "Factor loadings", "(factor =~ indicator)"
  )
  result <- add_row_header(
    result, table_to_print,
    "Inter-factor correlations", "(factor_x ~~ factor_y)"
  )
  result <- add_row_header(
    result, table_to_print,
    "R square", "(factor_x ~~ factor_x)"
  )
  result <- add_row_header(
    result, table_to_print,
    "Residual variances", "(indicator_x ~~ indicator_x)"
  )
  result <- add_row_header(
    result, table_to_print,
    "Error correlations", "(indicator_x ~~ indicator_y)"
  )

  if (isFALSE(simple)) {
    result <- kableExtra::add_header_above(
      result, c(
        "Relation" = 3, "Location" = 2, "Dispersion" = 4,
        "Parameter convergence" = 3
      )
    )
  }

  if (!is.null(save_html) && is.character(save_html)) {
    kableExtra::save_kable(result, file = save_html)
  }

  print(result)

  return(result)
}
