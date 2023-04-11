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
  n_obs <- paste0(object@data_list$Np, collapse = ", ")

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

  type_str <- ""
  if (object@data_list$meta == 1) {
    type_str <- paste0(
      type_hash(object@data_list$type, elaborate = TRUE),
      ", "
    )
  }

  table_to_print[
    table_to_print$group == "Latent regression coefficients", "group"
  ] <- "Latent regression coefficients (outcome ~ predictor)"

  caption_str <- paste0(
    "Parameter estimates (",
    type_str,
    "method = ", method_str,
    ", sample size(s) = ", n_obs, ")"
  )

  result <- huxtable::huxtable(table_to_print)
  result <- huxtable::column_to_header(result, "group")
  result <- huxtable::set_caption(result, caption_str)

  gof_row <- which(result$from == "Goodness of fit")
  disp_row <- which(result$from == "Dispersion between and within clusters")
  gof_disp_rows <- c(gof_row + 1:2, disp_row + 1:3)

  for (row_id in gof_disp_rows) {
    result <- huxtable::merge_cells(result, row_id, 1:3)
  }

  footnote_str <- ""

  if (object@data_list$sem_indicator == 1) {
    footnote_str <- paste0(
      "Loadings are fully standardized, error variances are not shown."
    )
  }

  if (isTRUE(simple) && object@data_list$method == 4) {
    footnote_str <- paste0(
      footnote_str, "\n",
      "Mean and SD RMSE are median and mad respectively ",
      "because RMSE is usually extremely right-skewed."
    )
  }

  result <- huxtable::add_footnote(result, text = footnote_str)

  if (isTRUE(simple)) {
    huxtable::number_format(result)[, 4:8] <- list(
      function(x) trimws(format(round(x, digits), nsmall = digits))
    )
    result <- huxtable::set_number_format(result, col = 9, value = 0)
  } else if (isFALSE(simple)) {
    huxtable::number_format(result)[, 4:10] <- list(
      function(x) trimws(format(round(x, digits), nsmall = digits))
    )
    result <- huxtable::set_number_format(result, col = 11:12, value = 0)
    result <- huxtable::insert_row(
      result,
      "Relation", rep("", 2),
      "Location", rep("", 1),
      "Dispersion", rep("", 3),
      "Parameter convergence", rep("", 2),
      after = 0
    )
    result <- huxtable::merge_cells(result, 1, 1:3)
    result <- huxtable::merge_cells(result, 1, 4:5)
    result <- huxtable::merge_cells(result, 1, 6:9)
    result <- huxtable::merge_cells(result, 1, 10:12)
    result <- huxtable::set_align(result, 1, huxtable::everywhere, "center")
    result <- huxtable::set_valign(result, 1, huxtable::everywhere, "bottom")
    result <- huxtable::set_tb_padding(result, 1, huxtable::everywhere, 10)
    result <- huxtable::set_header_rows(result, 1, TRUE)
    result <- huxtable::set_bottom_border(
      result, 1,
      col = huxtable::everywhere
    )
    result <- huxtable::set_right_border(
      result,
      row = huxtable::everywhere, col = c(3, 5, 9)
    )
  }

  header_rows <- rev(which(
    rownames(result) == "" | substr(rownames(result), 1, 1) == "."
  ))[-1]
  header_rows <- c(header_rows, header_rows - 1)
  result <- huxtable::set_bottom_border(
    result, header_rows, huxtable::everywhere
  )

  result <- huxtable::style_headers(result, bold = TRUE)

  if (!is.null(save_html) && is.character(save_html)) {
    huxtable::quick_html(result, file = save_html)
  }

  print(result)

  return()
}
