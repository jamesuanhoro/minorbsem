#' Pretty print model results
#'
#' @description Nice printing of model results, optionally produces HTML document
#' @param clean_results A model fitted with minorbsem.
#' @param digits (positive integer) Number of decimal places to print in table
#' @param simple (Logical) TRUE to produce table with less information about parameters;
#' FALSE: produces table with more information
#' @param save_html (string) Optional file name to save table as HTML
#' @returns NULL
#' @examples
#' fit <- minorbsem("F1 =~ x1 + x2 + x3
#'                   F2 =~ x4 + x5 + x6
#'                   F3 =~ x7 + x8 + x9", HS)
#' pretty_print_summary(fit)
#' pretty_print_summary(fit, simple = FALSE)
#' @export
pretty_print_summary <- function(
    clean_results, digits = 3, simple = TRUE,
    save_html = NULL) {
  table_to_print <- clean_results$major_parameters

  if (simple) {
    table_to_print <- table_to_print[
      c("group", "from", "op", "to", "mean", "sd", "q5", "q95", "rhat", "ess_bulk")
    ]
  }

  result <- kableExtra::kbl(
    table_to_print[, -1],
    caption = "Parameter estimates",
    digits = digits
  )
  result <- kableExtra::kable_paper(result)
  result <- kableExtra::kable_styling(result, full_width = FALSE)
  result <- kableExtra::pack_rows(result, "RMSE(residuals)", 1, 1)

  result <- add_row_header(
    result, table_to_print, "Latent regression coefficients", "(outcome ~ predictor)"
  )
  result <- add_row_header(
    result, table_to_print, "Factor loadings", "(factor =~ indicator)"
  )
  result <- add_row_header(
    result, table_to_print, "Inter-factor correlations", "(factor_x ~~ factor_y)"
  )
  result <- add_row_header(
    result, table_to_print, "Factor variances", "(factor_x ~~ factor_x)"
  )
  result <- add_row_header(
    result, table_to_print, "Residual variances", "(indicator_x ~~ indicator_x)"
  )
  result <- add_row_header(
    result, table_to_print, "Error correlations", "(indicator_x ~~ indicator_y)"
  )

  if (!simple) {
    result <- kableExtra::add_header_above(
      result, c(
        "Relation" = 3, "Location" = 2, "Dispersion" = 4,
        "Parameter convergence" = 3
      )
    )
  }

  if (!is.null(save_html) & is.character(save_html)) {
    kableExtra::save_kable(result, file = save_html)
  }

  print(result)
}
