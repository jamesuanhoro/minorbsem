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
  ) |>
    kableExtra::kable_paper() |>
    kableExtra::kable_styling(full_width = FALSE) |>
    kableExtra::pack_rows(
      "RMSE(residuals)",
      1, 1
    )

  result <- add_row_header(result, table_to_print, "Latent regression coefficients")
  result <- add_row_header(result, table_to_print, "Factor loadings")
  result <- add_row_header(result, table_to_print, "Inter-factor correlations")
  result <- add_row_header(result, table_to_print, "Factor variances")
  result <- add_row_header(result, table_to_print, "Residual variances")
  result <- add_row_header(result, table_to_print, "Error correlations")

  if (!simple) {
    result <- result |>
      kableExtra::add_header_above(c(
        "Relation" = 3, "Location" = 2, "Dispersion" = 4,
        "Parameter convergence" = 3
      ))
  }

  if (!is.null(save_html) & is.character(save_html)) {
    kableExtra::save_kable(result, file = save_html)
  }

  print(result)
}
