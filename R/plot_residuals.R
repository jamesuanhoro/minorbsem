#' Visualize standardized residual covariances
#'
#' @description Visualize distribution of standardized residual covariances clearly
#' @param clean_results A model fitted with minorbsem
#' @returns ggplot object
#' @examples
#' fit <- minorbsem("F1 =~ x1 + x2 + x3
#'                   F2 =~ x4 + x5 + x6
#'                   F3 =~ x7 + x8 + x9", HS)
#' plot_residuals(fit)
#' @export
plot_residuals <- function(clean_results) {
  parameter <- `50%` <- `5%` <- `95%` <- `2.5%` <- `97.5%` <- NULL

  clean_post_df <- prepare_stan_plot_data(clean_results$stan_fit, clean_results$data_list)
  clean_post_df <- clean_post_df[clean_post_df$param_class == "re", ]
  clean_post_df$parameter <- gsub("re: ", "", clean_post_df$parameter)
  plot_df <- stats::aggregate(value ~ parameter, clean_post_df, FUN = function(xs) {
    stats::quantile(xs, c(.5, .25, .75, .05, .95, .025, .975))
  })
  plot_df <- cbind(plot_df, plot_df$value)

  horiz_lines <- c(-.05, 0, .05)
  if (min(plot_df$`2.5%`) < -.1) {
    horiz_lines <- c(horiz_lines, -.1)
  }
  if (max(plot_df$`97.5%`) > .1) {
    horiz_lines <- c(horiz_lines, .1)
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(parameter, `50%`)) +
    ggplot2::geom_point(shape = 4) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = `5%`, ymax = `95%`), size = 3, alpha = .25) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = `2.5%`, ymax = `97.5%`), size = .25) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5)) +
    ggplot2::geom_hline(yintercept = horiz_lines, linetype = 2, alpha = .5) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Indicator-pairs",
      y = "90% & 95% credible intervals"
    )

  return(p)
}
