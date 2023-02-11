#' Visualize standardized residual covariances
#'
#' @description Visualize distribution of standardized residual covariances
#' @param object (mbsem) An object of \code{\link{mbsem-class}}
#' returned by \code{\link{minorbsem}}.
#' @param type (string) Either: "range" for lineranges (by default)
#' or "matrix" for a matrix with point estimates in the lower half.
#' @returns ggplot object
#' @examples
#' fit <- minorbsem("F1 =~ x1 + x2 + x3
#'                   F2 =~ x4 + x5 + x6
#'                   F3 =~ x7 + x8 + x9", HS)
#' plot_residuals(fit)
#' @export
plot_residuals <- function(object, type = "matrix") {
  parameter <- `50%` <- `5%` <- `95%` <- lo <- hi <- item_2 <- item_1 <- NULL

  if (!type %in% c("range", "matrix")) {
    stop("type must be either \"range\" or \"matrix\"")
  }

  clean_post_df <- prepare_stan_plot_data(object)
  clean_post_df <- clean_post_df[clean_post_df$param_class == "re", ]
  clean_post_df$parameter <- gsub("re: ", "", clean_post_df$parameter)
  plot_df <- stats::aggregate(
    value ~ parameter, clean_post_df,
    FUN = function(xs) {
      stats::quantile(xs, c(.5, .25, .75, .05, .95, .025, .975))
    }
  )
  plot_df <- cbind(plot_df, plot_df$value)
  plot_df$lo <- plot_df$`2.5%`
  plot_df$hi <- plot_df$`97.5%`

  p <- list()
  if (type == "range") {
    horiz_lines <- c(-.05, 0, .05)
    if (min(plot_df$`2.5%`) < -.1) {
      horiz_lines <- c(horiz_lines, -.1)
    }
    if (max(plot_df$`97.5%`) > .1) {
      horiz_lines <- c(horiz_lines, .1)
    }
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(parameter, `50%`)) +
      ggplot2::geom_point(shape = 4) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = `5%`, ymax = `95%`),
        size = 3, alpha = .25
      ) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = lo, ymax = hi),
        size = .25
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = .5)
      ) +
      ggplot2::geom_hline(
        yintercept = horiz_lines, linetype = 2, alpha = .5
      ) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(linewidth = .5)
      ) +
      ggplot2::labs(
        x = "Indicator-pairs",
        y = "90% & 95% credible intervals"
      )
  } else if (type == "matrix") {
    # Create items, use factors to ensure correct ordering
    plot_df$item_1 <- unlist(lapply(strsplit(plot_df$parameter, "~~"), "[[", 1))
    item_list_1 <- unique(unlist(lapply(strsplit(
      clean_post_df$parameter, "~~"
    ), "[[", 1)))
    plot_df$item_1 <- factor(plot_df$item_1, rev(item_list_1))
    plot_df$item_2 <- unlist(lapply(strsplit(plot_df$parameter, "~~"), "[[", 2))
    item_list_2 <- unique(unlist(lapply(strsplit(
      clean_post_df$parameter, "~~"
    ), "[[", 2)))
    plot_df$item_2 <- factor(plot_df$item_2, item_list_2)
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(item_2, item_1)) +
      ggplot2::geom_tile(
        ggplot2::aes(alpha = abs(`50%`)),
        fill = "#444444", col = 1
      ) +
      ggplot2::geom_text(ggplot2::aes(
        label = scales::number(`50%`, .01),
        alpha = abs(`50%`)
      )) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = .5)
      ) +
      ggplot2::guides(alpha = "none") +
      ggplot2::theme(
        axis.title = ggplot2::element_blank()
      )
  }

  return(p)
}
