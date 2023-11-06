#' Visualize standardized residual covariances
#'
#' @description Visualize distribution of standardized residual covariances
#' @param object (mbsem) An object of \code{\link{mbsem-class}}
#' returned by \code{\link{minorbsem}}.
#' @param type (string) Either: "range" for lineranges (by default)
#' or "matrix" for a matrix with point estimates in the lower half.
#' @returns ggplot object
#' @examples
#' \dontrun{
#' fit <- minorbsem("F1 =~ x1 + x2 + x3
#'                   F2 =~ x4 + x5 + x6
#'                   F3 =~ x7 + x8 + x9", HS)
#' plot_residuals(fit)
#' }
#' @export
plot_residuals <- function(object, type = "matrix") {
  parameter <- `50%` <- `5%` <- `95%` <- lo <- hi <- item_2_f <- item_1_f <-
    NULL

  if (!type %in% c("range", "matrix")) {
    stop("type must be either \"range\" or \"matrix\"")
  }

  if (object@data_list$method >= 90) {
    stop(paste0(
      "There are no residuals to plot when ",
      "method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
    ))
  }

  ind_names <- rownames(object@data_list$loading_pattern)
  plot_df <- posterior::summarise_draws(
    object@stan_fit$draws("Resid"),
    ~ stats::quantile(.x, c(.5, .25, .75, .05, .95, .025, .975), na.rm = TRUE)
  )
  plot_df$item_1 <- as.integer(gsub(
    "Resid\\[|,\\d+\\]", "", plot_df$variable
  ))
  plot_df$item_2 <- as.integer(gsub(
    "Resid\\[\\d+,|\\]", "", plot_df$variable
  ))
  plot_df <- plot_df[
    plot_df$item_1 < plot_df$item_2, ,
    drop = FALSE
  ]
  plot_df$item_1_f <- factor(
    plot_df$item_1, seq_along(ind_names), ind_names
  )
  plot_df$item_2_f <- factor(
    plot_df$item_2, rev(seq_along(ind_names)), rev(ind_names)
  )
  plot_df$parameter <- paste0(
    as.character(plot_df$item_1_f), "~~", as.character(plot_df$item_2_f)
  )

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
        linewidth = 3, alpha = .25
      ) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = lo, ymax = hi),
        linewidth = .25
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
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(item_1_f, item_2_f)) +
      ggplot2::geom_tile(
        ggplot2::aes(alpha = abs(`50%`)),
        fill = "#444444", col = 1
      ) +
      ggplot2::geom_text(ggplot2::aes(
        label = gsub(
          "0\\.", "\\.",
          trimws(format(
            round(`50%`, 3),
            scientific = FALSE, digits = 2, nsmall = 2
          ))
        ),
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
