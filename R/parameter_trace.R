#' Traceplots of parameter posterior distribution
#'
#' @description Produce traceplots of parameter posterior distribution,
#' option to limit plot by type of parameter.
#' @param object (mbsem) An object of \code{\link{mbsem-class}}
#' returned by \code{\link{minorbsem}}.
#' @param param_type (string vector) Choose from a list of options:
#' "all" = all structural parameters;
#' "rm" = Root Mean square error of standardized residual covariances;
#' "lo" = LOadings;
#' "ev" = Error Variances;
#' "rc" = Residual Correlations;
#' "fc" = Factor Correlations;
#' "fv" = Factor Variances;
#' "co" = latent regression COefficients;
#' "re" = standardized REsidual covariances
#' @returns ggplot object
#' @examples
#' fit <- minorbsem("F1 =~ x1 + x2 + x3
#'                   F2 =~ x4 + x5 + x6
#'                   F3 =~ x7 + x8 + x9", HS)
#' parameter_trace(fit)
#' parameter_trace(fit, param_type = "all")
#' parameter_trace(fit, param_type = c("rm", "lo", "fc"))
#' @export
parameter_trace <- function(
    object,
    param_type = c("rm", "co", "lo", "fc", "fv")) {
  .iteration <- value <- .chain <- NULL

  # param_type must in options
  validate_param_type(param_type)

  clean_post_df <- prepare_stan_plot_data(object)

  p <- list()
  if (any(param_type == "all")) {
    plot_df <- clean_post_df[
      clean_post_df$param_class != "re",
    ]
    p <- ggplot2::ggplot(
      plot_df, ggplot2::aes(.iteration, value, col = factor(.chain))
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_classic() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        legend.position = "top"
      ) +
      ggplot2::facet_wrap(~parameter, scales = "free_y") +
      ggplot2::labs(
        x = "Iteration",
        y = "Parameter estimates",
        col = "Chain"
      )
  } else {
    plot_df <- clean_post_df[
      clean_post_df$param_class %in% param_type,
    ]
    if (nrow(plot_df) == 0) {
      stop(paste0(
        "Selected param_type option(s) is/are not in the fitted model"
      ))
    }
    p <- ggplot2::ggplot(
      plot_df, ggplot2::aes(.iteration, value, col = factor(.chain))
    ) +
      ggplot2::geom_line() +
      ggplot2::theme_classic() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        legend.position = "top"
      ) +
      ggplot2::facet_wrap(~parameter, scales = "free_y") +
      ggplot2::labs(
        x = "Iteration",
        y = "Parameter estimates",
        col = "Chain"
      )
  }

  return(p)
}
