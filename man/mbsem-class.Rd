\name{mbsem-class}
\Rdversion{1.1}
\docType{class}
\alias{mbsem-class}
\alias{plot,mbsem-method}
\alias{logLik,mbsem-method}
\alias{residuals,mbsem-method}
\alias{show,mbsem-method}

\title{Class \code{"mbsem"}: For representing a model fitted with minorbsem}
\description{
  The \code{mbsem} class represents a fitted model.
  It contains summary statistics for structural parameters,
  summary statistics for standardized residual covariances,
  the data used to fit the model, and the fitted Stan model.
}
\section{Objects from the Class}{
  Objects can be created by calls to the \code{\link{minorbsem}} function.
}
\section{Slots}{
  \describe{
    \item{\code{major_parameters}:}{Summary statistics for structural parameters}
    \item{\code{minor_factor_matrix}:}{Summary statistics for standardized residual covariances}
    \item{\code{data_list}:}{Data used to fit the model}
    \item{\code{stan_fit}:}{Fitted CmdStan model}
    \item{\code{version}:}{Package version used to fit model}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{
      \code{signature(x = "mbsem", type = "hist", param_type = c("rm", "co", "lo", "fc", "fv"))}:
      Visualizes the posterior distribution of requested parameter types in
      \code{param_type = } argument.
      If \code{type = "hist"}, returns histogram(s).
      If \code{type = "trace"}, returns traceplot(s).
      See \code{\link{parameter_hist}} or \code{\link{parameter_trace}} for
      elaboration of \code{param_type} argument.
    }
    \item{logLik}{
      \code{signature(object = "mbsem", include_residuals = FALSE)}:
      Returns the casewise log-likelihood as long as the full data is available.
      \code{include_residuals = FALSE} ensures the minor factor residual
      covariances are not included when computing the log-likelihood.
      See \code{\link{casewise_log_likelihood}} for more details.
    }
    \item{residuals}{
      \code{signature(object = "mbsem", method = "matrix")}:
      If \code{method = "table"}, returns summary statistics for standardized residual covariances.
      If \code{method = "matrix"}, returns a plot of the matrix of standardized residual covariances.
      If \code{method = "range"}, returns a linerange plot of standardized residual covariances.
    }
    \item{show}{
      \code{signature(object = "mbsem")}:
      Pretty printing of model results.
      See \code{\link{pretty_print_summary}} for more printing options.
    }
  }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{parameter_hist}}, \code{\link{parameter_trace}},
\code{\link{casewise_log_likelihood}}, \code{\link{plot_residuals}}
\code{\link{pretty_print_summary}}
}
\examples{
  \dontrun{
    fit_1 <- minorbsem("F1 =~ x1 + x2 + x3
                        F2 =~ x4 + x5 + x6
                        F3 =~ x7 + x8 + x9", HS)
    plot(fit_1)
  }
}
\keyword{classes}
