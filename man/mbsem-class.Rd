\name{mbsem-class}
\Rdversion{1.1}
\docType{class}
\alias{mbsem-class}
\alias{fitted,mbsem-method}
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
    \item{\code{priors}:}{Priors used to fit the model}
    \item{\code{stan_fit}:}{Fitted CmdStan model}
    \item{\code{version}:}{Package version used to fit model}
  }
}
\section{Methods}{
  \describe{
    \item{logLik}{
      \code{signature(object = "mbsem")}:
      Returns the casewise log-likelihood as long as the full data is available.
    }
    \item{fitted}{
      \code{signature(object = "mbsem")}:
      Returns the posterior distribution of the model-implied covariance matrix
      as a #(iterations) by #(items ^ 2) matrix.
    }
    \item{residuals}{
      \code{signature(object = "mbsem", standardized = TRUE)}:
      Returns the posterior distribution of residual covariances reflecting the influences of
      minor factors as a #(iterations) by #(items ^ 2) matrix.
      If \code{standardized = TRUE}, returns standardized residual covariances.
      If \code{standardized = FALSE}, returns UNstandardized residual covariances.
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
\code{\link{plot_residuals}},
\code{\link{pretty_print_summary}}
}
\examples{
  \dontrun{
    fit_1 <- minorbsem("F1 =~ x1 + x2 + x3
                        F2 =~ x4 + x5 + x6
                        F3 =~ x7 + x8 + x9", HS)
    # Model-implied covariances
    mod_imp_cov <- fitted(fit_1)
    # Get average of model-implied covariance matrix
    matrix(colMeans(mod_imp_cov), nrow = fit_1@data_list$Ni)
  }
}
\keyword{classes}
