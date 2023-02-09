% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_trace.R
\name{parameter_trace}
\alias{parameter_trace}
\title{Traceplots of parameter posterior distribution}
\usage{
parameter_trace(clean_results, param_type = c("rm", "co", "lo", "fc", "fv"))
}
\arguments{
\item{clean_results}{A model fitted with minorbsem}

\item{param_type}{(string vector) Choose from a list of options:
"all" = all structural parameters;
"rm" = Root Mean square error of standardized residual covariances;
"lo" = LOadings;
"ev" = Error Variances;
"rc" = Residual Correlations;
"fc" = Factor Correlations;
"fv" = Factor Variances;
"co" = latent regression COefficients;
"re" = standardized REsidual covariances}
}
\value{
ggplot object
}
\description{
Produce traceplots of parameter posterior distribution,
option to limit plot by type of parameter.
}
\examples{
fit <- minorbsem("F1 =~ x1 + x2 + x3
                  F2 =~ x4 + x5 + x6
                  F3 =~ x7 + x8 + x9", HS)
parameter_trace(fit)
parameter_trace(fit, param_type = "all")
parameter_trace(fit, param_type = c("rm", "lo", "fc"))
}