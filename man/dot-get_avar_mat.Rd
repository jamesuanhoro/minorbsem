% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.get_avar_mat}
\alias{.get_avar_mat}
\title{Asymptotic variance matrix of lower half vector of log(correlation matrix)}
\usage{
.get_avar_mat(r_mat, n, acov_mat = NULL)
}
\arguments{
\item{r_mat}{(matrix) Correlation matrix}

\item{n}{(positive integer) Sample size}

\item{acov_mat}{(matrix) Asymptotic variance matrix of lower half vector
of correlation matrix.}
}
\value{
Asymptotic variance matrix of strict lower half vector
of log(correlation matrix).
}
\description{
Asymptotic variance matrix of lower half vector of log(correlation matrix)
}
\keyword{internal}
