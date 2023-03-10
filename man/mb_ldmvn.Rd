% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{mb_ldmvn}
\alias{mb_ldmvn}
\title{Multivariate normal density function}
\usage{
mb_ldmvn(x_mat, mu, s_mat)
}
\arguments{
\item{x_mat}{the data matrix}

\item{mu}{mean vector}

\item{s_mat}{covariance matrix}
}
\value{
sum of casewise log-likelihood
}
\description{
Multivariate normal density function
}
\keyword{internal}
