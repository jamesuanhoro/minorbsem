% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minorbpa.R
\name{minorbpa}
\alias{minorbpa}
\title{Fit Bayesian path analysis models}
\usage{
minorbpa(
  model = NULL,
  data = NULL,
  sample_cov = NULL,
  sample_nobs = NULL,
  data_list = NULL,
  method = "normal",
  orthogonal = FALSE,
  correlation = FALSE,
  centered = TRUE,
  seed = 12345,
  warmup = 1000,
  sampling = 1000,
  refresh = (warmup + sampling)/10,
  adapt_delta = 0.9,
  max_treedepth = 10,
  chains = 3,
  ncores = max(parallel::detectCores() - 2, 1),
  priors = new_mbsempriors(),
  show = TRUE,
  show_messages = TRUE,
  compute_ll = FALSE,
  acov_mat = NULL,
  ret_data_list = FALSE
)
}
\arguments{
\item{model}{A description of the user-specified model, lavaan syntax.}

\item{data}{An optional data frame containing the observed variables used in
the model.}

\item{sample_cov}{(matrix) sample variance-covariance matrix.
The rownames and/or colnames must contain the observed variable names.}

\item{sample_nobs}{(positive integer) Number of observations if the full
data frame
is missing and only sample covariance matrix is given.}

\item{data_list}{(list) A modified version of the data_list returned
by minorbsem. Can be used to modify specific priors, see example below.}

\item{method}{(character) One of "normal", "lasso", "logistic",
"GDP", "WB", "WB-cond", "WW", or "none". See details below.}

\item{orthogonal}{(logical) constrain factors orthogonal, must be TRUE to fit
bifactor models.}

\item{correlation}{(LOGICAL)
If TRUE: perform correlation structure analysis based on logarithm
of a matrix transformation
\insertCite{archakov_new_2021}{minorbsem};
If FALSE (default): perform covariance structure analysis.}

\item{centered}{(LOGICAL)
Only relevant for WB-cond and WW methods when \code{correlation = TRUE}.
If TRUE (default): Use a centered parameterization;
If FALSE: Use a non-centered parameterization.}

\item{seed}{(positive integer) seed, set to obtain replicable results.}

\item{warmup}{(positive integer) The number of warmup iterations to run per
chain.}

\item{sampling}{(positive integer) The number of post-warmup iterations to
run per chain, retained for inference.}

\item{refresh}{(positive integer) How often to print the status
of the sampler.}

\item{adapt_delta}{(real in (0, 1)) Increase to resolve divergent
transitions.}

\item{max_treedepth}{(positive integer) Increase to resolve problems with
maximum tree depth.}

\item{chains}{(positive integer) The number of Markov chains to run.}

\item{ncores}{(positive integer) The number of chains to run in parallel.}

\item{priors}{An object of \code{\link{mbsempriors-class}}.
See \code{\link{new_mbsempriors}} for more information.}

\item{show}{(Logical) If TRUE, show table of results, if FALSE, do not
show table of results. As an example, use FALSE for simulation studies.}

\item{show_messages}{(Logical) If TRUE, show messages from Stan sampler,
if FALSE, hide messages.}

\item{compute_ll}{(Logical) If TRUE, compute log-likelihood,
if FALSE, do not. This may be useful for cross-validation. This argument is
ignored when:
(i) the full dataset is not provided;
(ii) the method is WB, use WB-cond instead.}

\item{acov_mat}{(Optional) Asymptotic variance matrix of
lower triangular half (column-order) of the correlation matrix
to be used for correlation structure analysis.
This parameter is useful if importing polychoric or meta-analytic
SEM pooled correlation matrix.}

\item{ret_data_list}{(LOGICAL)
If TRUE, returns the \code{data_list} and \code{prior} objects, see example.
If FALSE (default), fits the model given user inputs.}
}
\value{
An object of \code{\link{mbsem-class}}
}
\description{
Fit Bayesian path anlaysis models with tests of
conditional independence.
}
\examples{
\dontrun{
minorbpa("x3 ~ x1 + x2\n x4 ~ x3 + x1", HS)
}
}
