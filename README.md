
# minorbsem

<!-- badges: start -->
<!-- badges: end -->

It is rare that structural equation models (SEMs) are able to reject the null
hypothesis that there is no misspecification in the model. One explanation for
this problem is that covariance structures are influenced by major factors which
we can hypothesize about and minor factors which we cannot predict a-priori,
e.g. MacCallum & Tucker, 1991.

The goal of minorbsem is to let you fit Bayesian SEMs that estimate the
influence of minor factors on the covariance matrix.

The package only fits a limited number of models:

- CFA, allowing cross-loadings, correlated errors terms, and fully
oblique or orthogonal factors (useful for fitting bifactor models)
- SEMs, allowing latent regressions (only), cross-loadings, and correlated error
terms.

One cannot fit MIMIC, multi-group models, multilevel models, or models with
specially constrained parameters (e.g. setting two parameters equal).

All data are assumed multivariate normal, i.e. no binary, ordinal models.
    
## Installation

The package requires that you **have CmdStan installed on your machine**,
see installation instructions here: https://mc-stan.org/cmdstanr/articles/cmdstanr.html

You can install minorbsem using devtools:

``` r
devtools::install_github("stonegold546/minorbsem")
```

## A reasonably complete demonstration

``` r
library(minorbsem)
# Basic Holzinger-Swineford model
# Expect a summary table output in Viewer
fit_1 <- minorbsem("F1 =~ x1 + x2 + x3
                    F2 =~ x4 + x5 + x6
                    F3 =~ x7 + x8 + x9", HS)

# Save output table to html file, see: ?pretty_print_summary for more options
pretty_print_summary(fit_1, save_html = "baseline_model.html")

# Histogram of parameters, see: ?parameter_hist for arguments
parameter_hist(fit_1)

# Traceplot of parameters, see: ?parameter_trace for arguments
parameter_trace(fit_1)

# Examine all standardized residual covariances
plot_residuals(fit_1)
plot_residuals(fit_1, method = "range")

## Fit a second model and compare both using loo
# A modified Holzinger-Swineford model
fit_2 <- minorbsem("F1 =~ x1 + x2 + x3 + x9
                    F2 =~ x4 + x5 + x6
                    F3 =~ x7 + x8 + x9
                    x2 ~~ x7", HS)

# Compute case wise log-likelihood, exclude minor factor residuals
ll_mat_1 <- casewise_log_likelihood(fit_1, include_residuals = FALSE)
ll_mat_2 <- casewise_log_likelihood(fit_2, include_residuals = FALSE)
chain_id <- posterior::as_draws_df(fit_1@stan_fit)$.chain

# loo for model 1
loo_1 <- loo::loo(
  ll_mat_1,
  r_eff = loo::relative_eff(ll_mat_1, chain_id = chain_id)
)
print(loo_1)

# loo for model 2
loo_2 <- loo::loo(
  ll_mat_2,
  r_eff = loo::relative_eff(ll_mat_2, chain_id = chain_id)
)
print(loo_2)

# Compare both models
print(loo::loo_compare(loo_1, loo_2), simplify = FALSE)
```

## Citations

MacCallum, R. C., & Tucker, L. R. (1991). Representing sources of error in the common-factor model: Implications for theory and practice. _Psychological Bulletin, 109_(3), 502–511. https://doi.org/10.1037/0033-2909.109.3.502
