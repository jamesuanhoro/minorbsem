---
title: 'minorbsem: An R package for structural equation models that account for the
  influence of minor factors'
tags:
- "Bayesian-statistics"
- "latent-variable-models"
- "structural-equation-modeling"
- "psychometrics"
- "meta-analytic-SEM"
date: "05/03/2022"
affiliations:
- name: Department of Educational Psychology, University or North Texas, USA
  index: 1
authors:
- name: James Ohisei Uanhoro
  orcid: "0000-0002-4843-927X"
  affiliation: 1
bibliography: inst/REFERENCES.bib
---

The minorbsem (**minor B**ayesian **s**tructural **e**quation **m**odels) R package is a solution to a common issue in structural equation modeling (SEM).
SEMs are hypothesized models that capture the inter-relations between factors (latent variables) and the
relations between these factors and observed items.
However, covariance structures in SEM are may be influenced by minor factors, which we cannot predict a priori [@maccallum_representing_1991].
This can result in SEMs being unable to reject the null hypothesis of no model misspecification. To address this issue, minorbsem estimates all residual covariances between observed items to account for the influence of minor factors [@uanhoro_modeling_2023].
minorbsem implements a variety of shrinkage options to estimate these residual covariances.
Additionally, the minorbsem returns the magnitude of the influence of minor factors.
By estimating both minor and major factors simultaneously, minorbsem provides a more accurate reflection of the uncertainty due to the influence of minor factors on the parameters of interest in SEMs.

## Available models and supported data types

minorbsem also fits random-effects meta-analytic confirmatory factor analysis (CFA) models that capture the influence of minor factors [@uanhoro_hierarchical_2022].

However, minorbsem has some limitations, such as assuming multivariate normal data and fitting a limited number of model configurations:

- CFA, allowing cross-loadings (which may be automatically estimated),
correlated errors terms, and fully oblique or orthogonal factors
(useful for fitting bifactor models)
- SEMs allowing latent regressions (only), cross-loadings, and correlated error
terms.

At the moment, minorbsem cannot be used to fit
multiple-indicator multiple-cause models, multi-group models,
multilevel models, or models with specially constrained parameters
(e.g. setting two parameters equal).

## Key dependencies

minorbsem relies on the Stan project [@Carpenter2017] for fitting Bayesian models,
and interfaces with Stan via the `cmdstanr` package.
Installation instructions for both Stan and minorbsem are available in the minorbsem [README file](README.md/#installation).

minorbsem relies on the `lavaan` package [@Rosseel2012] for data processing.
Moreover, minorbsem models are specified using lavaan-style syntax which allows for easy model specification.

minorbsem also relies on `Rcpp` to speed up costly operations [@rcpp],
the `kableExtra` package for printing results [@kable_ex],
the `posterior` package for accessing Stan model results [@posterior],
and `ggplot2` for plotting [@ggplot2].

## Software with similar functionality

The Mplus package [@muthen8] is also able to estimate Bayesian models that account for the
influence of minor factors. Mplus can estimate a large number of model configurations, but it requires users to manually specify the size of the influence of minor factors, which can lead to inadequate uncertainty estimates for parameters of interest [@uanhoro_modeling_2023].
Additionally, Mplus is not a free software option.

`LAWBL` [@LAWBL] is an R package with similar functionality to minorbsem. LAWBL accommodates more data types
than minorbsem e.g. binary data models are available. However, LAWBL only estimates CFA models and
its syntax is not a commonly used syntax. Additionally, LAWBL does not return the magnitude of the
influence of minor factors.

Finally, neither Mplus nor LAWBL allow the range of options for modeling the influence of minor factors
that are implemented in minorbsem. Nor do they offer options for fitting random-effects meta-analytic SEMs that also capture the
influence of minor factors.

# References
