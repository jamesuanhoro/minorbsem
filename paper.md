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
- name: Department of Educational Psychology, University of North Texas, USA
  index: 1
authors:
- name: James Ohisei Uanhoro
  orcid: "0000-0002-4843-927X"
  affiliation: 1
bibliography: inst/REFERENCES.bib
---

# Summary

minorbsem is an R package that allows users to fit structural equation models (SEMs) under the assumption that minor factors influence the relationship between observed indicators. SEMs are models that capture the inter-relations between factors (latent variables) and the relations between these factors and observed items. However, SEMs are rarely able to reject the null hypothesis that there is no misspecification in the hypothesized model. One explanation for this problem is that covariance structures are influenced by major factors that we can hypothesize about, and minor factors that we cannot predict a-priori [@maccallum_representing_1991].

minorbsem accounts for the influence of minor factors by estimating all residual covariances between the observed items [@uanhoro_modeling_2023]. These residual covariances are estimated with priors that shrink the residual covariances towards zero. Additionally, the model returns the magnitude of the influence of minor factors. By estimating the influence of minor factors simultaneously with major factors, the parameters of interest (related to major factors) in the SEM will reflect the uncertainty due to the influence of minor factors.

minorbsem also fits random-effects meta-analytic confirmatory factor analysis models (CFAs) that capture the influence of minor factors [@uanhoro_hierarchical_2022].

# Features

- Fits SEMs allowing latent regressions, cross-loadings, and correlated error terms.
- Fits CFAs allowing cross-loadings, correlated error terms, and fully oblique or orthogonal factors.
- Can estimate the size of the influence of minor factors on the SEM.
- Provides several priors to shrink residual covariances towards zero.
- Fits random-effects meta-analytic CFAs that capture the influence of minor factors.

## Limitations

minorbsem assumes multivariate normal data and only fits a limited number of model configurations. At the moment, it does not support multiple-indicator multiple-cause models, multi-group models, multilevel models, or models with specially constrained parameters.

# Statement of need

@Muthen2012 showed the flexibility of Bayesian SEMs,
and sparked interest in Bayesian SEMs that account for the influence of minor factors.
The work of @wu_quantifying_2015 in frequentist SEMs also addresses the same substantive problem from a different angle.
However, there are limited software implementations of these ideas.

Mplus [@muthen8] is a proprietary software package that is also able to estimate Bayesian models that account for the influence of minor factors. While Mplus can estimate a very large number of model configurations, the user has to specify the size of the influence of minor factors manually, which can lead to inadequate uncertainty estimates for parameters of interest [@uanhoro_modeling_2023]. Additionally, Mplus is not a free software option.

`LAWBL` [@LAWBL] is an R package with similar functionality to minorbsem. LAWBL accommodates more data types
than minorbsem e.g. binary data models are available. However, LAWBL only estimates CFA models and
its syntax is not a commonly used syntax. Additionally, LAWBL does not return the magnitude of the
influence of minor factors.

Finally, neither Mplus nor LAWBL allow the range of options for modeling the influence of minor factors
that are implemented in minorbsem. Nor do they offer options for fitting random-effects meta-analytic SEMs that also capture the
influence of minor factors.

# Dependencies

minorbsem relies on the Stan project [@Carpenter2017] for fitting Bayesian models,
and interfaces with Stan via the `cmdstanr` package.
Installation instructions for both Stan and minorbsem are available in the minorbsem [README](README.md/#installation).

minorbsem relies on the `lavaan` package [@Rosseel2012] for data processing.
Moreover, minorbsem models are specified using lavaan-style syntax which allows for easy model specification.

Additionally, it relies on `Rcpp` to speed up costly operations [@rcpp],
the `kableExtra` package for printing results [@kable_ex],
the `posterior` package for accessing Stan model results [@posterior],
and `ggplot2` for plotting [@ggplot2].

# References
