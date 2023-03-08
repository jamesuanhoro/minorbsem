---
title: 'minorbsem: An R package for structural equation models that account for the
  influence of minor factors'
tags:
- "Bayesian-statistics"
- "latent-variable-models"
- "structural-equation-modeling"
- psychometrics
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

# minorbsem: An R package for structural equation models that account for the influence of minor factors

minorbsem or **minor B**ayesian **s**tructural **e**quation **m**odels (SEMs)
is an R package that
allows the user fit SEMs under the assumption that minor factors influence the relation between observed
indicators.

SEMs are hypothesized models that capture the inter-relations between latent variables or factors and the
relations between latent variables and observed items.
SEMs are rarely ever able to reject the null hypothesis that there is no misspecification in the hypothesized model.
One explanation for this problem is that covariance structures are influenced by major factors which
we can hypothesize about and minor factors which we cannot predict a-priori, e.g. @maccallum_representing_1991.
The presence of these minor factors explains why SEMs are rarely able to reject the null hypothesis
of no model misspecification.

minorbsem is an implementation of SEMs that accounts for the influence of minor factors, following
the approach in @uanhoro_modeling_2023.
Briefly, minorbsem accounts for the influence of minor factors
by estimating all residual covariances between the observed items.
These residual covariances are estimated with priors that shrink the residual covariances towards zero.
Additionally, the model returns and the model returns the magnitude of the influence of minor factors.
Practically, by estimating the influence of minor factors simultaneously with major factors,
the parameters of interest in the SEM will reflect the uncertainty due to the influence of minor factors.

minorbsem also fits random-effects meta-analytic confirmatory factor analysis models (CFAs)
that capture the influence of minor factors following the approach in @uanhoro_hierarchical_2022.

## Available models and supported data types

minorbsem only fits a limited number of model configurations:

- CFA, allowing cross-loadings (which may be automatically estimated),
correlated errors terms, and fully oblique or orthogonal factors
(useful for fitting bifactor models)
- SEMs allowing latent regressions (only), cross-loadings, and correlated error
terms.

At the moment, one cannot use minorbsem to fit MIMIC, multi-group models,
multilevel models, or models with specially constrained parameters
(e.g. setting two parameters equal).

The meta-analysis models are only for the CFA configurations.
Finally, all data are assumed multivariate normal, i.e. no binary, ordinal models.

## Key dependencies

minorbsem relies on the Stan project [@Carpenter2017] for fitting Bayesian models.
And minorbsem interfaces with Stan via the `cmdstanr` package.
Installation instructions for both Stan and minorbsem are available [here](README/#installation).

minorbsem relies on the `lavaan` package [@Rosseel2012] for data processing.
Moreover, minorbsem models are specified using lavaan-style syntax which allows for easy model specification.

## Similar software



# References
