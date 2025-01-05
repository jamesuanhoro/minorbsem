# minorbsem (development version)

* Added a function (`minorbpa()`) for path analysis (observed variables only), that estimates SRCs based on assumptions of conditional independence within the path analysis model.
* Residual plots for lasso and GDP methods now show the mode as estimate (`modeest::hsm()`) instead of the median, since these methods shrink the mode to 0.
* There are some problems with residual correlations when `correlation = TRUE`

# minorbsem 0.2.15

* Added a function for residual network analysis, doi: 10.1007/s11336-017-9557-x
* Updated citations

# minorbsem 0.2.14

* Added the first CFA tutorial, to come: SEM, MIMIC, correlation structure analysis and modifying priors

# minorbsem 0.2.13

* Now possible to obtain `data_list` from Stan, useful for modifying specific priors.

# minorbsem 0.2.12

* Sort regression results by outcome variable
* Removed equality tests on real numbers in Stan
* Instantiate-forced changes

# minorbsem 0.2.11

* Fixed mistake in calculation of factor residual variance for latent regression models, by transforming `simsem::findFactorResidualVar()` into Stan code.
* Also removed 0 start values for latent regression coefficients.

# minorbsem 0.2.10

* Ensure asymptotic variance of log-correlation matrix is symmetric
* Minor improvement on residual plots

# minorbsem 0.2.9

* Fixed bug when `acov_mat` is supplied to ensure `acov_mat` is correctly ordered based on lavaan object

# minorbsem 0.2.8

* Users can now supply asymptotic variance matrix of correlation matrix as an input

# minorbsem 0.2.7

* Added correlation structure analysis via matrix logarithm transformation

# minorbsem 0.2.6

* Now using single Stan script for all models
  * There are now parameter constraints for latent variable regression models
* Latent regression coefficient parameters are now standardized

# minorbsem 0.2.5

* Eliminated the meta-analysis component, see [bayesianmasem](https://github.com/jamesuanhoro/bayesianmasem/) package
* Using only CmdStan via the instantiate package
* Permits equality constraints and fixed parameters for loadings, residual covariances and error variances.
* Code refactoring
* Log-likelihood now requested with minorbsem call
  * So no more Rcpp

# minorbsem 0.2.4

* JOSS paper added
* Added simple contribution guidelines in startup message and README.

# minorbsem 0.2.3

* Added a Wishart-Wishart method as an alternative to WB-cond.
* Added log-likelihood computations for meta-CFA
* Added control over prior for factor standard deviation parameters

# minorbsem 0.2.2

* Fixed error in implementation of RMSEA calculation in meta-analysis (marginal practical effect)
* Using huxtable to print tables instead of kableExtra

# minorbsem 0.2.1

* Added dependent-samples MASEM (experimental -- CmdStan only)
* Fixed bug that caused target = "cmdstan" to fail when method = "none"

# minorbsem 0.2.0

* RStan implementation now working and with tests
  * Simpler installation
  * Precursor to CRAN
* Fixed row order of variables printed
* Fixed an undiscovered error when user fits model with single chain

# minorbsem 0.1.1

* Added Wishart-based fixed-effects meta-analytic SEM (also permits minor factors)

# minorbsem 0.1.0

* First presentation of package
  * CFA/SEM with minor factors assumed
  * Random-effects meta-analytic permitting minor factor influences
* Added a `NEWS.md` file to track changes to the package.
