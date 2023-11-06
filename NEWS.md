# minorbsem 0.2.5

* Eliminated the meta-analysis component, see [bayesianmasem](https://github.com/jamesuanhoro/bayesianmasem/) package
* Using only CmdStan via the instantiate package
* Permits equality constraints and fixed parameters for loadings, residual covariances and error variances.

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
