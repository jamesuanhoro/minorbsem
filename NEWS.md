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
