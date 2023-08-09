
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tseriesEntropy

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/tseriesEntropy)](https://CRAN.R-project.org/package=tseriesEntropy)
[![CRAN_download](http://cranlogs.r-pkg.org/badges/tseriesEntropy)](https://cran.r-project.org/package=tseriesEntropy)
[![CRAN_download_total](http://cranlogs.r-pkg.org/badges/grand-total/tseriesEntropy)](https://cran.r-project.org/package=tseriesEntropy)
<!-- badges: end -->

The R package `tseriesEntropy` implements an entropy measure of
dependence based on the Bhattacharya-Hellinger-Matusita distance. It can
be used as a (nonlinear) autocorrelation/crosscorrelation function for
continuous and categorical time series. The package includes tests for
serial and cross dependence and nonlinearity based on it. Some routines
have a parallel version that can be used in a multicore/cluster
environment. The package makes use of S4 classes.

## Authors

- [Simone Giannerini, University of
  Bologna](https://www.simonegiannerini.net)

## References

Giannerini S., Maasoumi E., Bee Dagum E., (2015), [Entropy testing for
nonlinear serial dependence in time
series](https://doi.org/10.1093/biomet/asv007), *Biometrika*,
**102(3)**, 661–675.

Giannerini S, Goracci G. (2023) [Entropy-Based Tests for Complex
Dependence in Economic and Financial Time Series with the R Package
tseriesEntropy](https://doi.org/10.3390/math11030757), *Mathematics*,
**11(3):757**.

Granger C. W. J., Maasoumi E., Racine J., (2004) A dependence metric for
possibly nonlinear processes. *Journal of Time Series Analysis*,
**25(5)**, 649–669.

## Installation

You can install the stable version on
[CRAN](https://cran.r-project.org/package=tseriesEntropy):

``` r
install.packages('tseriesEntropy')
```

You can install the development version of tseriesEntropy from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sgiannerini/tseriesEntropy")
```

## License

This package is free and open source software, licensed under GPL.
