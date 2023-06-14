
# FitDynMix

The goal of FitDynMix is to estimate a dynamic lognormal - Generalized
Pareto mixture via the Approximate Maximum Likelihood and the
Cross-Entropy methods.

## Installation

You can install the development version of FitDynMix from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("marco-bee/FitDynMix")
```

## Example

This is a basic example of estimation via AMLE:

``` r
library(FitDynMix)
## basic example code
k <- 5000
epsilon <- .02
bootreps <- 2
res = AMLEfit(Metro2019, epsilon, k, bootreps)
```
