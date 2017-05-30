
<!-- README.md is generated from README.Rmd. Please edit that file -->
Solving and simulating a class of LRE models
============================================

This package solves for the recursive representation of the stable solution to a system of linear difference equations. Currently, this is done only for the autonomous case with the assumption of no exogenous shocks.

Inputs are two square matrices E and A and a natural number n where E and A are the coefficient matrices of the difference equation Ex(t+1) = Ax(t).

The output of the package are (currently) two functions g and h for the decision rule and the law of motion respectively.

To install run:

``` r
devtools::install_github("nuritovbek/lrem")
```

To use the package, call the lre\_auto function with the coefficient matrices of your model and a number of predetermined variables in n:

``` r
library(lrem)
lre_auto(A, E, n)
```
