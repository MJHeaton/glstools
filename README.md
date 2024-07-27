
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glstools

<!-- badges: start -->
<!-- badges: end -->

The goal of glsTools is to make analysis of correlated data easier by
providing tools to *correctly* predict and standardize *decorrelated*
residuals from `gls` objects in R. Additional functionality to calculate
Moran Basis functions is also provided.

## Installation

You can install the development version of glstools like so:

``` r
devtools::install_github('https://github.com/MJHeaton/glstools')
```

## Example Usage

``` r
library(glstools)
library(nlme)

## Fit a GLS model
my_model <- gls(y ~ x1 + x2, data = data, correlation = corSpatial(form = ~ x1 + x2))

## Predict new values
new_data <- data.frame(x1 = c(1, 2, 3), x2 = c(4, 5, 6))
my_preds <- predictgls(glsobj = my_model, newdframe = new_data)

## Calculate standardized (and decorrelated) residuals
decorr_resid <- stdres.gls(my_model)
```
