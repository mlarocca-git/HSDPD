
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HSDPD

<!-- badges: start -->

<!-- badges: end -->

The goal of HSDPD is to …

## Public API

The primary workflow is:

- `build_sdpd_model()` to define an SDP-D model;
- `build_sdpd_series()`, `read_data_from_dataframe()`, or
  `read_data_from_raster()` to prepare data;
- `fit_sdpd_model()` to estimate the model;
- `test_sdpd_model()` and the `plot_sdpd_*()` functions for testing and
  diagnostics.

Advanced low-level helpers such as `fit_sdpd_covs()`,
`fit_sdpd_coefficients()`, `fit_second_stage()`, `fit_sdpd_series()`, and
`fit_sdpd_mean_equation_model()` remain exported for numeric workflows and
internal pipeline compatibility. Most users should prefer `fit_sdpd_model()`.

Deprecated dot-name aliases such as `build.sdpd.series()` and
`fit.sdpd.model()` are kept for backward compatibility and forward to the
snake_case API.

## Installation

You can install the development version of HSDPD from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("mlarocca-git/HSDPD")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(HSDPD)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" alt="" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
