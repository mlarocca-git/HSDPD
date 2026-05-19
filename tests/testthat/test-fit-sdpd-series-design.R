make_series_design_fixture <- function(covariates = 0) {
  rr_y <- matrix(
    c(1, 2, 4, 7, 11, 16,
      2, 4, 7, 11, 16, 22,
      1, 3, 6, 10, 15, 21),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), paste0("t", 1:6))
  )

  ww_index <- matrix(
    c(2, 3,
      1, 3,
      1, 2),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), NULL)
  )
  ww_values <- matrix(
    0.5,
    nrow = 3,
    ncol = 2,
    dimnames = list(c("1", "2", "3"), NULL)
  )

  if (identical(covariates, 0)) {
    rr_xx <- NULL
  } else {
    rr_xx <- array(
      c(2, 4, 6, 8, 10, 12,
        1, 3, 5, 7, 9, 11,
        3, 5, 7, 9, 11, 13),
      dim = c(1, 3, 6),
      dimnames = list(
        "x1",
        c("1", "2", "3"),
        paste0("t", 1:6)
      )
    )
  }

  model <- build_sdpd_model(
    lambda_0 = FALSE,
    lambda_1 = TRUE,
    lambda_2 = FALSE,
    covariates = covariates,
    fixed_effects = FALSE
  )

  series <- read_data_from_dataframe(
    px = c("1", "2", "3"),
    rr_y = rr_y,
    rr_xx = rr_xx,
    model = model,
    ww_index = ww_index,
    ww_values = ww_values,
    px_neighbors = list(index = ww_index)
  )

  as_sdpd_design <- getFromNamespace(".as_sdpd_design", "HSDPD")
  as_sdpd_design(series = series, model = model)
}

make_series_coefficients <- function(design) {
  fit_sdpd_covs_from_design <- getFromNamespace(
    ".fit_sdpd_covs_from_design",
    "HSDPD"
  )
  fit_sdpd_coefficients_from_design <- getFromNamespace(
    ".fit_sdpd_coefficients_from_design",
    "HSDPD"
  )

  covs <- fit_sdpd_covs_from_design(design)
  fit_sdpd_coefficients_from_design(design, covs)$coeff_hat
}

expect_series_design_equal_numeric_path <- function(design,
                                                    coeff_hat,
                                                    resids = NULL,
                                                    markovian = TRUE,
                                                    num_steps = 10) {
  fit_sdpd_series_from_design <- getFromNamespace(
    ".fit_sdpd_series_from_design",
    "HSDPD"
  )

  design_fit <- fit_sdpd_series_from_design(
    design = design,
    coeff_hat = coeff_hat,
    resids = resids,
    markovian = markovian,
    num_steps = num_steps
  )
  numeric_fit <- fit_sdpd_series(
    data_series = design$series,
    ww = design$ww,
    x_centered = design$x,
    model = design$model,
    coeff_hat = coeff_hat,
    time_effects = design$time_effects,
    resids = resids,
    markovian = markovian,
    num_steps = num_steps
  )

  expect_equal(names(design_fit), names(numeric_fit))
  expect_equal(design_fit$series, numeric_fit$series)
  expect_equal(design_fit$fitted, numeric_fit$fitted)
  expect_equal(design_fit$resid, numeric_fit$resid)
}

test_that("design-aware fitted-series helper matches numeric path without covariates", {
  design <- make_series_design_fixture(covariates = 0)
  coeff_hat <- make_series_coefficients(design)

  expect_series_design_equal_numeric_path(design, coeff_hat)
})

test_that("design-aware fitted-series helper matches numeric path with one covariate", {
  design <- make_series_design_fixture(covariates = "x1")
  coeff_hat <- make_series_coefficients(design)

  expect_series_design_equal_numeric_path(design, coeff_hat)
})

test_that("design-aware fitted-series helper matches numeric residual update path", {
  design <- make_series_design_fixture(covariates = 0)
  coeff_hat <- make_series_coefficients(design)
  resids <- matrix(
    c(0.1, 0.2, 0.1, 0.2, 0.1, 0.2,
      0.2, 0.1, 0.2, 0.1, 0.2, 0.1,
      0.1, 0.1, 0.2, 0.2, 0.1, 0.1),
    nrow = 3,
    byrow = TRUE,
    dimnames = dimnames(design$series)
  )

  expect_series_design_equal_numeric_path(
    design = design,
    coeff_hat = coeff_hat,
    resids = resids,
    markovian = FALSE,
    num_steps = 2
  )
})
