make_mean_equation_design_fixture <- function() {
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
  model <- build_sdpd_model(
    lambda_0 = FALSE,
    lambda_1 = TRUE,
    lambda_2 = FALSE,
    covariates = 0,
    fixed_effects = FALSE
  )

  series <- read_data_from_dataframe(
    px = c("1", "2", "3"),
    rr_y = rr_y,
    model = model,
    ww_index = ww_index,
    ww_values = ww_values,
    px_neighbors = list(index = ww_index)
  )

  as_sdpd_design <- getFromNamespace(".as_sdpd_design", "HSDPD")
  as_sdpd_design(series = series, model = model)
}

make_mean_equation_fitted_result <- function(design) {
  fit_sdpd_covs_from_design <- getFromNamespace(
    ".fit_sdpd_covs_from_design",
    "HSDPD"
  )
  fit_sdpd_coefficients_from_design <- getFromNamespace(
    ".fit_sdpd_coefficients_from_design",
    "HSDPD"
  )
  fit_sdpd_series_from_design <- getFromNamespace(
    ".fit_sdpd_series_from_design",
    "HSDPD"
  )

  covs <- fit_sdpd_covs_from_design(design)
  coeff_hat <- fit_sdpd_coefficients_from_design(design, covs)$coeff_hat
  fit_sdpd_series_from_design(design, coeff_hat)
}

expect_mean_equation_design_equal_numeric_path <- function(design,
                                                          fitted_result) {
  fit_sdpd_mean_equation_from_design <- getFromNamespace(
    ".fit_sdpd_mean_equation_from_design",
    "HSDPD"
  )

  design_mean <- fit_sdpd_mean_equation_from_design(
    design = design,
    fitted_result = fitted_result
  )
  numeric_mean <- fit_sdpd_mean_equation_model(
    result = fitted_result,
    ww = design$ww,
    time_weights = design$time_weights,
    index_weights = design$index_weights
  )

  expect_equal(names(design_mean), names(numeric_mean))
  expect_equal(design_mean, numeric_mean)
}

test_that("design-aware mean-equation helper matches numeric path", {
  design <- make_mean_equation_design_fixture()
  fitted_result <- make_mean_equation_fitted_result(design)

  expect_mean_equation_design_equal_numeric_path(design, fitted_result)
})

test_that("design-aware mean-equation helper matches numeric path with index weights", {
  design <- make_mean_equation_design_fixture()
  design$index_weights <- c("1" = 1, "2" = 2, "3" = 3)
  fitted_result <- make_mean_equation_fitted_result(design)

  expect_mean_equation_design_equal_numeric_path(design, fitted_result)
})

test_that("design-aware mean-equation helper matches numeric path with numeric time weights", {
  design <- make_mean_equation_design_fixture()
  design$time_weights <- seq_len(design$nn)
  fitted_result <- make_mean_equation_fitted_result(design)

  expect_mean_equation_design_equal_numeric_path(design, fitted_result)
})

test_that("design-aware mean-equation helper matches numeric path with factor time weights", {
  design <- make_mean_equation_design_fixture()
  design$index_weights <- c("1" = 1, "2" = 2, "3" = 3)
  design$time_weights <- factor(
    c("early", "early", "mid", "mid", "late", "late"),
    levels = c("early", "mid", "late")
  )
  fitted_result <- make_mean_equation_fitted_result(design)

  expect_mean_equation_design_equal_numeric_path(design, fitted_result)
})
