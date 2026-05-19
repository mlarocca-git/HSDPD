make_coefficients_design_fixture <- function(covariates = 0,
                                             fixed_effects = FALSE) {
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
    fixed_effects = fixed_effects
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

expect_design_coefficients_equal_numeric_path <- function(design) {
  fit_sdpd_covs_from_design <- getFromNamespace(
    ".fit_sdpd_covs_from_design",
    "HSDPD"
  )
  fit_sdpd_coefficients_from_design <- getFromNamespace(
    ".fit_sdpd_coefficients_from_design",
    "HSDPD"
  )

  covs <- fit_sdpd_covs_from_design(design)
  design_fit <- fit_sdpd_coefficients_from_design(design, covs = covs)
  numeric_fit <- fit_sdpd_coefficients(
    ww = design$ww,
    covs = covs,
    mu = design$mu,
    model = design$model
  )

  expect_equal(names(design_fit), names(numeric_fit))
  expect_equal(design_fit$coeff_hat, numeric_fit$coeff_hat)
  expect_equal(rownames(design_fit$coeff_hat), design$unit_index)
}

test_that("design-aware coefficient helper matches numeric path without covariates", {
  design <- make_coefficients_design_fixture(
    covariates = 0,
    fixed_effects = FALSE
  )

  expect_design_coefficients_equal_numeric_path(design)
})

test_that("design-aware coefficient helper matches numeric path with one covariate", {
  design <- make_coefficients_design_fixture(
    covariates = "x1",
    fixed_effects = FALSE
  )

  expect_design_coefficients_equal_numeric_path(design)
})

test_that("design-aware coefficient helper matches numeric path with fixed effects", {
  design <- make_coefficients_design_fixture(
    covariates = 0,
    fixed_effects = TRUE
  )

  expect_design_coefficients_equal_numeric_path(design)
  expect_true("fixed_effects" %in% colnames(
    getFromNamespace(".fit_sdpd_coefficients_from_design", "HSDPD")(
      design,
      covs = getFromNamespace(".fit_sdpd_covs_from_design", "HSDPD")(design)
    )$coeff_hat
  ))
})
