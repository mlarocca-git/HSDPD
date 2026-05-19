make_covs_design_fixture <- function(covariates = 0) {
  rr_y <- matrix(
    c(1, 2, 4, 7, 11, 16,
      2, 3, 5, 8, 12, 17,
      3, 4, 6, 9, 13, 18),
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

test_that("design-aware covariance helper matches numeric path without covariates", {
  fit_sdpd_covs_from_design <- getFromNamespace(
    ".fit_sdpd_covs_from_design",
    "HSDPD"
  )
  design <- make_covs_design_fixture(covariates = 0)

  design_covs <- fit_sdpd_covs_from_design(design)
  numeric_covs <- fit_sdpd_covs(
    series = design$series,
    x = design$x,
    kk = design$kk,
    px_neighbors = design$px_neighbors,
    nn = design$nn,
    pp = design$pp
  )

  expect_equal(names(design_covs), names(numeric_covs))
  expect_equal(design_covs$index, numeric_covs$index)
  expect_equal(design_covs$cov11, numeric_covs$cov11)
  expect_equal(design_covs$cov12, numeric_covs$cov12)
  expect_equal(design_covs$cov_x, numeric_covs$cov_x)
})

test_that("design-aware covariance helper matches numeric path with one covariate", {
  fit_sdpd_covs_from_design <- getFromNamespace(
    ".fit_sdpd_covs_from_design",
    "HSDPD"
  )
  design <- make_covs_design_fixture(covariates = "x1")

  design_covs <- fit_sdpd_covs_from_design(design)
  numeric_covs <- fit_sdpd_covs(
    series = design$series,
    x = design$x,
    kk = design$kk,
    px_neighbors = design$px_neighbors,
    nn = design$nn,
    pp = design$pp
  )

  expect_equal(names(design_covs), names(numeric_covs))
  expect_equal(design_covs$index, numeric_covs$index)
  expect_equal(design_covs$cov11, numeric_covs$cov11)
  expect_equal(design_covs$cov12, numeric_covs$cov12)
  expect_equal(design_covs$cov_x, numeric_covs$cov_x)
})

test_that("design-aware covariance helper passes na_covs through", {
  fit_sdpd_covs_from_design <- getFromNamespace(
    ".fit_sdpd_covs_from_design",
    "HSDPD"
  )
  design <- make_covs_design_fixture(covariates = "x1")

  design_covs <- fit_sdpd_covs_from_design(
    design,
    na_covs = "complete.obs"
  )
  numeric_covs <- fit_sdpd_covs(
    series = design$series,
    x = design$x,
    kk = design$kk,
    px_neighbors = design$px_neighbors,
    nn = design$nn,
    pp = design$pp,
    na_covs = "complete.obs"
  )

  expect_equal(design_covs, numeric_covs)
})
