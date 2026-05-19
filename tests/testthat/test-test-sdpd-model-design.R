make_test_sdpd_model_checked_data <- function() {
  series <- matrix(
    c(1, 2, 4, 7, 11, 16,
      2, 4, 7, 11, 16, 22,
      1, 3, 6, 10, 15, 21),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), paste0("t", 1:6))
  )
  ww <- matrix(
    c(0, 0.5, 0.5,
      0.5, 0, 0.5,
      0.5, 0.5, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), c("1", "2", "3"))
  )
  px_neighbors <- list(
    index = matrix(
      c(2, 3,
        1, 3,
        1, 2),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("1", "2", "3"), NULL)
    ),
    series_boundary = NULL
  )

  list(
    series = series,
    xx = NULL,
    ww = ww,
    px_neighbors = px_neighbors,
    pp = nrow(series),
    nn = ncol(series),
    kk = 0,
    mu = rowMeans(series)
  )
}

make_test_sdpd_model_result <- function() {
  data <- make_test_sdpd_model_checked_data()
  ww_index <- data$px_neighbors$index
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
  model$ww_index <- ww_index
  model$ww_values <- ww_values

  coeff_hat <- matrix(
    c(0.1, 0.15, 0.2),
    nrow = 3,
    dimnames = list(c("1", "2", "3"), "lambda_1")
  )
  fitted <- fit_sdpd_series(
    data_series = data$series,
    ww = data$ww,
    x_centered = data$xx,
    model = model,
    coeff_hat = coeff_hat
  )

  list(
    model = model,
    coeff_hat = coeff_hat,
    fitted = fitted$fitted,
    resid = fitted$resid,
    data = data
  )
}

test_that("test_sdpd_model design bridge builds a valid design object", {
  test_sdpd_design_from_checked_data <- getFromNamespace(
    ".test_sdpd_design_from_checked_data",
    "HSDPD"
  )
  validate_sdpd_design <- getFromNamespace(".validate_sdpd_design", "HSDPD")
  data <- make_test_sdpd_model_checked_data()
  model <- build_sdpd_model(
    lambda_0 = FALSE,
    lambda_1 = TRUE,
    lambda_2 = FALSE,
    covariates = 0,
    fixed_effects = FALSE
  )

  design <- test_sdpd_design_from_checked_data(data, model)
  design <- validate_sdpd_design(design)

  expect_s3_class(design, "sdpd_design")
  expect_equal(design$errors, character())
  expect_equal(design$series, data$series)
  expect_null(design$x)
  expect_equal(design$unit_index, rownames(data$series))
  expect_equal(design$time_index, colnames(data$series))
  expect_equal(design$pp, data$pp)
  expect_equal(design$nn, data$nn)
  expect_equal(design$kk, data$kk)
  expect_equal(design$ww, data$ww)
  expect_equal(design$px_neighbors, data$px_neighbors)
  expect_equal(design$mu, data$mu)
  expect_identical(design$model, model)
})

test_that("test_sdpd_model design bridge does not mutate checked data", {
  test_sdpd_design_from_checked_data <- getFromNamespace(
    ".test_sdpd_design_from_checked_data",
    "HSDPD"
  )
  data <- make_test_sdpd_model_checked_data()
  data_before <- data
  model <- build_sdpd_model(
    lambda_0 = FALSE,
    lambda_1 = TRUE,
    lambda_2 = FALSE,
    covariates = 0,
    fixed_effects = FALSE
  )

  invisible(test_sdpd_design_from_checked_data(data, model))

  expect_equal(data, data_before)
})

test_that("test_sdpd_model public bootstrap result structure remains stable", {
  set.seed(20240519)
  res_fit <- make_test_sdpd_model_result()

  result <- test_sdpd_model(
    res_fit = res_fit,
    px = rownames(res_fit$data$series),
    n_boot = 21,
    h0 = "zero",
    boot_options = list(
      markovian = FALSE,
      resid = "normal",
      sigma_resid = 0.01,
      boot_plot = FALSE,
      folder = "",
      y_limits = NULL,
      label_index = NULL
    )
  )

  expect_named(
    result,
    c(
      "pvalue",
      "n_boot",
      "h0",
      "diagnostics_model",
      "coeff_hat",
      "diagnostics_coeff_boot",
      "diagnostics_sdevs_tsboot",
      "warnings"
    )
  )
  expect_equal(result$n_boot, 21)
  expect_equal(result$h0, "zero")
  expect_equal(result$coeff_hat, res_fit$coeff_hat)
  expect_equal(dim(result$pvalue), dim(res_fit$coeff_hat))
  expect_equal(rownames(result$pvalue), rownames(res_fit$coeff_hat))
  expect_equal(colnames(result$pvalue), colnames(res_fit$coeff_hat))
  expect_true(all(c(
    "max_mod_eigen_1",
    "max_mod_eigen_2",
    "max_mod_eigen_a"
  ) %in% names(result$diagnostics_model)))
  expect_equal(
    dim(result$diagnostics_coeff_boot),
    c(4, dim(res_fit$coeff_hat))
  )
  expect_equal(
    dimnames(result$diagnostics_coeff_boot)[[2]],
    rownames(res_fit$coeff_hat)
  )
  expect_equal(
    dimnames(result$diagnostics_coeff_boot)[[3]],
    colnames(res_fit$coeff_hat)
  )
  expect_equal(
    dim(result$diagnostics_sdevs_tsboot),
    c(nrow(res_fit$data$series), 2)
  )
})
