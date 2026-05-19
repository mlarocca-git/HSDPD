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
