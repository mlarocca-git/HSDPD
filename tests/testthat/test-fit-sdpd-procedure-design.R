make_fit_procedure_series_fixture <- function() {
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
  model <- build_sdpd_model(
    lambda_0 = FALSE,
    lambda_1 = TRUE,
    lambda_2 = FALSE,
    covariates = 0,
    fixed_effects = FALSE
  )

  list(
    model = model,
    series = read_data_from_dataframe(
      px = c("1", "2", "3"),
      rr_y = rr_y,
      model = model,
      ww_index = ww_index,
      ww_values = ww_values,
      px_neighbors = list(index = ww_index)
    )
  )
}

test_that("fit_sdpd_procedure keeps the public result structure", {
  fixture <- make_fit_procedure_series_fixture()

  result <- fit_sdpd_procedure(
    series = fixture$series,
    model = fixture$model
  )

  expect_named(
    result,
    c(
      "px",
      "lon",
      "lat",
      "group",
      "coeff_hat",
      "fitted",
      "resid",
      "mean_equation",
      "time_effects",
      "model",
      "warnings",
      "diagnostics"
    )
  )
  expect_equal(result$px, c("1", "2", "3"))
  expect_true(is.matrix(result$coeff_hat))
  expect_true(is.matrix(result$fitted))
  expect_true(is.matrix(result$resid))
  expect_equal(rownames(result$fitted), c("1", "2", "3"))
  expect_equal(colnames(result$fitted), paste0("t", 1:6))
  expect_equal(result$diagnostics, FALSE)
})

test_that("fit_sdpd_procedure returns errors and warnings for invalid input", {
  fixture <- make_fit_procedure_series_fixture()
  invalid_series <- fixture$series
  invalid_series$px_neighbors <- NULL

  expect_output(
    result <- fit_sdpd_procedure(
      series = invalid_series,
      model = fixture$model
    ),
    regexp = "There are errors in the series:.*px_neighbors"
  )

  expect_named(result, c("errors", "warnings"))
  expect_match(paste(result$errors, collapse = " "), "px_neighbors")
  expect_type(result$warnings, "character")
})

test_that("fit_sdpd_procedure completes with diagnostics enabled", {
  fixture <- make_fit_procedure_series_fixture()

  result <- fit_sdpd_procedure(
    series = fixture$series,
    model = fixture$model,
    check = TRUE
  )

  expect_named(
    result,
    c(
      "px",
      "lon",
      "lat",
      "group",
      "coeff_hat",
      "fitted",
      "resid",
      "mean_equation",
      "time_effects",
      "model",
      "warnings",
      "diagnostics"
    )
  )
  expect_length(result$diagnostics, 1)
  expect_true(is.numeric(result$diagnostics) || is.na(result$diagnostics))
})
