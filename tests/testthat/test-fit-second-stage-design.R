make_second_stage_design_fixture <- function(fixed_effects = FALSE) {
  rr_y <- matrix(
    c(1, 3, 6, 10, 15, 21, 28, 36,
      2, 5, 9, 14, 20, 27, 35, 44,
      1, 4, 8, 13, 19, 26, 34, 43),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), paste0("t", 1:8))
  )
  rr_xx <- array(
    c(1, 2, 4, 7, 11, 16, 22, 29,
      2, 4, 7, 11, 16, 22, 29, 37,
      1, 3, 6, 10, 15, 21, 28, 36,
      2, 5, 10, 17, 26, 37, 50, 65,
      3, 7, 13, 21, 31, 43, 57, 73,
      1, 4, 9, 16, 25, 36, 49, 64),
    dim = c(2, 3, 8),
    dimnames = list(
      c("x1", "x2"),
      c("1", "2", "3"),
      paste0("t", 1:8)
    )
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
    covariates = c("x1", "x2"),
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

expect_second_stage_design_equal_numeric_path <- function(design) {
  fit_second_stage_from_design <- getFromNamespace(
    ".fit_second_stage_from_design",
    "HSDPD"
  )

  coeff_names <- c(
    names(design$model$lambda_coeffs)[design$model$lambda_coeffs],
    names(design$model$beta_coeffs)[design$model$beta_coeffs],
    "fixed_effects"[design$model$fixed_effects]
  )
  coeff_hat <- matrix(
    0,
    nrow = design$pp,
    ncol = length(coeff_names),
    dimnames = list(design$unit_index, coeff_names)
  )
  coeff_hat[, "lambda_1"] <- c(0.1, 0.2, 0.15)

  design_fit <- fit_second_stage_from_design(
    design = design,
    coeff_hat = coeff_hat
  )
  numeric_fit <- fit_second_stage(
    data_series = design$series,
    x_centered = design$x,
    ww = design$ww,
    model = design$model,
    coeff_hat = coeff_hat
  )

  expect_equal(design_fit, numeric_fit)
  expect_equal(rownames(design_fit), design$unit_index)
}

test_that("design-aware second-stage helper matches numeric path", {
  design <- make_second_stage_design_fixture(fixed_effects = FALSE)

  expect_second_stage_design_equal_numeric_path(design)
})

test_that("design-aware second-stage helper matches numeric path with fixed effects", {
  design <- make_second_stage_design_fixture(fixed_effects = TRUE)

  expect_second_stage_design_equal_numeric_path(design)

  fit_second_stage_from_design <- getFromNamespace(
    ".fit_second_stage_from_design",
    "HSDPD"
  )

  coeff_hat <- matrix(
    0,
    nrow = design$pp,
    ncol = 4,
    dimnames = list(design$unit_index, c("lambda_1", "x1", "x2", "fixed_effects"))
  )
  coeff_hat[, "lambda_1"] <- c(0.1, 0.2, 0.15)
  design_fit <- fit_second_stage_from_design(design, coeff_hat)

  expect_true("fixed_effects" %in% colnames(design_fit))
})
