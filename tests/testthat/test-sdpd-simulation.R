make_partial_simulation_sum <- function() {
  px <- c("1", "2")
  coefficient_names <- c("lambda_1", "x1")
  residual_names <- c("mean", "sd", "pvalue.LB", "pvalue.JB", "max.eigenA")

  coefficient_zero <- matrix(
    0,
    nrow = length(px),
    ncol = length(coefficient_names),
    dimnames = list(px, coefficient_names)
  )
  residual_zero <- matrix(
    0,
    nrow = length(px),
    ncol = length(residual_names),
    dimnames = list(px, residual_names)
  )

  list(
    estimation.bias = data.frame(coefficient_zero, check.names = FALSE),
    estimation.mse = data.frame(coefficient_zero, check.names = FALSE),
    resid.sim = data.frame(residual_zero, check.names = FALSE),
    contatore1 = data.frame(coefficient_zero, check.names = FALSE),
    contatore2 = data.frame(residual_zero, check.names = FALSE)
  )
}

test_that("increment_partial_sum updates coefficient sums and counters", {
  partial_sum <- make_partial_simulation_sum()
  fit_results <- data.frame(
    coeff_hat = data.frame(
      lambda_1 = c(0.3, NA),
      x1 = c(1.1, 1.8),
      row.names = c("1", "2"),
      check.names = FALSE
    ),
    mean_resid = c(0.1, NA),
    sd_resid = c(1.2, 1.4),
    pvalue_lb = c(0.4, 0.5),
    pvalue_jb = c(0.6, NA),
    max_eigen_a = c(0.9, 0.8),
    row.names = c("1", "2")
  )
  simulated_coefficients <- data.frame(
    lambda_1 = c(0.1, 0.2),
    x1 = c(1.0, 2.0),
    row.names = c("1", "2"),
    check.names = FALSE
  )

  result <- increment_partial_sum(
    partial_sum = partial_sum,
    fit_results = fit_results,
    simulated_coefficients = simulated_coefficients
  )

  expect_equal(result$estimation.bias["1", "lambda_1"], 0.2)
  expect_equal(result$estimation.bias["1", "x1"], 0.1)
  expect_equal(result$estimation.bias["2", "x1"], -0.2)
  expect_equal(result$estimation.bias["2", "lambda_1"], 0)
  expect_equal(result$estimation.mse["1", "lambda_1"], 0.04)
  expect_equal(result$estimation.mse["2", "x1"], 0.04)
  expect_equal(result$contatore1["1", "lambda_1"], 1)
  expect_equal(result$contatore1["2", "lambda_1"], 0)
  expect_equal(result$contatore1["2", "x1"], 1)
})

test_that("increment_partial_sum updates residual diagnostics with legacy names", {
  partial_sum <- make_partial_simulation_sum()
  fit_results <- data.frame(
    coeff_hat = data.frame(
      lambda_1 = c(0.3, 0.4),
      x1 = c(1.1, 1.8),
      row.names = c("1", "2"),
      check.names = FALSE
    ),
    mean_resid = c(0.1, NA),
    sd_resid = c(1.2, 1.4),
    pvalue_lb = c(0.4, 0.5),
    pvalue_jb = c(0.6, NA),
    max_eigen_a = c(0.9, 0.8),
    row.names = c("1", "2")
  )
  simulated_coefficients <- data.frame(
    lambda_1 = c(0.1, 0.2),
    x1 = c(1.0, 2.0),
    row.names = c("1", "2"),
    check.names = FALSE
  )

  result <- increment_partial_sum(
    partial_sum = partial_sum,
    fit_results = fit_results,
    simulated_coefficients = simulated_coefficients
  )

  expect_equal(result$resid.sim["1", "mean"], 0.1)
  expect_equal(result$resid.sim["2", "mean"], 0)
  expect_equal(result$resid.sim["2", "sd"], 1.4)
  expect_equal(result$resid.sim["1", "pvalue.LB"], 0.4)
  expect_equal(result$resid.sim["1", "pvalue.JB"], 0.6)
  expect_equal(result$resid.sim["2", "max.eigenA"], 0.8)
  expect_equal(result$contatore2["1", "mean"], 1)
  expect_equal(result$contatore2["2", "mean"], 0)
  expect_equal(result$contatore2["2", "pvalue.JB"], 0)
})

test_that("summarise_sdpd_simulation returns compact column means", {
  partial_sum <- make_partial_simulation_sum()
  partial_sum$estimation.bias[1, "lambda_1"] <- 0.2
  partial_sum$estimation.mse[2, "x1"] <- 0.04
  partial_sum$resid.sim[1, "mean"] <- 0.1

  result <- summarise_sdpd_simulation(partial_sum)

  expect_named(result, c("component", "statistic", "value"))
  expect_true(all(c("estimation.bias", "estimation.mse", "resid.sim") %in% result$component))
  expect_equal(
    result$value[result$component == "estimation.bias" & result$statistic == "lambda_1"],
    0.1
  )
})
