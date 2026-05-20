test_that("build_sdpd_series stores nn for simulated checked series", {
  set.seed(123)

  sdpd_model <- build_sdpd_model(
    endogenous = "yy",
    covariates = c("x1", "x2"),
    sim = TRUE,
    pp = 20,
    d_i = 5,
    check = FALSE
  )

  sdpd_series <- expect_no_error(
    build_sdpd_series(
      model = sdpd_model,
      sim = TRUE,
      nn = 10,
      check = TRUE
    )
  )

  expect_null(sdpd_series$error)
  expect_equal(sdpd_series$nn, 10)

  series_check <- expect_no_error(
    check_sdpd_series(series = sdpd_series, model = sdpd_model)
  )

  expect_length(series_check$errors, 0)
  expect_equal(series_check$nn, 10)
})

test_that("deprecated build.sdpd.series wrapper keeps simulated nn behavior", {
  set.seed(123)

  sdpd_model <- build_sdpd_model(
    endogenous = "yy",
    covariates = "x1",
    sim = TRUE,
    pp = 10,
    d_i = 5,
    check = FALSE
  )

  expect_warning(
    sdpd_series <- build.sdpd.series(
      model = sdpd_model,
      SIM = TRUE,
      nn = 8,
      check = TRUE
    ),
    "deprecated"
  )

  expect_null(sdpd_series$error)
  expect_equal(sdpd_series$nn, 8)
})
