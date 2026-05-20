test_that("primary workflow functions remain exported", {
  exports <- getNamespaceExports("HSDPD")

  expect_setequal(
    intersect(
      c(
        "build_sdpd_model",
        "build_sdpd_series",
        "read_data_from_dataframe",
        "read_data_from_raster",
        "fit_sdpd_model",
        "test_sdpd_model",
        "plot_sdpd_residuals_stat",
        "increment_partial_sum",
        "summarise_sdpd_simulation"
      ),
      exports
    ),
    c(
      "build_sdpd_model",
      "build_sdpd_series",
      "read_data_from_dataframe",
      "read_data_from_raster",
      "fit_sdpd_model",
      "test_sdpd_model",
      "plot_sdpd_residuals_stat",
      "increment_partial_sum",
      "summarise_sdpd_simulation"
    )
  )
})

test_that("low-level estimation helpers remain exported", {
  exports <- getNamespaceExports("HSDPD")

  expect_setequal(
    intersect(
      c(
        "fit_sdpd_covs",
        "fit_sdpd_coefficients",
        "fit_second_stage",
        "fit_sdpd_series",
        "fit_sdpd_mean_equation_model"
      ),
      exports
    ),
    c(
      "fit_sdpd_covs",
      "fit_sdpd_coefficients",
      "fit_second_stage",
      "fit_sdpd_series",
      "fit_sdpd_mean_equation_model"
    )
  )
})

test_that("deprecated compatibility wrappers remain exported", {
  exports <- getNamespaceExports("HSDPD")

  expect_setequal(
    intersect(
      c(
        "build.sdpd.model",
        "build.sdpd.series",
        "fit.sdpd.model",
        "fit.sdpd.covs",
        "fit.sdpd.coefficients",
        "fit.2nd.stage",
        "read.data.from.dataframe",
        "read.data.from.raster",
        "test.sdpd.model"
      ),
      exports
    ),
    c(
      "build.sdpd.model",
      "build.sdpd.series",
      "fit.sdpd.model",
      "fit.sdpd.covs",
      "fit.sdpd.coefficients",
      "fit.2nd.stage",
      "read.data.from.dataframe",
      "read.data.from.raster",
      "test.sdpd.model"
    )
  )
})

test_that("registered plot S3 methods remain available", {
  expect_type(getS3method("plot", "sdpd.model", optional = TRUE), "closure")
  expect_type(getS3method("plot", "sdpd.series", optional = TRUE), "closure")
  expect_type(getS3method("plot", "sdpd.estimates", optional = TRUE), "closure")
  expect_type(getS3method("plot", "sdpd.boot.series", optional = TRUE), "closure")
  expect_type(getS3method("plot", "sdpd.series.FITs", optional = TRUE), "closure")
})
