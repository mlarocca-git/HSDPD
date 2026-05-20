make_residual_plot_results <- function() {
  list(
    resid = data.frame(
      t0 = c(0, 0, 0),
      t1 = c(-1, 0, 1),
      t2 = c(1, 2, 3),
      check.names = FALSE
    ),
    lon = c(10, 11, 12),
    lat = c(45, 46, 47)
  )
}

test_that("plot_sdpd_residuals_stat is exported", {
  expect_true("plot_sdpd_residuals_stat" %in% getNamespaceExports("HSDPD"))
  expect_type(HSDPD::plot_sdpd_residuals_stat, "closure")
})

test_that("plot_sdpd_residuals_stat returns a ggplot object", {
  results <- make_residual_plot_results()

  plot <- plot_sdpd_residuals_stat(
    results = results,
    main = "Temporal mean of residuals",
    statistic = mean,
    mid_value = 0
  )

  expect_s3_class(plot, "ggplot")
})

test_that("legacy residual statistic plot wrapper still works", {
  results <- make_residual_plot_results()

  expect_warning(
    plot <- fun.plot.stat.discrete.RESIDs(
      results,
      main = "Temporal mean of residuals",
      statistic = mean,
      mid_value = 0
    ),
    "deprecated"
  )

  expect_s3_class(plot, "ggplot")
})

test_that("new and legacy residual statistic plot wrappers are equivalent", {
  results <- make_residual_plot_results()

  new_plot <- plot_sdpd_residuals_stat(
    results = results,
    main = "Temporal mean of residuals",
    statistic = mean,
    mid_value = 0
  )

  legacy_plot <- suppressWarnings(
    fun.plot.stat.discrete.RESIDs(
      results,
      main = "Temporal mean of residuals",
      statistic = mean,
      mid_value = 0
    )
  )

  expect_equal(class(new_plot), class(legacy_plot))
  expect_equal(ggplot2::ggplot_build(new_plot)$data, ggplot2::ggplot_build(legacy_plot)$data)
  expect_equal(new_plot$labels, legacy_plot$labels)
})
