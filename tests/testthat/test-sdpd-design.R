make_dataframe_series_fixture <- function(px = c("1", "2")) {
  rr_y <- data.frame(
    latitude = c(40.1, 40.2, 40.3),
    longitude = c(14.1, 14.2, 14.3),
    t1 = c(1, 2, 3),
    t2 = c(2, 3, 4),
    t3 = c(3, 4, 5),
    row.names = c("1", "2", "3")
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
  groups <- data.frame(
    COD = c(1, 1, 2),
    LABEL = c("north", "north", "south"),
    row.names = c("1", "2", "3")
  )
  model <- build_sdpd_model(covariates = 0, fixed_effects = TRUE)

  list(
    series = read_data_from_dataframe(
      px = px,
      rr_y = rr_y,
      rr_groups = groups,
      model = model,
      ww_index = ww_index,
      ww_values = ww_values,
      px_neighbors = list(index = ww_index)
    ),
    model = model
  )
}

test_that("sdpd_design is built from dataframe-based input", {
  as_sdpd_design <- getFromNamespace(".as_sdpd_design", "HSDPD")
  fixture <- make_dataframe_series_fixture()

  design <- as_sdpd_design(
    series = fixture$series,
    model = fixture$model
  )

  expect_s3_class(design, "sdpd_design")
  expect_equal(design$errors, character())
  expect_true(is.matrix(design$series))
  expect_true(is.numeric(design$series))
  expect_equal(design$unit_index, c("1", "2", "3"))
  expect_equal(design$time_index, c("t1", "t2", "t3"))
  expect_equal(design$pp, 3)
  expect_equal(design$nn, 3)
  expect_equal(design$kk, 0)
  expect_true(is.matrix(design$ww))
  expect_true(is.numeric(design$ww))
  expect_equal(rownames(design$ww), design$unit_index)
  expect_equal(colnames(design$ww), design$unit_index)
  expect_equal(length(design$mu), design$pp)
  expect_equal(length(design$time_effects), design$nn)
  expect_equal(design$source_px, c("1", "2"))
})

test_that("sdpd_design stores centered dataframe covariates", {
  as_sdpd_design <- getFromNamespace(".as_sdpd_design", "HSDPD")

  rr_y <- matrix(
    c(1, 2, 3,
      2, 3, 4,
      3, 4, 5),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), c("t1", "t2", "t3"))
  )
  rr_xx <- array(
    c(2, 4, 6,
      1, 3, 5,
      3, 5, 7),
    dim = c(1, 3, 3),
    dimnames = list(
      "x1",
      c("1", "2", "3"),
      c("t1", "t2", "t3")
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
    covariates = "x1",
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

  design <- as_sdpd_design(series = series, model = model)

  expect_equal(design$errors, character())
  expect_equal(design$kk, 1)
  expect_true(is.matrix(design$x))
  expect_equal(as.numeric(rowMeans(design$x)), rep(0, design$pp))
})

test_that("sdpd_design reports invalid construction inputs", {
  as_sdpd_design <- getFromNamespace(".as_sdpd_design", "HSDPD")
  new_sdpd_design <- getFromNamespace(".new_sdpd_design", "HSDPD")
  validate_sdpd_design <- getFromNamespace(".validate_sdpd_design", "HSDPD")

  rr_y <- matrix(
    c(1, 2, 3,
      2, 3, 4,
      3, 4, 5),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), c("t1", "t2", "t3"))
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
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE)

  design <- as_sdpd_design(
    series = list(
      series = rr_y,
      ww_index = ww_index,
      ww_values = ww_values,
      px = c("1", "2", "3")
    ),
    model = model
  )

  expect_s3_class(design, "sdpd_design")
  expect_match(paste(design$errors, collapse = " "), "px_neighbors")

  malformed_design <- new_sdpd_design(
    series = matrix(1:4, nrow = 2),
    pp = 3,
    nn = 2,
    kk = 0
  )
  malformed_design <- validate_sdpd_design(malformed_design)

  expect_match(
    paste(malformed_design$errors, collapse = " "),
    "row and column names|pp"
  )
})
