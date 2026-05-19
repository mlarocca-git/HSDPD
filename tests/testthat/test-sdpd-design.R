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

make_manual_design <- function(model = build_sdpd_model(
                                 covariates = 0,
                                 fixed_effects = FALSE
                               ),
                               kk = sum(model$beta_coeffs),
                               x = NULL,
                               mu = NULL,
                               time_effects = NULL) {
  new_sdpd_design <- getFromNamespace(".new_sdpd_design", "HSDPD")
  series <- matrix(
    c(1, 2, 3,
      2, 3, 4,
      3, 4, 5),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), c("t1", "t2", "t3"))
  )

  if (is.null(mu)) {
    mu <- stats::setNames(numeric(nrow(series)), rownames(series))
  }

  if (is.null(time_effects)) {
    time_effects <- stats::setNames(numeric(ncol(series)), colnames(series))
  }

  new_sdpd_design(
    series = series,
    x = x,
    unit_index = rownames(series),
    time_index = colnames(series),
    pp = nrow(series),
    nn = ncol(series),
    kk = kk,
    mu = mu,
    time_effects = time_effects,
    model = model
  )
}

validate_manual_design <- function(design) {
  validate_sdpd_design <- getFromNamespace(".validate_sdpd_design", "HSDPD")
  validate_sdpd_design(design)
}

test_that("sdpd_design validates unit and time index alignment", {
  design <- make_manual_design()
  design$unit_index <- c("a", "b", "c")

  design <- validate_manual_design(design)

  expect_match(
    paste(design$errors, collapse = " "),
    "unit_index must match rownames"
  )
})

test_that("sdpd_design validates kk against model covariates", {
  model <- build_sdpd_model(
    covariates = "x1",
    fixed_effects = FALSE
  )
  design <- make_manual_design(model = model, kk = 0)

  design <- validate_manual_design(design)

  expect_match(
    paste(design$errors, collapse = " "),
    "kk value must match the model covariate count"
  )
})

test_that("sdpd_design validates zero-covariate x is NULL", {
  design <- make_manual_design(
    x = matrix(
      1,
      nrow = 3,
      ncol = 3,
      dimnames = list(c("1", "2", "3"), c("t1", "t2", "t3"))
    )
  )

  design <- validate_manual_design(design)

  expect_match(
    paste(design$errors, collapse = " "),
    "zero-covariate design must store x as NULL"
  )
})

test_that("sdpd_design validates covariate x is present", {
  model <- build_sdpd_model(
    covariates = "x1",
    fixed_effects = FALSE
  )
  design <- make_manual_design(model = model, kk = 1, x = NULL)

  design <- validate_manual_design(design)

  expect_match(
    paste(design$errors, collapse = " "),
    "covariate design must include x"
  )
})

test_that("sdpd_design validates one-covariate x dimensions", {
  model <- build_sdpd_model(
    covariates = "x1",
    fixed_effects = FALSE
  )
  x <- matrix(
    1,
    nrow = 2,
    ncol = 3,
    dimnames = list(c("1", "2"), c("t1", "t2", "t3"))
  )
  design <- make_manual_design(model = model, kk = 1, x = x)

  design <- validate_manual_design(design)

  expect_match(
    paste(design$errors, collapse = " "),
    "one-covariate design x matrix must have dimensions pp by nn"
  )
})

test_that("sdpd_design validates multi-covariate x dimnames", {
  model <- build_sdpd_model(
    covariates = c("x1", "x2"),
    fixed_effects = FALSE
  )
  x <- array(
    1,
    dim = c(2, 3, 3),
    dimnames = list(
      c("x2", "x1"),
      c("a", "b", "c"),
      c("t1", "t2", "t3")
    )
  )
  design <- make_manual_design(model = model, kk = 2, x = x)

  design <- validate_manual_design(design)
  errors <- paste(design$errors, collapse = " ")

  expect_match(errors, "location names must match unit_index")
  expect_match(errors, "covariate names must match the model covariates")
})

test_that("sdpd_design validates mu and time_effects metadata", {
  design <- make_manual_design(
    mu = c("1" = 0, "2" = 0, "bad" = 0),
    time_effects = c("t1" = "0", "t2" = "0", "t3" = "0")
  )

  design <- validate_manual_design(design)
  errors <- paste(design$errors, collapse = " ")

  expect_match(errors, "mu names must match unit_index")
  expect_match(errors, "time_effects vector must be numeric")
})
