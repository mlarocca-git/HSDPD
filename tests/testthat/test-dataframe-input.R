test_that("read_data_from_dataframe accepts valid data-frame inputs", {
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
  px_neighbors <- list(index = ww_index)
  groups <- data.frame(
    COD = c(1, 1, 2),
    LABEL = c("north", "north", "south"),
    row.names = c("1", "2", "3")
  )
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE)

  result <- read_data_from_dataframe(
    px = c("1", "2"),
    rr_y = rr_y,
    rr_groups = groups,
    model = model,
    ww_index = ww_index,
    ww_values = ww_values,
    px_neighbors = px_neighbors
  )

  expect_null(result$error)
  expect_type(result, "list")
  expect_named(
    result,
    c(
      "series",
      "xx",
      "ww_index",
      "ww_values",
      "px_neighbors",
      "px",
      "lon",
      "lat",
      "group"
    )
  )
  expect_equal(result$px, c("1", "2"))
  expect_equal(names(result$lat), c("1", "2"))
  expect_equal(names(result$lon), c("1", "2"))
  expect_equal(rownames(result$series), c("1", "2", "3"))
  expect_equal(colnames(result$series), c("t1", "t2", "t3"))
  expect_equal(rownames(result$group), c("1", "2"))
})

test_that("read_data_from_dataframe accepts explicit coordinate vectors", {
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

  result <- read_data_from_dataframe(
    px = c("2", "3"),
    lat = c("1" = 40.1, "2" = 40.2, "3" = 40.3),
    lon = c("1" = 14.1, "2" = 14.2, "3" = 14.3),
    rr_y = rr_y,
    model = model,
    ww_index = ww_index,
    ww_values = ww_values,
    px_neighbors = list(index = ww_index)
  )

  expect_null(result$error)
  expect_equal(result$lat, c("2" = 40.2, "3" = 40.3))
  expect_equal(result$lon, c("2" = 14.2, "3" = 14.3))
})

test_that("read_data_from_dataframe reports invalid data-frame inputs", {
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE)
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
  valid_rr_y <- matrix(
    c(1, 2, 3,
      2, 3, 4,
      3, 4, 5),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("1", "2", "3"), c("t1", "t2", "t3"))
  )

  no_row_names <- unname(valid_rr_y)

  expect_match(
    read_data_from_dataframe(
      px = c("1", "2", "3"),
      rr_y = no_row_names,
      model = model,
      ww_index = ww_index,
      ww_values = ww_values,
      px_neighbors = list(index = ww_index)
    )$error,
    "Row names"
  )

  expect_match(
    read_data_from_dataframe(
      px = c("1", "4"),
      rr_y = valid_rr_y,
      model = model,
      ww_index = ww_index,
      ww_values = ww_values,
      px_neighbors = list(index = ww_index)
    )$error,
    "Some px values"
  )

  expect_match(
    read_data_from_dataframe(
      px = c("1", "2", "3"),
      lat = c(40.1, 40.2),
      lon = c(14.1, 14.2),
      rr_y = valid_rr_y,
      model = model,
      ww_index = ww_index,
      ww_values = ww_values,
      px_neighbors = list(index = ww_index)
    )$error,
    "Latitude and longitude"
  )

  expect_match(
    read_data_from_dataframe(
      px = c("1", "2", "3"),
      rr_y = valid_rr_y,
      rr_groups = data.frame(group = c(1, 1, 2)),
      model = model,
      ww_index = ww_index,
      ww_values = ww_values,
      px_neighbors = list(index = ww_index)
    )$error,
    "COD and LABEL"
  )
})
