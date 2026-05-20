make_test_spatraster <- function(values_start = 1, varname = "y") {
  rr <- terra::rast(
    nrows = 3,
    ncols = 3,
    nlyrs = 3,
    xmin = 0,
    xmax = 3,
    ymin = 0,
    ymax = 3
  )
  terra::values(rr) <- matrix(values_start:(values_start + 26), nrow = 9, ncol = 3)
  terra::time(rr) <- as.Date("2020-01-01") + 0:2
  terra::varnames(rr) <- varname
  rr
}

make_raster_vec_options <- function() {
  list(px_core = 1, px_neighbors = 0, na_rm = TRUE)
}

make_raster_vec_options_with_neighbors <- function() {
  list(px_core = 1, px_neighbors = 1, na_rm = TRUE)
}

test_that("read_data_from_raster returns expected top-level fields", {
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)

  result <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  expect_null(result$error)
  expect_named(
    result,
    c(
      "series",
      "xx",
      "ww_index",
      "ww_values",
      "px_neighbors",
      "na_summary",
      "px",
      "lon",
      "lat",
      "group"
    )
  )
})

test_that("read_data_from_raster returns stable no-covariate dimensions and names", {
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)

  result <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  expect_null(result$error)
  expect_null(result$xx)
  expect_equal(dim(result$series), c(9L, 3L))
  expect_equal(
    rownames(result$series),
    c("5", "1", "2", "3", "4", "6", "7", "8", "9")
  )
  expect_equal(
    colnames(result$series),
    as.character(as.Date("2020-01-01") + 0:2)
  )
  expect_equal(dim(result$ww_index), c(9L, 9L))
  expect_equal(rownames(result$ww_index), rownames(result$series))
  expect_equal(rownames(result$ww_values), rownames(result$series))
})

test_that("read_data_from_raster preserves selected px after current filtering", {
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)

  result <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  expect_null(result$error)
  expect_equal(result$px, 5)
  expect_equal(names(result$lat), "5")
  expect_equal(names(result$lon), "5")
  expect_true(as.character(result$px) %in% rownames(result$series))
})

test_that("read_data_from_raster returns na_summary", {
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)

  result <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  expect_null(result$error)
  expect_s3_class(result$na_summary, "data.frame")
  expect_equal(rownames(result$na_summary), rownames(result$series))
  expect_true(all(c("lat", "lon", "y", "y.XORGroups") %in% names(result$na_summary)))
  expect_true("isolated_points" %in% names(result$na_summary))
})

test_that("read_data_from_raster accepts one covariate SpatRaster", {
  rr_y <- make_test_spatraster()
  rr_x <- make_test_spatraster(values_start = 101, varname = "x1")
  model <- build_sdpd_model(
    covariates = "x1",
    fixed_effects = FALSE,
    check = FALSE
  )

  result <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    rr_xx = list(x1 = rr_x),
    model = model,
    vec_options = make_raster_vec_options()
  )

  expect_null(result$error)
  expect_true(is.array(result$xx))
  expect_equal(dim(result$xx), c(1L, 9L, 3L))
  expect_equal(dimnames(result$xx)[[1]], "x1")
  expect_equal(dimnames(result$xx)[[2]], rownames(result$series))
  expect_equal(dimnames(result$xx)[[3]], colnames(result$series))
})

test_that("read_data_from_raster returns current error list for invalid type_w", {
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)

  result <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options(),
    type_w = "bad"
  )

  expect_type(result, "list")
  expect_named(result, "error")
  expect_match(result$error, "type_w")
})

test_that("raster components bridge prepares no-covariate dataframe args", {
  raster_components_to_dataframe_args <- getFromNamespace(
    ".raster_components_to_dataframe_args",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  args <- raster_components_to_dataframe_args(components)

  expect_null(components$error)
  expect_true(is.matrix(args$rr_y))
  expect_equal(args$rr_y, components$series)
  expect_equal(args$px, as.character(components$px))
  expect_true(is.matrix(args$ww_index))
  expect_true(is.matrix(args$ww_values))
  expect_null(args$rr_xx)
  expect_true(is.matrix(args$rr_groups) || is.data.frame(args$rr_groups))
  expect_true(all(c("COD", "LABEL") %in% colnames(args$rr_groups)))
})

test_that("raster components bridge preserves one-covariate array metadata", {
  raster_components_to_dataframe_args <- getFromNamespace(
    ".raster_components_to_dataframe_args",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  rr_x <- make_test_spatraster(values_start = 101, varname = "x1")
  model <- build_sdpd_model(
    covariates = "x1",
    fixed_effects = FALSE,
    check = FALSE
  )
  components <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    rr_xx = list(x1 = rr_x),
    model = model,
    vec_options = make_raster_vec_options()
  )

  args <- raster_components_to_dataframe_args(components)

  expect_null(components$error)
  expect_true(is.array(args$rr_xx))
  expect_equal(args$rr_xx, components$xx)
  expect_equal(dimnames(args$rr_xx), dimnames(components$xx))
  expect_equal(dimnames(args$rr_xx)[[1]], "x1")
})

test_that("raster components bridge leaves na_summary raster-only", {
  raster_components_to_dataframe_args <- getFromNamespace(
    ".raster_components_to_dataframe_args",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  args <- raster_components_to_dataframe_args(components)

  expect_null(components$error)
  expect_s3_class(components$na_summary, "data.frame")
  expect_false("na_summary" %in% names(args))
})

test_that("raster components bridge preserves NULL px_neighbors", {
  raster_components_to_dataframe_args <- getFromNamespace(
    ".raster_components_to_dataframe_args",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  args <- raster_components_to_dataframe_args(components)

  expect_null(components$error)
  expect_null(components$px_neighbors)
  expect_null(args$px_neighbors)
})

test_that("build_sdpd_series forwards explicit px for SpatRaster input", {
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)

  result <- build_sdpd_series(
    df_obj = list(px = 1),
    px = 5,
    rr_y = rr_y,
    model = model,
    check = FALSE,
    vec_options = make_raster_vec_options()
  )

  expect_null(result$error)
  expect_equal(result$px, 5)
  expect_equal(names(result$lat), "5")
  expect_equal(names(result$lon), "5")
  expect_true("5" %in% rownames(result$series))
})

test_that("raster dataframe compatibility is TRUE for safe all-selected components", {
  raster_components_are_dataframe_compatible <- getFromNamespace(
    ".raster_components_are_dataframe_compatible",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- read_data_from_raster(
    px = seq_len(9),
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options_with_neighbors()
  )

  expect_null(components$error)
  expect_true(raster_components_are_dataframe_compatible(components))
})

test_that("raster dataframe compatibility is FALSE for subset px with NULL px_neighbors", {
  raster_components_are_dataframe_compatible <- getFromNamespace(
    ".raster_components_are_dataframe_compatible",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  expect_null(components$error)
  expect_null(components$px_neighbors)
  expect_false(raster_components_are_dataframe_compatible(components))
})

test_that("raster dataframe compatibility is TRUE for safe all-selected NULL px_neighbors", {
  raster_components_are_dataframe_compatible <- getFromNamespace(
    ".raster_components_are_dataframe_compatible",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- read_data_from_raster(
    px = seq_len(9),
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  expect_null(components$error)
  expect_null(components$px_neighbors)
  expect_true(raster_components_are_dataframe_compatible(components))
})

test_that("raster dataframe compatibility is FALSE for subset px group mismatch", {
  raster_components_are_dataframe_compatible <- getFromNamespace(
    ".raster_components_are_dataframe_compatible",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )
  components$px_neighbors <- list(index = components$ww_index)

  expect_null(components$error)
  expect_false(
    identical(rownames(components$group), rownames(components$series))
  )
  expect_false(raster_components_are_dataframe_compatible(components))
})

test_that("raster dataframe compatibility is FALSE for covariate components", {
  raster_components_are_dataframe_compatible <- getFromNamespace(
    ".raster_components_are_dataframe_compatible",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  rr_x <- make_test_spatraster(values_start = 101, varname = "x1")
  model <- build_sdpd_model(
    covariates = "x1",
    fixed_effects = FALSE,
    check = FALSE
  )
  components <- read_data_from_raster(
    px = seq_len(9),
    rr_y = rr_y,
    rr_xx = list(x1 = rr_x),
    model = model,
    vec_options = make_raster_vec_options_with_neighbors()
  )

  expect_null(components$error)
  expect_true(is.array(components$xx))
  expect_false(raster_components_are_dataframe_compatible(components))
})

test_that("raster dataframe helper preserves safe component behavior", {
  raster_components_via_dataframe <- getFromNamespace(
    ".raster_components_via_dataframe",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- read_data_from_raster(
    px = seq_len(9),
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options_with_neighbors()
  )

  result <- raster_components_via_dataframe(components, model)

  expect_null(components$error)
  expect_named(
    result,
    c(
      "series",
      "xx",
      "ww_index",
      "ww_values",
      "px_neighbors",
      "na_summary",
      "px",
      "lon",
      "lat",
      "group"
    )
  )
  expect_equal(result$series, components$series)
  expect_equal(result$group, components$group)
  expect_equal(result$ww_index, components$ww_index)
  expect_equal(result$ww_values, components$ww_values)
  expect_equal(result$na_summary, components$na_summary)
  expect_equal(result$px, components$px)
})

test_that("raster dataframe helper preserves one-covariate xx array", {
  raster_input_to_dataframe_components <- getFromNamespace(
    ".raster_input_to_dataframe_components",
    "HSDPD"
  )
  raster_components_via_dataframe <- getFromNamespace(
    ".raster_components_via_dataframe",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  rr_x <- make_test_spatraster(values_start = 101, varname = "x1")
  model <- build_sdpd_model(
    covariates = "x1",
    fixed_effects = FALSE,
    check = FALSE
  )
  components <- raster_input_to_dataframe_components(
    px = seq_len(9),
    rr_y = rr_y,
    rr_xx = list(x1 = rr_x),
    model = model,
    vec_options = make_raster_vec_options()
  )

  result <- raster_components_via_dataframe(components, model)

  expect_null(components$error)
  expect_true(is.array(result$xx))
  expect_equal(dim(result$xx), dim(components$xx))
  expect_equal(dimnames(result$xx), dimnames(components$xx))
  expect_equal(result$xx, components$xx)
})

test_that("read_data_from_raster routes safe all-selected components through dataframe path", {
  raster_input_to_dataframe_components <- getFromNamespace(
    ".raster_input_to_dataframe_components",
    "HSDPD"
  )
  raster_components_via_dataframe <- getFromNamespace(
    ".raster_components_via_dataframe",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- raster_input_to_dataframe_components(
    px = seq_len(9),
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options_with_neighbors()
  )
  expected <- raster_components_via_dataframe(components, model)

  result <- read_data_from_raster(
    px = seq_len(9),
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options_with_neighbors()
  )

  expect_null(result$error)
  expect_named(result, names(expected))
  expect_equal(result$series, expected$series)
  expect_equal(result$group, expected$group)
  expect_equal(result$ww_index, expected$ww_index)
  expect_equal(result$ww_values, expected$ww_values)
  expect_equal(result$na_summary, components$na_summary)
  expect_equal(result$px, components$px)
})

test_that("read_data_from_raster keeps NULL px_neighbors case on legacy path", {
  raster_input_to_dataframe_components <- getFromNamespace(
    ".raster_input_to_dataframe_components",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  expected <- raster_input_to_dataframe_components(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  result <- read_data_from_raster(
    px = 5,
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  expect_null(result$error)
  expect_null(result$px_neighbors)
  expect_equal(result, expected)
})

test_that("read_data_from_raster routes safe all-selected NULL px_neighbors", {
  raster_input_to_dataframe_components <- getFromNamespace(
    ".raster_input_to_dataframe_components",
    "HSDPD"
  )
  raster_components_via_dataframe <- getFromNamespace(
    ".raster_components_via_dataframe",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE, check = FALSE)
  components <- raster_input_to_dataframe_components(
    px = seq_len(9),
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )
  expected <- raster_components_via_dataframe(components, model)

  result <- read_data_from_raster(
    px = seq_len(9),
    rr_y = rr_y,
    model = model,
    vec_options = make_raster_vec_options()
  )

  expect_null(result$error)
  expect_named(result, names(expected))
  expect_null(result$px_neighbors)
  expect_equal(result$series, expected$series)
  expect_equal(result$ww_index, expected$ww_index)
  expect_equal(result$ww_values, expected$ww_values)
  expect_equal(result$na_summary, components$na_summary)
  expect_equal(result$px, components$px)
})

test_that("read_data_from_raster keeps covariate raster case on legacy path", {
  raster_input_to_dataframe_components <- getFromNamespace(
    ".raster_input_to_dataframe_components",
    "HSDPD"
  )
  rr_y <- make_test_spatraster()
  rr_x <- make_test_spatraster(values_start = 101, varname = "x1")
  model <- build_sdpd_model(
    covariates = "x1",
    fixed_effects = FALSE,
    check = FALSE
  )
  expected <- raster_input_to_dataframe_components(
    px = seq_len(9),
    rr_y = rr_y,
    rr_xx = list(x1 = rr_x),
    model = model,
    vec_options = make_raster_vec_options_with_neighbors()
  )

  result <- read_data_from_raster(
    px = seq_len(9),
    rr_y = rr_y,
    rr_xx = list(x1 = rr_x),
    model = model,
    vec_options = make_raster_vec_options_with_neighbors()
  )

  expect_null(result$error)
  expect_true(is.array(result$xx))
  expect_equal(result, expected)
})
