#' Read SDP-D Data from Raster Objects
#'
#' Extracts and organizes SDP-D input data from raster objects.
#'
#' The function builds a spatio-temporal series from raster data, optionally
#' extracts exogenous raster regressors, builds spatial-weight matrices, removes
#' missing or isolated pixels, and constructs the neighbor-boundary structure
#' used by local SDP-D estimation.
#'
#' @param px Optional vector. Raster cell identifiers to include. If `NULL` and
#'   coordinates are not supplied, all cells are used.
#' @param lat Optional numeric vector. Latitude coordinates used to select cells
#'   when `px` is `NULL`.
#' @param lon Optional numeric vector. Longitude coordinates used to select cells
#'   when `px` is `NULL`.
#' @param rr_y `SpatRaster`. Endogenous raster time series.
#' @param rr_xx Optional `SpatRaster` or list of `SpatRaster` objects.
#'   Exogenous raster regressors.
#' @param rr_groups Optional `SpatRaster`. Raster object defining group
#'   membership.
#' @param label_groups Optional named vector. Group labels indexed by group code.
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()].
#' @param vec_options List. Vectorization options. Expected elements include
#'   `px_core`, `px_neighbors`, and `na_rm`.
#' @param type_w Character scalar. Type of spatial weights. Currently only
#'   `"distance"` is supported.
#'
#' @return A list containing:
#' \describe{
#'   \item{series}{Endogenous series matrix.}
#'   \item{xx}{Optional regressor array.}
#'   \item{ww_index}{Spatial-neighbor index matrix.}
#'   \item{ww_values}{Spatial-weight value matrix.}
#'   \item{px_neighbors}{Pixel-neighbor structure with `index` and
#'   `series_boundary`.}
#'   \item{na_summary}{Data frame summarizing missing values and isolated
#'   pixels.}
#'   \item{px}{Selected pixel identifiers.}
#'   \item{lon}{Longitude vector.}
#'   \item{lat}{Latitude vector.}
#'   \item{group}{Group data frame with columns `COD` and `LABEL`.}
#' }
#'
#' If input validation fails, a list with an `error` element is returned.
#'
#' @details
#' If `px = "all"`, all raster cells are considered and cells with invalid data
#' are removed according to `vec_options$na_rm`.
#'
#' If `px` or coordinates are supplied, the function extracts only the selected
#' pixels and adds the required core and boundary neighbors.
#'
#' The returned object follows the renamed package API: `xx`, `ww_index`,
#' `ww_values`, `px_neighbors`, and `series_boundary`.
#'
#' @seealso [build_sdpd_series()], [read_data_from_dataframe()]
#'
#' @export
read_data_from_raster <- function(px = NULL,
                                  lat = NULL,
                                  lon = NULL,
                                  rr_y,
                                  rr_xx = NULL,
                                  rr_groups = NULL,
                                  label_groups = NULL,
                                  model,
                                  vec_options,
                                  type_w = "distance") {
  components <- .raster_input_to_dataframe_components(
    px = px,
    lat = lat,
    lon = lon,
    rr_y = rr_y,
    rr_xx = rr_xx,
    rr_groups = rr_groups,
    label_groups = label_groups,
    model = model,
    vec_options = vec_options,
    type_w = type_w
  )

  if (!is.null(components$error)) {
    return(components)
  }

  if (.raster_components_are_dataframe_compatible(components)) {
    return(.raster_components_via_dataframe(components, model))
  }

  components
}

.raster_input_to_dataframe_components <- function(px = NULL,
                                                  lat = NULL,
                                                  lon = NULL,
                                                  rr_y,
                                                  rr_xx = NULL,
                                                  rr_groups = NULL,
                                                  label_groups = NULL,
                                                  model,
                                                  vec_options,
                                                  type_w = "distance") {
  # Extract and build a spatio-temporal series from raster variables.
  # The px argument is kept first so the function can be used in parallelized
  # workflows.

  if (is.null(px)) {
    if (is.null(lat) || is.null(lon)) {
      px <- "all"
    } else {
      if (length(lat) == length(lon)) {
        n_px <- length(lat)
      } else {
        return(list(error = "Latitude and longitude have different lengths."))
      }

      px <- terra::cellFromXY(rr_y, cbind(lon, lat))
    }
  }

  if (is.numeric(px)) {
    n_px <- length(px)

    coordinates_temp <- terra::xyFromCell(rr_y, px)
    lat <- coordinates_temp[, "y"]
    lon <- coordinates_temp[, "x"]

    names(lat) <- px
    names(lon) <- px
  } else if (is.character(px)) {
    n_px <- length(px)

    if (px[1] == "all") {
      px <- seq_len(dim(terra::values(rr_y))[1])
    } else {
      px <- as.numeric(px)
    }

    coordinates_temp <- terra::xyFromCell(rr_y, px)
    lat <- coordinates_temp[, "y"]
    lon <- coordinates_temp[, "x"]

    names(lat) <- px
    names(lon) <- px
  }

  if (n_px == 0) {
    return(list(error = "There are no series to analyze."))
  }

  # Initialize series including close neighbors.
  if (vec_options$px_core == 0) {
    close_neighbors <- NULL
    ww_index <- NULL
    ww_values <- NULL
    px_neighbors <- NULL
  } else {
    core_matrix <- matrix(
      1,
      ncol = 2 * vec_options$px_core + 1,
      nrow = 2 * vec_options$px_core + 1
    )

    core_matrix[
      vec_options$px_core + 1,
      vec_options$px_core + 1
    ] <- 0

    close_neighbors <- terra::adjacent(
      rr_y,
      px,
      pairs = TRUE,
      directions = core_matrix,
      include = TRUE,
      symmetrical = FALSE
    )
  }

  indices <- stats::na.exclude(unique(c(px, close_neighbors)))
  time_index <- terra::time(rr_y)

  series <- terra::values(rr_y)[indices, ]
  dimnames(series)[[1]] <- indices
  dimnames(series)[[2]] <- as.character(time_index)

  tt <- length(time_index)
  pp <- dim(series)[1]

  # Create regressors.
  coordinates <- terra::xyFromCell(rr_y, indices)

  if (is.null(rr_xx)) {
    kk <- 0
    n_regressors <- 0
    xx <- NULL
    regressor_names <- NULL

    if (sum(model$beta_coeffs) > 0) {
      return(list(error = "The model includes exogenous regressors, but rr_xx was not supplied."))
    }
  } else if (!is.list(rr_xx)) {
    kk <- 1
    n_regressors <- 1

    if (sum(model$beta_coeffs) > 1) {
      return(list(error = "The model includes more than one exogenous regressor, but rr_xx does not contain enough variables."))
    }

    rr_xx <- list(rr_xx)
  } else {
    kk <- length(rr_xx)

    if (sum(model$beta_coeffs) > kk) {
      return(list(error = paste(
        "The model includes more than",
        kk,
        "exogenous regressors, but rr_xx does not contain enough variables."
      )))
    }

    if (!is.null(names(model$beta_coeffs))) {
      required_regressors <- names(model$beta_coeffs)[model$beta_coeffs]

      if (sum(required_regressors %in% names(rr_xx)) < sum(model$beta_coeffs)) {
        return(list(error = "The model includes exogenous regressors that are not contained in rr_xx."))
      }
    }
  }

  if (kk > 0) {
    n_regressors <- 0
    regressor_names <- character(kk)
    xx <- array(0, dim = c(kk, pp, tt))

    for (regressor_index in seq_len(dim(xx)[1])) {
      if (!is.null(rr_xx[[regressor_index]])) {
        n_regressors <- n_regressors + 1

        if (is.character(rr_xx[[regressor_index]])) {
          if (rr_xx[[regressor_index]] == "trend") {
            xx[n_regressors, , ] <- matrix(
              rep(seq_len(tt) / tt, pp),
              byrow = TRUE,
              nrow = pp
            )

            regressor_names[n_regressors] <- "trend"
          } else {
            n_regressors <- n_regressors - 1
          }
        } else if (!is.null(names(model$beta_coeffs))) {
          if (names(rr_xx)[regressor_index] %in% names(model$beta_coeffs)) {
            regressor_names[n_regressors] <- names(rr_xx)[regressor_index]

            p_xx <- terra::cellFromXY(rr_xx[[regressor_index]], coordinates)
            xx[n_regressors, , ] <- terra::values(rr_xx[[regressor_index]])[p_xx, ]
          } else {
            n_regressors <- n_regressors - 1
          }
        } else {
          regressor_names[n_regressors] <- names(rr_xx)[regressor_index]

          p_xx <- terra::cellFromXY(rr_xx[[regressor_index]], coordinates)
          xx[n_regressors, , ] <- terra::values(rr_xx[[regressor_index]])[p_xx, ]
        }
      }
    }

    if (n_regressors > 0) {
      kk <- n_regressors
      xx <- xx[seq_len(kk), , , drop = FALSE]
      regressor_names <- regressor_names[seq_len(kk)]

      dimnames(xx) <- list(
        regressor_names,
        indices,
        as.character(time_index)
      )
    }
  }

  # Define groups.
  if (is.null(rr_groups)) {
    groups <- rep(1, pp)
    group_labels <- rep("Common group", pp)

    groups <- cbind(COD = groups, LABEL = group_labels)
  } else {
    groups <- terra::values(rr_groups)[indices, 1]
    group_labels <- label_groups[as.character(groups)]

    groups <- data.frame(
      COD = groups,
      LABEL = group_labels
    )
  }

  dimnames(groups)[[1]] <- indices

  # Remove pixels with missing values.
  n_na_y <- apply(series, 1, FUN = function(x) sum(is.na(x)))
  n_na_group <- is.na(groups[, 1])
  n_na_group_xor <- xor(n_na_y > 0, n_na_group > 0)

  if (kk > 1) {
    n_na_x <- t(apply(xx, c(1, 2), FUN = function(x) sum(is.na(x))))
    n_na_xor <- apply(n_na_x, 2, FUN = function(x, y) xor(x > 0, y > 0), y = n_na_y)
    total_na_x <- apply(n_na_x, 1, sum)
  } else if (kk == 1) {
    n_na_x <- total_na_x <- apply(xx, 1, FUN = function(x) sum(is.na(x)))
    n_na_xor <- xor(n_na_x > 0, n_na_y > 0)
  } else {
    n_na_x <- n_na_xor <- total_na_x <- rep(0, pp)
  }

  if (vec_options$na_rm) {
    keep_indices <- n_na_y == 0 & n_na_group == 0 & total_na_x == 0
  } else {
    keep_indices <- n_na_y < tt & n_na_group == 0 & total_na_x < tt
  }

  na_indices <- indices[!keep_indices]

  if (sum(keep_indices) == 0) {
    return(list(error = "There are no series to analyze."))
  }

  y_name <- terra::varnames(rr_y)[1]

  if (kk > 0) {
    na_summary <- cbind(
      coordinates,
      n_na_y,
      n_na_x,
      n_na_xor,
      n_na_group_xor
    )

    dimnames(na_summary) <- list(
      indices,
      c(
        "lat",
        "lon",
        y_name,
        regressor_names,
        paste(y_name, "XOR", regressor_names, sep = ""),
        paste(y_name, "XORGroups")
      )
    )
  } else {
    na_summary <- cbind(
      coordinates,
      n_na_y,
      n_na_group_xor
    )

    dimnames(na_summary) <- list(
      indices,
      c("lat", "lon", y_name, paste(y_name, "XORGroups"))
    )
  }

  # Create spatial-weight vectors and check isolated points.
  if (vec_options$px_core > 0) {
    coordinates_1 <- terra::xyFromCell(rr_y, close_neighbors[, 1])
    coordinates_2 <- terra::xyFromCell(rr_y, close_neighbors[, 2])

    ww_index <- ww_values <- matrix(
      0,
      nrow = sum(keep_indices),
      ncol = (2 * vec_options$px_core + 1)^2
    )

    dimnames(ww_index)[[1]] <- dimnames(ww_values)[[1]] <- indices[keep_indices]

    if (type_w == "distance") {
      ww <- terra::distance(
        x = coordinates_1,
        y = coordinates_2,
        lonlat = TRUE,
        pairwise = TRUE
      )

      ww <- ifelse(ww < 0.01, 0, 1 / ww)
    } else {
      return(list(error = "The value set for type_w is not allowed."))
    }

    isolated_points <- logical(pp)
    names(isolated_points) <- indices

    for (ii in seq_len(sum(keep_indices))) {
      pixel <- indices[keep_indices][ii]

      # First check the first column of close_neighbors.
      ww_1 <- close_neighbors[, 1] == pixel
      ww_2 <- close_neighbors[ww_1, 2]

      to_exclude <- ww_2 %in% na_indices
      ww_2 <- ww_2[!to_exclude]

      if (sum(ww_1) - sum(to_exclude) > 0) {
        ww_index[ii, seq_along(ww_2)] <- ww_2
        ww_values[ii, seq_along(ww_2)] <- ww[ww_1][!to_exclude]
      } else {
        # Also check additional points in the second column of close_neighbors.
        ww_1 <- close_neighbors[, 2] == pixel
        ww_2 <- close_neighbors[ww_1, 1]

        to_exclude <- ww_2 %in% na_indices
        ww_2 <- ww_2[!to_exclude]

        if (sum(ww_1) - sum(to_exclude) > 0) {
          ww_index[ii, seq_along(ww_2)] <- ww_2
          ww_values[ii, seq_along(ww_2)] <- ww[ww_1][!to_exclude]
        }
      }

      if (sum(abs(ww_values[ii, ])) > 0) {
        ww_values[ii, ] <- ww_values[ii, ] / sum(abs(ww_values[ii, ]))
      } else {
        isolated_points[keep_indices][ii] <- TRUE
      }
    }

    ww_values <- ww_values[!(isolated_points[keep_indices]), , drop = FALSE]
    ww_index <- ww_index[!(isolated_points[keep_indices]), , drop = FALSE]

    keep_indices <- keep_indices & !isolated_points
    na_summary <- cbind(na_summary, isolated_points = isolated_points)
  }

  # Redefine objects after removing isolated points and missing values.
  indices <- as.character(indices[keep_indices])

  series <- series[indices, , drop = FALSE]

  if (!is.null(xx)) {
    xx <- xx[, indices, , drop = FALSE]
  }

  selected_px <- as.character(px) %in% indices

  lat <- lat[selected_px]
  lon <- lon[selected_px]
  px <- px[selected_px]

  n_px <- length(px)

  if (n_px == 0) {
    return(list(error = "All selected series have NA values in Y, covariates, or neighborhoods."))
  }

  if (!is.null(groups)) {
    groups <- groups[as.character(px), , drop = FALSE]
  }

  # Create index matrix for far-neighbor series.
  if (vec_options$px_neighbors > 0) {
    neighborhood_matrix <- matrix(
      1,
      ncol = 2 * vec_options$px_neighbors + 1,
      nrow = 2 * vec_options$px_neighbors + 1
    )

    neighborhood_matrix[
      vec_options$px_neighbors + 1,
      vec_options$px_neighbors + 1
    ] <- 0

    px_neighbors_index <- matrix(
      NA,
      nrow = length(indices),
      ncol = (2 * vec_options$px_neighbors + 1)^2 - 1
    )

    dimnames(px_neighbors_index)[[1]] <- indices

    for (ii in seq_len(n_px)) {
      temp_neighbors <- terra::adjacent(
        rr_y,
        cells = px[ii],
        directions = neighborhood_matrix
      )

      temp_neighbors <- temp_neighbors[temp_neighbors > 0]

      px_neighbors_index[
        as.character(px[ii]),
        seq_along(temp_neighbors)
      ] <- temp_neighbors
    }

    neighbor_indices <- stats::na.exclude(unique(as.vector(px_neighbors_index)))
    remaining_indices <- indices[!(as.numeric(indices) %in% px)]

    for (ii in remaining_indices) {
      temp_neighbors <- terra::adjacent(
        rr_y,
        cells = as.numeric(ii),
        directions = neighborhood_matrix
      )

      temp_neighbors <- temp_neighbors[temp_neighbors > 0]
      temp_neighbors <- temp_neighbors[temp_neighbors %in% neighbor_indices]

      px_neighbors_index[
        as.character(ii),
        seq_along(temp_neighbors)
      ] <- temp_neighbors
    }

    # Remove far-neighbor series with NA values.
    neighbor_indices <- stats::na.exclude(unique(as.vector(px_neighbors_index)))

    neighbor_series <- terra::values(rr_y)[neighbor_indices, ]
    dimnames(neighbor_series)[[1]] <- neighbor_indices
    dimnames(neighbor_series)[[2]] <- as.character(time_index)

    neighbor_na <- apply(neighbor_series, 1, FUN = function(x) sum(is.na(x)))

    px_neighbors_index <- t(apply(
      px_neighbors_index,
      1,
      FUN = function(x, ind) {
        ifelse(x %in% ind, NA, x)
      },
      ind = neighbor_indices[neighbor_na > 0]
    ))

    # Extract far-neighbor boundary series for covariance computation.
    neighbor_indices <- stats::na.exclude(unique(as.vector(px_neighbors_index)))

    external_indices <- neighbor_indices[
      !(neighbor_indices %in% indices)
    ]

    neighbor_series <- neighbor_series[as.character(external_indices), , drop = FALSE]

    px_neighbors <- list(
      index = px_neighbors_index,
      series_boundary = neighbor_series
    )
  } else {
    px_neighbors <- NULL
  }

  # Output.
  list(
    series = series,
    xx = xx,
    ww_index = ww_index,
    ww_values = ww_values,
    px_neighbors = px_neighbors,
    na_summary = data.frame(na_summary),
    px = px,
    lon = lon,
    lat = lat,
    group = groups
  )
}

.raster_components_to_dataframe_args <- function(components) {
  list(
    px = as.character(components$px),
    lat = components$lat,
    lon = components$lon,
    rr_y = components$series,
    rr_xx = components$xx,
    rr_groups = components$group,
    ww_index = components$ww_index,
    ww_values = components$ww_values,
    px_neighbors = components$px_neighbors
  )
}

.raster_components_are_dataframe_compatible <- function(components) {
  if (!is.null(components$error) ||
      !is.matrix(components$series) ||
      !is.matrix(components$ww_index) ||
      !is.matrix(components$ww_values) ||
      !identical(as.character(components$px), rownames(components$series)) ||
      is.null(components$group) ||
      !identical(rownames(components$group), rownames(components$series))) {
    return(FALSE)
  }

  if (!is.null(components$xx)) {
    return(
      is.array(components$xx) &&
        length(dim(components$xx)) == 3 &&
        identical(dimnames(components$xx)[[2]], rownames(components$series)) &&
        identical(dimnames(components$xx)[[3]], colnames(components$series))
    )
  }

  TRUE
}

.raster_components_via_dataframe <- function(components, model) {
  args <- .raster_components_to_dataframe_args(components)
  result <- do.call(read_data_from_dataframe, c(args, list(model = model)))

  if (!is.null(result$error)) {
    return(components)
  }

  list(
    series = result$series,
    xx = result$xx,
    ww_index = result$ww_index,
    ww_values = result$ww_values,
    px_neighbors = result$px_neighbors,
    na_summary = components$na_summary,
    px = components$px,
    lon = components$lon,
    lat = components$lat,
    group = result$group
  )
}
