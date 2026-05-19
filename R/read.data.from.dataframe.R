.validate_sdpd_dataframe_input <- function(px,
                                           lat = NULL,
                                           lon = NULL,
                                           rr_y,
                                           rr_xx = NULL,
                                           rr_groups = NULL) {
  if (!is.data.frame(rr_y) && !is.matrix(rr_y)) {
    return(list(error = "The rr_y argument must be a data frame or matrix."))
  }

  rr_y <- as.matrix(rr_y)

  if (is.null(dimnames(rr_y)[[1]])) {
    return(list(error = "Row names, i.e. locations, must be defined in rr_y."))
  }

  if (is.null(px)) {
    px <- dimnames(rr_y)[[1]]
  }

  px <- as.character(px)
  indices_all <- dimnames(rr_y)[[1]]

  if (sum(px %in% indices_all) < length(px)) {
    return(list(error = "Some px values are not contained in the rr_y data frame."))
  }

  if (is.null(dimnames(rr_y)[[2]])) {
    latitude_column <- numeric(0)
    longitude_column <- numeric(0)
  } else {
    latitude_column <- which(dimnames(rr_y)[[2]] == "latitude")
    longitude_column <- which(dimnames(rr_y)[[2]] == "longitude")
  }

  coordinates_from_columns <- FALSE

  if (length(latitude_column) != 0 && length(longitude_column) != 0) {
    lat <- rr_y[px, latitude_column]
    lon <- rr_y[px, longitude_column]

    names(lat) <- px
    names(lon) <- px

    rr_y <- rr_y[, -c(latitude_column, longitude_column), drop = FALSE]
    coordinates_from_columns <- TRUE
  } else if (length(latitude_column) != 0 || length(longitude_column) != 0) {
    rr_y <- rr_y[, -c(latitude_column, longitude_column), drop = FALSE]
  }

  if (is.null(lat) || is.null(lon)) {
    lat <- NULL
    lon <- NULL
  } else if (coordinates_from_columns) {
    lat <- as.numeric(lat)
    lon <- as.numeric(lon)

    names(lat) <- px
    names(lon) <- px
  } else if (is.numeric(lat) && is.numeric(lon)) {
    if (length(lat) == 1 && length(lon) == 1) {
      temp_lat <- rr_y[px, lat]
      temp_lon <- rr_y[px, lon]

      names(temp_lat) <- px
      names(temp_lon) <- px

      rr_y <- rr_y[, -c(lat, lon), drop = FALSE]

      lat <- temp_lat
      lon <- temp_lon
    } else if (length(lat) == dim(rr_y)[1] &&
               length(lon) == dim(rr_y)[1]) {
      if (is.null(names(lat))) {
        names(lat) <- dimnames(rr_y)[[1]]
      }

      if (is.null(names(lon))) {
        names(lon) <- dimnames(rr_y)[[1]]
      }

      lat <- lat[px]
      lon <- lon[px]
    } else {
      return(list(error = "Latitude and longitude values do not have the expected format or length."))
    }
  } else {
    return(list(error = "Latitude and longitude values do not have the expected format or length."))
  }

  if (!is.null(rr_groups)) {
    if (!is.data.frame(rr_groups) && !is.matrix(rr_groups)) {
      return(list(error = "The rr_groups argument must be a data frame or matrix."))
    }

    if (dim(rr_groups)[1] != dim(rr_y)[1]) {
      return(list(error = "The rr_groups object must have the same number of rows as rr_y."))
    }

    if (sum(c("COD", "LABEL") %in% dimnames(rr_groups)[[2]]) < 2) {
      return(list(error = "The rr_groups object must contain the columns COD and LABEL."))
    }
  }

  .new_sdpd_data(
    series = rr_y,
    px = px,
    lat = lat,
    lon = lon,
    rr_xx = rr_xx,
    rr_groups = rr_groups
  )
}

#' Read SDP-D Data from a Data Frame or Matrix
#'
#' Extracts and organizes SDP-D input data from a data frame or matrix.
#'
#' The function reads the endogenous series, optional coordinates, groups,
#' spatial-weight components, pixel-neighbor information, and optional
#' regressors. It also builds the boundary series needed for local spatial
#' estimation.
#'
#' @param px Optional vector. Spatial-unit identifiers to include. If `NULL`,
#'   all row names of `rr_y` are used.
#' @param lat Optional numeric vector or column index identifying latitudes.
#' @param lon Optional numeric vector or column index identifying longitudes.
#' @param rr_y Data frame or matrix. Endogenous series. Rows must identify
#'   locations and columns must identify time points, except optional coordinate
#'   columns named `latitude` and `longitude`.
#' @param rr_xx Optional list or array. Exogenous regressors.
#' @param rr_groups Optional data frame or matrix. Group information. Must
#'   contain columns `COD` and `LABEL`.
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()].
#' @param ww_index Matrix. Spatial-neighbor index matrix.
#' @param ww_values Matrix. Spatial-weight value matrix.
#' @param px_neighbors List. Pixel or proximity-neighbor structure. Expected
#'   elements are `index` and optionally `series_boundary`.
#'
#' @return A list containing:
#' \describe{
#'   \item{series}{Endogenous series for selected pixels and required close
#'   neighbors.}
#'   \item{xx}{Optional regressor array.}
#'   \item{ww_index}{Spatial-neighbor index matrix restricted to required
#'   locations.}
#'   \item{ww_values}{Spatial-weight value matrix restricted to required
#'   locations.}
#'   \item{px_neighbors}{Pixel-neighbor structure with `index` and
#'   `series_boundary`.}
#'   \item{px}{Selected spatial-unit identifiers.}
#'   \item{lon}{Longitude vector, if available.}
#'   \item{lat}{Latitude vector, if available.}
#'   \item{group}{Group matrix with columns `COD` and `LABEL`.}
#' }
#'
#' If input validation fails, a list with an `error` element is returned.
#'
#' @details
#' The function adds close spatial neighbors required by the spatial weight
#' matrix and builds a boundary series for farther neighbors found in
#' `px_neighbors$index`.
#'
#' For consistency with the renamed package API, the returned regressor component
#' is named `xx`, not `X`, and the boundary component is named
#' `series_boundary`, not `seriesBoundary`.
#'
#' @examples
#' rr_y <- matrix(
#'   c(1, 2, 3,
#'     2, 3, 4,
#'     3, 4, 5),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#' rownames(rr_y) <- colnames(rr_y) <- c("1", "2", "3")
#'
#' ww_index <- matrix(
#'   c(2, 3,
#'     1, 3,
#'     1, 2),
#'   nrow = 3,
#'   byrow = TRUE,
#'   dimnames = list(c("1", "2", "3"), NULL)
#' )
#' ww_values <- matrix(
#'   0.5,
#'   nrow = 3,
#'   ncol = 2,
#'   dimnames = list(c("1", "2", "3"), NULL)
#' )
#' px_neighbors <- list(index = ww_index)
#' model <- build_sdpd_model(covariates = 0, fixed_effects = FALSE)
#'
#' series_object <- read_data_from_dataframe(
#'   px = c("1", "2", "3"),
#'   rr_y = rr_y,
#'   model = model,
#'   ww_index = ww_index,
#'   ww_values = ww_values,
#'   px_neighbors = px_neighbors
#' )
#'
#' names(series_object)
#' @seealso [build_sdpd_series()]
#'
#' @export
read_data_from_dataframe <- function(px,
                                     lat = NULL,
                                     lon = NULL,
                                     rr_y,
                                     rr_xx = NULL,
                                     rr_groups = NULL,
                                     model,
                                     ww_index,
                                     ww_values,
                                     px_neighbors) {
  # Extract and organize matrix/data-frame input data.
  # The function returns the spatio-temporal series, close-neighbor information,
  # spatial weights, and optional regressors.

  input_data <- .validate_sdpd_dataframe_input(
    px = px,
    lat = lat,
    lon = lon,
    rr_y = rr_y,
    rr_xx = rr_xx,
    rr_groups = rr_groups
  )

  if (!is.null(input_data$error)) {
    return(input_data)
  }

  rr_y <- input_data$series
  px <- input_data$px
  lat <- input_data$lat
  lon <- input_data$lon
  rr_xx <- input_data$rr_xx
  rr_groups <- input_data$rr_groups

  # Derive close-neighbor series.
  if (is.null(ww_index) || is.null(ww_values)) {
    return(list(error = "The spatial matrix cannot be built without ww_index and ww_values."))
  }

  if (is.matrix(ww_index)) {
    close_neighbors <- as.numeric(ww_index[as.character(px), ])
    close_neighbors <- close_neighbors[close_neighbors > 0]
    close_neighbors <- unique(close_neighbors)

    indices <- as.character(unique(c(px, close_neighbors)))

    if (length(indices) > 1) {
      series <- as.matrix(rr_y[indices, ])
    } else if (length(indices) == 1) {
      series <- matrix(rr_y[indices, ], nrow = 1)
      dimnames(series)[[1]] <- indices
    } else {
      return(list(error = "There are no data in the series matrix."))
    }

    time_index <- dimnames(series)[[2]]
    tt <- dim(series)[2]
    pp <- length(indices)
  } else {
    return(list(error = "The ww_index argument must be a matrix."))
  }

  # Derive groups.
  if (is.null(rr_groups) && is.null(model$groups)) {
    groups <- rep(1, length(px))
    group_labels <- rep("group_1", length(px))

    groups <- cbind(COD = groups, LABEL = group_labels)
    dimnames(groups)[[1]] <- px
  } else if (!is.null(rr_groups)) {
    if (!is.data.frame(rr_groups) && !is.matrix(rr_groups)) {
      return(list(error = "The rr_groups argument must be a data frame or matrix."))
    }

    if (dim(rr_groups)[1] != dim(rr_y)[1]) {
      return(list(error = "The rr_groups object must have the same number of rows as rr_y."))
    }

    if (sum(c("COD", "LABEL") %in% dimnames(rr_groups)[[2]]) < 2) {
      return(list(error = "The rr_groups object must contain the columns COD and LABEL."))
    }

    groups <- as.matrix(rr_groups[as.character(px), c("COD", "LABEL")])
  } else if (!is.null(model$groups)) {
    groups <- as.matrix(model$groups[indices, ])
  }

  # Derive spatial-weight and far-neighbor matrices.
  if (is.matrix(ww_values)) {
    ww_values <- ww_values[indices, ]
    ww_index <- ww_index[indices, ]
  } else {
    return(list(error = "The ww_values argument must be a matrix."))
  }

  # Correct spatial-matrix and neighbor-matrix indices and weights by adding
  # boundary effects.
  if (is.null(px_neighbors)) {
    return(list(error = "The px_neighbors matrix is required."))
  }

  neighbor_indices_around <- unique(stats::na.omit(as.vector(
    px_neighbors$index[as.character(px), ]
  )))

  neighbors <- as.matrix(px_neighbors$index[indices, ])

  remaining_indices <- indices[!(indices %in% px)]

  for (row_id in remaining_indices) {
    present_in_px <- ww_index[row_id, ] %in% px

    ww_index[row_id, !present_in_px] <- 0
    ww_values[row_id, !present_in_px] <- 0

    ww_values_sum <- sum(abs(ww_values[row_id, ]))

    if (ww_values_sum > 0) {
      ww_values[row_id, ] <- ww_values[row_id, ] / ww_values_sum
    }

    present_in_neighbors <- neighbors[row_id, ] %in% as.numeric(neighbor_indices_around)
    neighbors[row_id, !present_in_neighbors] <- NA
  }

  # Derive far-neighbor boundary series.
  neighbor_indices_around <- as.character(unique(stats::na.omit(as.vector(neighbors))))

  external_indices <- neighbor_indices_around[
    !(neighbor_indices_around %in% indices)
  ]

  residual_indices <- external_indices[
    !(external_indices %in% dimnames(rr_y)[[1]])
  ]

  external_indices <- external_indices[
    external_indices %in% dimnames(rr_y)[[1]]
  ]

  if (length(external_indices) > 1) {
    series_boundary_first <- as.matrix(rr_y[as.character(external_indices), ])
  } else if (length(external_indices) == 1) {
    series_boundary_first <- matrix(
      rr_y[as.character(external_indices), ],
      nrow = 1
    )

    dimnames(series_boundary_first)[[1]] <- external_indices
  } else {
    series_boundary_first <- NULL
  }

  series_boundary_second <- NULL

  if (length(residual_indices) > 0 &&
      is.null(px_neighbors$series_boundary)) {
    temp_dimnames <- dimnames(neighbors)

    neighbors <- t(apply(
      neighbors,
      1,
      FUN = function(x, ind) {
        ifelse(x %in% ind, NA, x)
      },
      ind = as.numeric(residual_indices)
    ))

    dimnames(neighbors) <- temp_dimnames
  } else if (length(residual_indices) > 1) {
    series_boundary_second <- as.matrix(
      px_neighbors$series_boundary[as.character(residual_indices), ]
    )
  } else if (length(residual_indices) == 1) {
    series_boundary_second <- matrix(
      px_neighbors$series_boundary[as.character(residual_indices), ],
      nrow = 1
    )

    dimnames(series_boundary_second)[[1]] <- residual_indices
  }

  series_boundary <- rbind(series_boundary_first, series_boundary_second)

  # Create regressors.
  if (is.null(rr_xx) || sum(model$beta_coeffs) == 0) {
    kk <- 0
    xx <- NULL
    regressor_names <- NULL

    if (sum(model$beta_coeffs) > 0) {
      return(list(error = "The model includes exogenous regressors, but rr_xx was not supplied."))
    }
  } else {
    kk <- sum(model$beta_coeffs)
  }

  regressor_names <- names(model$beta_coeffs)[model$beta_coeffs]

  if (kk > 0 && is.list(rr_xx)) {
    xx <- array(0, dim = c(kk, pp, tt))
    dimnames(xx) <- list(regressor_names, indices, time_index)

    for (regressor in regressor_names) {
      if (regressor == "trend") {
        xx[regressor, , ] <- matrix(
          rep(seq_len(tt) / tt, pp),
          byrow = TRUE,
          nrow = pp
        )
      } else if (regressor %in% names(rr_xx)) {
        if (is.null(rr_xx[[regressor]])) {
          return(list(error = "There is a NULL object in the rr_xx covariate list."))
        }

        xx[regressor, , ] <- as.matrix(
          rr_xx[[regressor]][indices, as.character(time_index)]
        )
      } else {
        return(list(error = "The rr_xx argument does not contain all covariates required by the SDP-D model."))
      }
    }
  } else if (kk > 0 && is.array(rr_xx)) {
    xx <- rr_xx[regressor_names, indices, ]
  }

  # Output.
  list(
    series = series,
    xx = xx,
    ww_index = ww_index[indices, ],
    ww_values = ww_values[indices, ],
    px_neighbors = list(
      index = neighbors,
      series_boundary = series_boundary
    ),
    px = px,
    lon = lon,
    lat = lat,
    group = groups
  )
}
