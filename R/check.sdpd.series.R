#' Check an SDP-D Series Object
#'
#' Checks the consistency of an SDP-D series object with an SDP-D model object.
#'
#' The function receives either an `sdpd_series` object and an `sdpd_model`
#' object, typically created with [build_sdpd_series()] and
#' [build_sdpd_model()], or the individual components of an SDP-D series.
#'
#' It validates the data series, regressors, spatial matrix components,
#' coordinates, groups, weights, fixed effects, and time effects. It returns the
#' checked objects together with vectors of errors and warnings.
#'
#' @param series Matrix, data frame, or list. The endogenous series, or an
#'   `sdpd_series`-like list containing the series components.
#' @param xx Optional matrix, data frame, or three-dimensional array. Regressor
#'   data. Defaults to `NULL`.
#' @param model Optional SDP-D model object, typically created with
#'   [build_sdpd_model()].
#' @param ww_index Optional matrix. Spatial-neighbor index matrix.
#' @param ww_values Optional matrix. Spatial-weight value matrix.
#' @param px_neighbors Optional object defining pixel or proximity neighbors.
#' @param px Optional object defining the spatial units or pixel structure.
#' @param lat Optional named numeric vector of latitude coordinates.
#' @param lon Optional named numeric vector of longitude coordinates.
#' @param group Optional matrix or data frame defining regional groups. It must
#'   contain two columns named `COD` and `LABEL`.
#' @param index_weights Optional vector, matrix, or data frame of spatial-unit
#'   weights.
#' @param time_weights Optional vector of time weights.
#'
#' @return A list containing the checked series components:
#' \describe{
#'   \item{series}{Checked and possibly transformed endogenous series.}
#'   \item{ww_index}{Spatial-neighbor index matrix.}
#'   \item{ww_values}{Spatial-weight value matrix.}
#'   \item{xx}{Checked and possibly transformed regressors.}
#'   \item{px_neighbors}{Pixel or proximity-neighbor information.}
#'   \item{time_effects}{Estimated or empty time effects.}
#'   \item{nn}{Number of time observations.}
#'   \item{pp}{Number of spatial units.}
#'   \item{kk}{Number of covariates.}
#'   \item{na}{Missing-value counts.}
#'   \item{px}{Spatial-unit or pixel structure.}
#'   \item{lat}{Latitude vector.}
#'   \item{lon}{Longitude vector.}
#'   \item{group}{Group structure.}
#'   \item{index_weights}{Spatial-unit weights.}
#'   \item{time_weights}{Time weights.}
#'   \item{mu}{Fixed-effect means.}
#'   \item{errors}{Character vector of validation errors.}
#'   \item{warnings}{Character vector of validation warnings.}
#' }
#'
#' @details
#' If `series` is a list, the function assumes that all series components are
#' passed through that list. This supports parallelized workflows.
#'
#' If the model includes time effects, the function removes time means from the
#' series and stores them in `time_effects`.
#'
#' If the model includes fixed effects, the function computes spatial-unit means
#' and stores them in `mu`.
#'
#' Regressors are mean-centered when supplied.
#'
#' @examples
#' \dontrun{
#' checked_series <- check_sdpd_series(
#'   series = series,
#'   model = model
#' )
#' }
#'
#' @seealso [build_sdpd_series()], [build_sdpd_model()]
#'
#' @export
check_sdpd_series <- function(series,
                              xx = NULL,
                              model = NULL,
                              ww_index = NULL,
                              ww_values = NULL,
                              px_neighbors = NULL,
                              px = NULL,
                              lat = NULL,
                              lon = NULL,
                              group = NULL,
                              index_weights = NULL,
                              time_weights = NULL) {

  ## This function receives an sdpd_series object and an sdpd_model object
  ## built with build_sdpd_series() and build_sdpd_model(), respectively,
  ## or receives the individual components of an sdpd_series object.
  ##
  ## It checks the consistency of the series components with the model and
  ## returns the checked objects together with errors and warnings.

  error_vector <- character(20)
  warning_vector <- character(20)

  n_error <- 0
  n_warning <- 0

  model_obj <- model
  ww_index_temp <- NULL
  ww_values_temp <- NULL
  nn <- NULL
  pp <- NULL

  # Check validity of the series.
  if (is.list(series)) {
    # In this case, all objects are passed through a list, as in parallelized
    # workflows.
    data_series <- series$series
    px <- series$px
    lon <- series$lon
    lat <- series$lat
    group <- series$group
    px_neighbors <- series$px_neighbors
    ww_index_temp <- series$ww_index
    ww_values_temp <- series$ww_values
    nn <- series$nn

    if (is.null(nn) && !is.null(series$metadata$nn)) {
      nn <- series$metadata$nn
    }

    if (is.null(xx)) {
      xx <- series$xx
    }
  } else {
    data_series <- series
  }

  if (is.matrix(data_series) || is.data.frame(data_series)) {
    data_series <- as.matrix(data_series)

    data_nn <- dim(data_series)[2]
    pp <- dim(data_series)[1]

    if (is.null(nn)) {
      nn <- data_nn
    } else if (!identical(as.integer(nn), as.integer(data_nn))) {
      n_error <- n_error + 1
      error_vector[n_error] <- "The series nn value must match ncol(series)."
    }

    if (is.null(dimnames(data_series)[[1]]) ||
        is.null(dimnames(data_series)[[2]])) {
      n_error <- n_error + 1
      error_vector[n_error] <- "\nMissing names for data_series. They should be numbers for compatibility with plot functions."
    }

    if (!is.null(time_weights)) {
      if (length(time_weights) != nn) {
        time_weights <- NULL
        n_warning <- n_warning + 1
        warning_vector[n_warning] <- "\nThere is a problem with time_weights; they have been set to NULL."
      }
    }

    if (!is.null(index_weights)) {
      if ((is.matrix(index_weights) || is.data.frame(index_weights)) &&
          dim(index_weights)[2] > 1) {
        if (is.null(dimnames(index_weights)[[1]])) {
          dimnames(index_weights)[[1]] <- dimnames(data_series)[[1]]
        }

        spatial_weights <- index_weights[dimnames(data_series)[[1]], 2]
        names(spatial_weights) <- index_weights[dimnames(data_series)[[1]], 1]
        index_weights <- spatial_weights
      } else {
        index_weights <- NULL
        n_warning <- n_warning + 1
        warning_vector[n_warning] <- "\nThere is a problem with index_weights; it has been set to NULL."
      }
    }
  } else {
    n_error <- n_error + 1
    error_vector[n_error] <- "The series is not a matrix or data frame."

    if (is.null(nn)) {
      n_error <- n_error + 1
      error_vector[n_error] <- "The series object is missing nn and it cannot be inferred from series$series."
      nn <- 0L
    }

    pp <- 0L
    data_series <- matrix(numeric(0), nrow = pp, ncol = nn)
    dimnames(data_series) <- list(character(pp), character(nn))
  }

  # Check validity of regressors.
  if (is.null(model_obj)) {
    n_error <- n_error + 1
    error_vector[n_error] <- "The model is missing. Please check."

    covariate_names <- NULL
  } else {
    model_obj$kk <- sum(model_obj$beta_coeffs)
    covariate_names <- names(model$beta_coeffs)[model$beta_coeffs]
  }

  if (is.null(xx)) {
    if (!is.null(model_obj$kk) && model_obj$kk > 0) {
      n_error <- n_error + 1
      error_vector[n_error] <- "The model has exogenous covariates, but no data were passed through xx."
    }

    kk <- 0
  } else if (model_obj$kk == 0) {
    n_warning <- n_warning + 1
    warning_vector[n_warning] <- "The model has no exogenous covariates, so the object passed through xx has been ignored."

    kk <- 0
  } else if (is.matrix(xx) || is.data.frame(xx)) {
    if (dim(xx)[2] != nn) {
      n_error <- n_error + 1
      error_vector[n_error] <- "The regressor must have the same number of observations, i.e. columns, as the series."
    }

    if (dim(xx)[1] != pp) {
      n_error <- n_error + 1
      error_vector[n_error] <- "The regressor must have the same number of locations, i.e. rows, as the series."
    }

    xx <- as.matrix(xx)
    xx <- apply(xx, 1, FUN = function(x) x - mean(x))
    xx <- t(xx)

    n_warning <- n_warning + 1
    warning_vector[n_warning] <- "The regressor has been mean-centered."

    if (is.null(dimnames(xx)[[1]]) || is.null(dimnames(xx)[[2]])) {
      n_error <- n_error + 1
      error_vector[n_error] <- "Dimnames for xx are missing."
    }

    kk <- 1
  } else if (is.array(xx) && length(dim(xx)) == 3) {
    if (dim(xx)[3] != nn) {
      n_error <- n_error + 1
      error_vector[n_error] <- paste(
        "The regressors must have the same number of observations as the series: ",
        nn,
        " instead of ",
        dim(xx)[3],
        ".",
        sep = ""
      )
    }

    if (dim(xx)[2] != pp) {
      n_error <- n_error + 1
      error_vector[n_error] <- paste(
        "The regressors must have the same number of locations as the series: ",
        pp,
        " instead of ",
        dim(xx)[2],
        ".",
        sep = ""
      )
    }

    xx <- apply(xx, c(1, 2), FUN = function(x) x - mean(x))
    xx <- aperm(xx, c(2, 3, 1))

    n_warning <- n_warning + 1
    warning_vector[n_warning] <- "The regressors have been mean-centered."

    if (is.null(dimnames(xx)[[1]]) ||
        is.null(dimnames(xx)[[2]]) ||
        is.null(dimnames(xx)[[3]])) {
      n_error <- n_error + 1
      error_vector[n_error] <- "Dimnames for xx are missing."
    }

    kk <- dim(xx)[1]
  } else {
    n_error <- n_error + 1
    error_vector[n_error] <- "Something is wrong with regressor xx. It must be a matrix, data frame, or three-dimensional array."

    kk <- 0
  }

  if (kk > 0 && length(dim(xx)) == 3) {
    covariate_index <- dimnames(xx)[[1]] %in% covariate_names
    xx <- xx[covariate_index, , , drop = FALSE]
    kk <- dim(xx)[1]

    if (kk != model_obj$kk) {
      n_error <- n_error + 1
      error_vector[n_error] <- "Some covariates in the model are missing from xx."
    }
  }

  if (kk == 1 && length(dim(xx)) == 3) {
    xx_dimnames <- dimnames(xx)
    xx <- xx[1, , , drop = TRUE]
    dimnames(xx) <- xx_dimnames[2:3]
  }

  if (kk == 1 && length(dim(xx)) == 2) {
    if (is.null(names(model_obj$beta_coeffs))) {
      n_error <- n_error + 1
      error_vector[n_error] <- "The name of the covariate in the model is missing."
    }
  }

  # Check validity of the spatial matrix.
  if (is.null(ww_index)) {
    ww_index <- ww_index_temp
  }

  if (is.null(ww_values)) {
    ww_values <- ww_values_temp
  }

  if (!is.null(model) && (is.null(ww_index) || is.null(ww_values))) {
    ww_index <- model$ww_index
    ww_values <- model$ww_values
  }

  if (is.null(ww_index) || is.null(ww_values)) {
    n_error <- n_error + 1
    error_vector[n_error] <- "The spatial matrix components are missing."
  } else if (is.null(dimnames(ww_index)[[1]]) ||
             is.null(dimnames(ww_values)[[1]])) {
    n_error <- n_error + 1
    error_vector[n_error] <- "Names for the spatial matrix are missing."
  }

  # Check missing values.
  na_values <- sum(is.na(data_series))

  if (na_values > 0) {
    n_error <- n_error + 1
    error_vector[n_error] <- "The series has NA values."
  }

  if (kk == 0) {
    na_x <- 0
  } else if (kk == 1) {
    na_x <- sum(is.na(xx))
  } else if (kk > 1) {
    na_x <- apply(xx, 1, FUN = function(x) sum(is.na(x)))
  }

  if (sum(na_x) != 0) {
    n_error <- n_error + 1
    error_vector[n_error] <- "There are missing values in the covariates xx."

    na_values <- c(na_values, na_x)
    names(na_values) <- c("na_y", paste("na_x", seq_len(kk), sep = ""))
  }

  if (model_obj$time_effects) {
    time_effects <- apply(data_series, 2, mean)
    data_series <- data_series - matrix(
      rep(time_effects, pp),
      nrow = pp,
      byrow = TRUE
    )
  } else {
    time_effects <- numeric(nn)
  }

  # Check lon, lat, and group.
  if (!is.null(lon)) {
    if (!is.vector(lon) || !is.numeric(lon) || is.null(names(lon))) {
      n_error <- n_error + 1
      error_vector[n_error] <- "Something is wrong with the longitude vector. Are names missing?"
    }
  }

  if (!is.null(lat)) {
    if (!is.vector(lat) || !is.numeric(lat) || is.null(names(lat))) {
      n_error <- n_error + 1
      error_vector[n_error] <- "Something is wrong with the latitude vector. Are names missing?"
    }
  }

  if (!is.null(group)) {
    if (length(dim(group)) != 2 ||
        dim(group)[2] != 2 ||
        is.null(dimnames(group)[[1]])) {
      n_error <- n_error + 1
      error_vector[n_error] <- "Something is wrong with the group object. Are names missing?"
    } else if (sum(c("COD", "LABEL") %in% dimnames(group)[[2]]) < 2) {
      n_error <- n_error + 1
      error_vector[n_error] <- "The group object must contain the columns COD and LABEL."
    }
  }

  if (model_obj$fixed_effects) {
    mu <- apply(data_series, 1, mean)
  } else {
    mu <- numeric(length(px))
    names(mu) <- dimnames(data_series)[[1]]
  }

  # Return checked data structure.
  list(
    series = data_series,
    ww_index = ww_index,
    ww_values = ww_values,
    xx = xx,
    px_neighbors = px_neighbors,
    time_effects = time_effects,
    nn = nn,
    pp = pp,
    kk = kk,
    na = na_values,
    px = px,
    lat = lat,
    lon = lon,
    group = group,
    index_weights = index_weights,
    time_weights = time_weights,
    mu = mu,
    errors = error_vector[error_vector != ""],
    warnings = warning_vector[warning_vector != ""]
  )
}
