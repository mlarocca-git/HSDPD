.new_sdpd_design <- function(series = NULL,
                             x = NULL,
                             unit_index = NULL,
                             time_index = NULL,
                             pp = NULL,
                             nn = NULL,
                             kk = NULL,
                             ww_index = NULL,
                             ww_values = NULL,
                             ww = NULL,
                             px_neighbors = NULL,
                             mu = NULL,
                             time_effects = NULL,
                             lat = NULL,
                             lon = NULL,
                             group = NULL,
                             index_weights = NULL,
                             time_weights = NULL,
                             source_px = NULL,
                             model = NULL,
                             errors = character(),
                             warnings = character()) {
  structure(
    list(
      series = series,
      x = x,
      unit_index = unit_index,
      time_index = time_index,
      pp = pp,
      nn = nn,
      kk = kk,
      ww_index = ww_index,
      ww_values = ww_values,
      ww = ww,
      px_neighbors = px_neighbors,
      mu = mu,
      time_effects = time_effects,
      lat = lat,
      lon = lon,
      group = group,
      index_weights = index_weights,
      time_weights = time_weights,
      source_px = source_px,
      model = model,
      errors = errors,
      warnings = warnings
    ),
    class = "sdpd_design"
  )
}

.validate_sdpd_design <- function(design) {
  errors <- design$errors
  warnings <- design$warnings
  model_beta_names <- NULL
  model_kk <- NULL

  if (!inherits(design, "sdpd_design")) {
    errors <- c(errors, "The design object must inherit from class sdpd_design.")
    design$errors <- errors
    design$warnings <- warnings
    return(design)
  }

  if (!is.null(design$model) && !is.null(design$model$beta_coeffs)) {
    model_kk <- sum(design$model$beta_coeffs)
    model_beta_names <- names(design$model$beta_coeffs)[design$model$beta_coeffs]
  }

  if (!is.matrix(design$series) || !is.numeric(design$series)) {
    errors <- c(errors, "The design series must be a numeric matrix.")
  } else {
    if (is.null(rownames(design$series)) || is.null(colnames(design$series))) {
      errors <- c(errors, "The design series must have row and column names.")
    }

    if (!identical(nrow(design$series), design$pp)) {
      errors <- c(errors, "The design pp value must match nrow(series).")
    }

    if (!identical(ncol(design$series), design$nn)) {
      errors <- c(errors, "The design nn value must match ncol(series).")
    }

    if (is.null(design$unit_index)) {
      errors <- c(errors, "The design unit_index must not be NULL.")
    } else if (!identical(rownames(design$series), design$unit_index)) {
      errors <- c(errors, "The design unit_index must match rownames(series).")
    }

    if (is.null(design$time_index)) {
      errors <- c(errors, "The design time_index must not be NULL.")
    } else if (!identical(colnames(design$series), design$time_index)) {
      errors <- c(errors, "The design time_index must match colnames(series).")
    }
  }

  if (!is.null(model_kk) && !isTRUE(design$kk == model_kk)) {
    errors <- c(errors, "The design kk value must match the model covariate count.")
  }

  if (!is.null(design$kk) && length(design$kk) == 1) {
    if (design$kk == 0 && !is.null(design$x)) {
      errors <- c(errors, "A zero-covariate design must store x as NULL.")
    }

    if (design$kk > 0 && is.null(design$x)) {
      errors <- c(errors, "A covariate design must include x.")
    }

    if (design$kk == 1 && !is.null(design$x)) {
      if (!is.matrix(design$x)) {
        errors <- c(errors, "A one-covariate design must store x as a matrix.")
      } else {
        if (!all(dim(design$x) == c(design$pp, design$nn))) {
          errors <- c(errors, "A one-covariate design x matrix must have dimensions pp by nn.")
        }

        if (!identical(rownames(design$x), design$unit_index)) {
          errors <- c(errors, "A one-covariate design x row names must match unit_index.")
        }

        if (!identical(colnames(design$x), design$time_index)) {
          errors <- c(errors, "A one-covariate design x column names must match time_index.")
        }
      }
    }

    if (design$kk > 1 && !is.null(design$x)) {
      if (!is.array(design$x) || length(dim(design$x)) != 3) {
        errors <- c(errors, "A multi-covariate design must store x as a three-dimensional array.")
      } else {
        if (!all(dim(design$x)[2:3] == c(design$pp, design$nn))) {
          errors <- c(errors, "A multi-covariate design x array must have dimensions covariates by pp by nn.")
        }

        if (!identical(dimnames(design$x)[[2]], design$unit_index)) {
          errors <- c(errors, "A multi-covariate design x location names must match unit_index.")
        }

        if (!identical(dimnames(design$x)[[3]], design$time_index)) {
          errors <- c(errors, "A multi-covariate design x time names must match time_index.")
        }

        if (!is.null(model_beta_names) &&
            !is.null(dimnames(design$x)[[1]]) &&
            !identical(dimnames(design$x)[[1]], model_beta_names)) {
          errors <- c(errors, "A multi-covariate design x covariate names must match the model covariates.")
        }
      }
    }
  } else if (!is.null(design$x)) {
    errors <- c(errors, "The design kk value must be defined when x is supplied.")
  }

  if (!is.null(design$ww)) {
    if (!is.matrix(design$ww) || !is.numeric(design$ww)) {
      errors <- c(errors, "The design spatial weights object must be a numeric matrix.")
    } else {
      if (!all(dim(design$ww) == c(design$pp, design$pp))) {
        errors <- c(errors, "The design spatial weights object must have dimensions pp by pp.")
      }

      if (!identical(rownames(design$ww), design$unit_index) ||
          !identical(colnames(design$ww), design$unit_index)) {
        errors <- c(errors, "The design spatial weights names must match unit_index.")
      }
    }
  }

  if (!is.null(design$ww_index) || !is.null(design$ww_values)) {
    if (is.null(design$ww_index) || is.null(design$ww_values)) {
      errors <- c(errors, "The design ww_index and ww_values objects must both be present.")
    } else if (!is.matrix(design$ww_index) || !is.matrix(design$ww_values)) {
      errors <- c(errors, "The design ww_index and ww_values objects must be matrices.")
    } else {
      if (!identical(dim(design$ww_index), dim(design$ww_values))) {
        errors <- c(errors, "The design ww_index and ww_values objects must have matching dimensions.")
      }

      if (!is.null(rownames(design$ww_index)) &&
          !identical(rownames(design$ww_index), design$unit_index)) {
        errors <- c(errors, "The design ww_index row names must match unit_index.")
      }

      if (!is.null(rownames(design$ww_values)) &&
          !identical(rownames(design$ww_values), design$unit_index)) {
        errors <- c(errors, "The design ww_values row names must match unit_index.")
      }

      if (!is.null(colnames(design$ww_index)) &&
          !is.null(colnames(design$ww_values)) &&
          !identical(colnames(design$ww_index), colnames(design$ww_values))) {
        errors <- c(errors, "The design ww_index and ww_values column names must match.")
      }
    }
  }

  if (!is.null(design$px_neighbors)) {
    if (!is.list(design$px_neighbors)) {
      errors <- c(errors, "The design px_neighbors object must be a list.")
    } else if (is.null(design$px_neighbors$index)) {
      errors <- c(errors, "The design px_neighbors object must contain an index matrix.")
    } else if (!is.matrix(design$px_neighbors$index)) {
      errors <- c(errors, "The design px_neighbors index must be a matrix.")
    } else {
      if (!is.null(rownames(design$px_neighbors$index)) &&
          !identical(rownames(design$px_neighbors$index), design$unit_index)) {
        errors <- c(errors, "The design px_neighbors index row names must match unit_index.")
      }
    }

    if (is.list(design$px_neighbors) &&
        !is.null(design$px_neighbors$series_boundary)) {
      if (!is.matrix(design$px_neighbors$series_boundary) ||
          !is.numeric(design$px_neighbors$series_boundary)) {
        errors <- c(errors, "The design px_neighbors series_boundary must be a numeric matrix.")
      } else {
        boundary_rownames <- rownames(design$px_neighbors$series_boundary)

        if (!is.null(boundary_rownames) && any(!nzchar(boundary_rownames))) {
          errors <- c(errors, "The design px_neighbors series_boundary row names must be non-empty.")
        }

        if (!is.null(colnames(design$px_neighbors$series_boundary)) &&
            !identical(colnames(design$px_neighbors$series_boundary), design$time_index)) {
          errors <- c(errors, "The design px_neighbors series_boundary column names must match time_index.")
        }
      }
    }
  }

  if (!is.null(design$mu)) {
    if (!is.numeric(design$mu)) {
      errors <- c(errors, "The design mu vector must be numeric.")
    }

    if (length(design$mu) != design$pp) {
      errors <- c(errors, "The design mu vector must have length pp.")
    }

    if (!is.null(names(design$mu)) &&
        !identical(names(design$mu), design$unit_index)) {
      errors <- c(errors, "The design mu names must match unit_index.")
    }
  }

  if (!is.null(design$time_effects)) {
    if (!is.numeric(design$time_effects)) {
      errors <- c(errors, "The design time_effects vector must be numeric.")
    }

    if (length(design$time_effects) != design$nn) {
      errors <- c(errors, "The design time_effects vector must have length nn.")
    }

    if (!is.null(names(design$time_effects)) &&
        !identical(names(design$time_effects), design$time_index)) {
      errors <- c(errors, "The design time_effects names must match time_index.")
    }
  }

  design$errors <- unique(errors)
  design$warnings <- unique(warnings)
  design
}

.as_sdpd_design <- function(series,
                            model,
                            index_weights = NULL,
                            time_weights = NULL) {
  checked <- check_sdpd_series(
    series = series,
    model = model,
    index_weights = index_weights,
    time_weights = time_weights
  )

  if (length(checked$errors) > 0) {
    return(.new_sdpd_design(
      model = model,
      errors = checked$errors,
      warnings = checked$warnings
    ))
  }

  ww <- build_spatial_matrix(
    ww_index = checked$ww_index,
    ww_values = checked$ww_values
  )

  errors <- checked$errors

  if (is.list(ww) && !is.null(ww$error)) {
    errors <- c(errors, ww$error)
    ww <- NULL
  }

  design <- .new_sdpd_design(
    series = checked$series,
    x = checked$xx,
    unit_index = rownames(checked$series),
    time_index = colnames(checked$series),
    pp = checked$pp,
    nn = checked$nn,
    kk = checked$kk,
    ww_index = checked$ww_index,
    ww_values = checked$ww_values,
    ww = ww,
    px_neighbors = checked$px_neighbors,
    mu = checked$mu,
    time_effects = checked$time_effects,
    lat = checked$lat,
    lon = checked$lon,
    group = checked$group,
    index_weights = checked$index_weights,
    time_weights = checked$time_weights,
    source_px = checked$px,
    model = model,
    errors = errors,
    warnings = checked$warnings
  )

  .validate_sdpd_design(design)
}
