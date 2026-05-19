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

  if (!inherits(design, "sdpd_design")) {
    errors <- c(errors, "The design object must inherit from class sdpd_design.")
    design$errors <- errors
    design$warnings <- warnings
    return(design)
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
  }

  if (!is.null(design$x)) {
    if (design$kk == 1 && !is.matrix(design$x)) {
      errors <- c(errors, "A one-covariate design must store x as a matrix.")
    }

    if (design$kk > 1 &&
        (!is.array(design$x) || length(dim(design$x)) != 3)) {
      errors <- c(errors, "A multi-covariate design must store x as a three-dimensional array.")
    }
  }

  if (!is.null(design$ww)) {
    if (!is.matrix(design$ww) || !is.numeric(design$ww)) {
      errors <- c(errors, "The design spatial weights object must be a numeric matrix.")
    } else if (!identical(rownames(design$ww), design$unit_index) ||
               !identical(colnames(design$ww), design$unit_index)) {
      errors <- c(errors, "The design spatial weights names must match unit_index.")
    }
  }

  if (!is.null(design$mu) && length(design$mu) != design$pp) {
    errors <- c(errors, "The design mu vector must have length pp.")
  }

  if (!is.null(design$time_effects) &&
      length(design$time_effects) != design$nn) {
    errors <- c(errors, "The design time_effects vector must have length nn.")
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
    design <- .new_sdpd_design(
      model = model,
      errors = checked$errors,
      warnings = checked$warnings
    )

    return(.validate_sdpd_design(design))
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
