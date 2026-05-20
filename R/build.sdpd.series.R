#' Build an SDP-D Series Object
#'
#' Builds an `sdpd_series` object, which is the data object expected by
#' subsequent estimation functions in the package.
#'
#' The function reads observed data from a data frame, matrix, or raster object.
#' If `sim = TRUE`, the series is simulated using the model structure supplied
#' through `model` and the simulation settings supplied through `vec_options`.
#'
#' @param df_obj Optional list-like object. Used mainly to support parallelized
#'   calls. If `px` is `NULL`, `px` is taken from `df_obj$px`.
#' @param px Optional object defining the spatial units or pixel structure.
#' @param rr_y Data frame, matrix, or `SpatRaster`. Endogenous variable series.
#'   Required unless `sim = TRUE`.
#' @param rr_xx Optional data frame, matrix, or raster-like object containing
#'   exogenous covariate series.
#' @param lon Optional vector or column name identifying longitude coordinates.
#' @param lat Optional vector or column name identifying latitude coordinates.
#' @param rr_groups Optional matrix or data frame defining regional groups.
#' @param label_groups Optional group labels, used mainly when reading raster
#'   data.
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()].
#' @param check Logical scalar. Whether the resulting series should be checked
#'   using [check_sdpd_series()]. Defaults to `TRUE`.
#' @param ww_index Optional object defining the index structure of the spatial
#'   weight matrix.
#' @param ww_values Optional object defining the values of the spatial weight
#'   matrix.
#' @param px_neighbors Optional object defining pixel or proximity neighbors.
#' @param sim Logical scalar. Whether the data series should be simulated.
#'   Defaults to `FALSE`.
#' @param nn Optional integer. Temporal dimension or number of observations used
#'   when simulating the series.
#' @param vec_options List. Options used when reading or simulating vectorized
#'   data. Expected elements include `px_core`, `px_neighbors`, `na_rm`,
#'   `na_covs`, `covariates_sim_model`, `markovian`, and `num_steps`.
#'
#' @return An SDP-D series object. The exact structure depends on the input data
#'   source and on whether the object is read from observed data or simulated.
#'   If `check = TRUE`, the returned object may also contain a `warnings`
#'   element. If validation fails, a list with `errors` and `warnings` is
#'   returned.
#'
#' @details
#' If `sim = FALSE`, `rr_y` must be supplied.
#'
#' If `sim = TRUE`, the function calls [generate_sdpd_series()] and uses the
#' spatial structure stored in `model`, including `ww_index`, `ww_values`,
#' `px_neighbors`, and `groups`.
#'
#' Data are read using [read_data_from_dataframe()] when `rr_y` is a data frame
#' or matrix, and using [read_data_from_raster()] when `rr_y` inherits from
#' class `SpatRaster`.
#'
#' @examples
#' \dontrun{
#' series <- build_sdpd_series(
#'   px = px,
#'   rr_y = rr_y,
#'   rr_xx = rr_xx,
#'   model = model
#' )
#'
#' simulated_series <- build_sdpd_series(
#'   model = model,
#'   sim = TRUE,
#'   nn = 50
#' )
#' }
#'
#' @seealso
#' [build_sdpd_model()],
#' [generate_sdpd_series()],
#' [check_sdpd_series()],
#' [read_data_from_dataframe()],
#' [read_data_from_raster()]
#'
#' @export
build_sdpd_series <- function(df_obj = NULL,
                              px = NULL,
                              rr_y = NULL,
                              rr_xx = NULL,
                              lon = NULL,
                              lat = NULL,
                              rr_groups = NULL,
                              label_groups = NULL,
                              model,
                              check = TRUE,
                              ww_index = NULL,
                              ww_values = NULL,
                              px_neighbors = NULL,
                              sim = FALSE,
                              nn = NULL,
                              vec_options = list(
                                px_core = 1,
                                px_neighbors = 6,
                                na_rm = TRUE,
                                na_covs = "pairwise.complete.obs",
                                covariates_sim_model = list(ar = c(0.8), sd = 1),
                                markovian = TRUE,
                                num_steps = 10
                              )) {

  ## This function builds the sdpd_series object expected by the subsequent
  ## estimation functions.
  ##
  ## Data are simulated when they are not supplied as input, provided that
  ## sim = TRUE.
  ##
  ## Note that the only argument without a default value is the sdpd_model object,
  ## built with build_sdpd_model().
  ##
  ## The first object, df_obj, is used to allow this function to be called in
  ## parallelized workflows.
  ##
  ## If sim = TRUE, any rr_y series supplied by the user is used only to
  ## initialize the simulated series.

  if (is.null(px)) {
    px <- df_obj$px
  }

  if (!sim && is.null(rr_y)) {
    return(list(error = "The rr_y object is required unless sim = TRUE"))
  }

  if (sim) {
    new_data <- generate_sdpd_series(
      nn = nn,
      rr_y = rr_y,
      rr_xx = rr_xx,
      model = model,
      markovian = vec_options$markovian,
      num_steps = vec_options$num_steps,
      covariates_sim_model = vec_options$covariates_sim_model
    )

    rr_y <- new_data$rr_y
    rr_xx <- new_data$rr_xx
    ww_index <- model$ww_index
    ww_values <- model$ww_values
    px_neighbors <- model$px_neighbors
    rr_groups <- model$groups
  }

  if (is.data.frame(rr_y) || is.matrix(rr_y)) {
    result <- read_data_from_dataframe(
      px = px,
      rr_y = rr_y,
      rr_xx = rr_xx,
      model = model,
      lon = lon,
      lat = lat,
      rr_groups = rr_groups,
      ww_index = ww_index,
      ww_values = ww_values,
      px_neighbors = px_neighbors
    )
  } else if (inherits(rr_y, "SpatRaster")) {
    result <- read_data_from_raster(
      px = px,
      rr_y = rr_y,
      rr_xx = rr_xx,
      lon = lon,
      lat = lat,
      rr_groups = rr_groups,
      model = model,
      label_groups = label_groups,
      vec_options = vec_options
    )
  } else {
    return(list(error = "The rr_y and rr_xx objects have an invalid format: they must be a data frame, matrix, or raster"))
  }

  if (!is.null(result$error)) {
    return(result)
  }

  if (is.null(nn) && !is.null(result$series)) {
    nn <- ncol(result$series)
  }

  result$nn <- nn

  # Diagnostics.
  if (check) {
    series_check <- check_sdpd_series(series = result, model = model)

    if (length(series_check$errors) > 0) {
      return(list(
        errors = series_check$errors,
        warnings = series_check$warnings
      ))
    }

    result$warnings <- series_check$warnings
  }

  # Output.
  result
}
