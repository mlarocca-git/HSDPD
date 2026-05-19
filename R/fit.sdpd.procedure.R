#' Run the SDP-D Estimation Procedure
#'
#' Runs the full SDP-D estimation procedure on a single SDP-D series object.
#'
#' The function checks the input series, computes covariance components, builds
#' the spatial weight matrix, estimates model coefficients, optionally applies
#' the second-stage estimator, computes fitted values and residuals, estimates
#' mean-equation summaries, and optionally checks stationarity diagnostics.
#'
#' @param series SDP-D series object, typically created with
#'   [build_sdpd_series()].
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()].
#' @param check Logical scalar. Whether stationarity diagnostics should be
#'   computed for the estimated model. Defaults to `FALSE`.
#' @param two_stage Logical scalar. Whether the second-stage estimator should be
#'   applied after the first-stage coefficient estimation. Defaults to `FALSE`.
#' @param na_covs Character scalar. Missing-value handling method passed to
#'   covariance computations. Defaults to `"pairwise.complete.obs"`.
#'
#' @return A list containing estimation outputs, including spatial metadata,
#'   coefficient estimates, fitted values, residuals, mean-equation summaries,
#'   time effects, model object, warnings, and diagnostics. If validation fails,
#'   a list with `errors` and `warnings` is returned.
#'
#' @details
#' This function is the core single-series estimation routine used by
#' [fit_sdpd_model()]. For parallelized workflows, it is called separately on
#' each group-specific SDP-D series object.
#'
#' @examples
#' \dontrun{
#' result <- fit_sdpd_procedure(
#'   series = series,
#'   model = model
#' )
#'
#' result_checked <- fit_sdpd_procedure(
#'   series = series,
#'   model = model,
#'   check = TRUE,
#'   two_stage = TRUE
#' )
#'}
#' 
#' @seealso
#' [fit_sdpd_model()],
#' [check_sdpd_series()],
#' [fit_sdpd_covs()],
#' [build_spatial_matrix()],
#' [fit_sdpd_coefficients()],
#' [fit_second_stage()],
#' [fit_sdpd_series()],
#' [fit_sdpd_mean_equation_model()],
#' [check_sdpd_model()]
#'
#' @export
fit_sdpd_procedure <- function(series,
                               model,
                               check = FALSE,
                               two_stage = FALSE,
                               na_covs = "pairwise.complete.obs") {

  ## series can be an sdpd_series object created by build_sdpd_series(), or it
  ## can contain the data series matrix, where rows are spatial units and
  ## columns are time observations.

  if (is.null(series)) {
    cat("\nError: data is NULL.")
    return(NULL)
  }

  # Check validity of data and construct the internal numeric design.
  design <- .as_sdpd_design(series = series, model = model)

  if (length(design$errors) > 0) {
    cat("\nThere are errors in the series:\n", design$errors)
    return(list(
      errors = design$errors,
      warnings = design$warnings
    ))
  }

  covs <- fit_sdpd_covs(
    series = design$series,
    x = design$x,
    px_neighbors = design$px_neighbors,
    kk = design$kk,
    nn = design$nn,
    pp = design$pp,
    na_covs = na_covs
  )

  # Estimate model parameters.
  fit <- fit_sdpd_coefficients(
    ww = design$ww,
    covs = covs,
    mu = design$mu,
    model = model
  )

  if (sum(is.na(fit$coeff_hat)) > 0) {
    design$warnings <- c(
      design$warnings,
      "There are NA values in the estimated coefficients."
    )
    cat("\nThere are NA values in the estimated coefficients.\n")
  }

  # Second-stage estimation.
  if (two_stage) {
    fit$coeff_hat <- fit_second_stage(
      data_series = design$series,
      ww = design$ww,
      x_centered = design$x,
      model = model,
      coeff_hat = fit$coeff_hat
    )
  }

  # Estimate fitted values.
  fitted_result <- fit_sdpd_series(
    data_series = design$series,
    ww = design$ww,
    x_centered = design$x,
    model = model,
    coeff_hat = fit$coeff_hat,
    time_effects = design$time_effects
  )

  # Estimate the mean-equation meta-model.
  mean_equation <- fit_sdpd_mean_equation_model(
    result = fitted_result,
    ww = design$ww,
    time_weights = design$time_weights,
    index_weights = design$index_weights
  )

  result <- list(
    px = design$source_px,
    lon = design$lon,
    lat = design$lat,
    group = design$group,
    coeff_hat = fit$coeff_hat,
    fitted = fitted_result$fitted,
    resid = fitted_result$resid,
    mean_equation = mean_equation,
    time_effects = design$time_effects,
    model = model,
    warnings = design$warnings
  )

  # Check stationarity conditions for the estimated model.
  diagnostics <- check

  if (check) {
    diagnostics <- check_sdpd_model(
      res_fit = result,
      ww_index = design$ww_index,
      ww_values = design$ww_values
    )$diagnostics[1, "max_mod_eigen_a"]
  }

  # Output.
  result$diagnostics <- diagnostics
  attr(result, "sdpd_design") <- design
  result
}
