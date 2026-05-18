#' Fit an SDP-D Model
#'
#' Fits an SDP-D model to an SDP-D series object.
#'
#' The function estimates H-SDPD model parameters either on the full series or,
#' when `parallelize = TRUE`, on parallelized groups of spatial units. Each
#' group-specific result is assembled into a standardized output and then row
#' bound into a single result object.
#'
#' @param series SDP-D series object, typically created with
#'   [build_sdpd_series()].
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()].
#' @param check Logical scalar. Whether model checks should be performed during
#'   fitting. Defaults to `FALSE`.
#' @param two_stage Logical scalar. Whether the second-stage estimator should be
#'   used. Defaults to `FALSE`.
#' @param na_covs Character scalar. Missing-value handling method passed to
#'   covariance computations. Defaults to `"pairwise.complete.obs"`.
#' @param parallelize Logical scalar. Whether the series should be split into
#'   parallelized spatial groups before estimation. Defaults to `FALSE`.
#'
#' @return A data frame or tibble containing assembled SDP-D estimation results.
#'
#' @details
#' If `parallelize = TRUE`, the function first calls
#' [parallelize_sdpd_data()] to split the SDP-D series into group-specific
#' objects. Otherwise, the full `series` object is wrapped in a single-element
#' list.
#'
#' Estimation is performed by [fit_sdpd_procedure()] and each group result is
#' assembled by [fit_sdpd_assemble()].
#'
#' @examples
#' \dontrun{
#' fit <- fit_sdpd_model(
#'   series = series,
#'   model = model
#' )
#'
#' fit_parallel <- fit_sdpd_model(
#'   series = series,
#'   model = model,
#'   parallelize = TRUE
#' )
#'}
#' @seealso
#' [build_sdpd_series()],
#' [build_sdpd_model()],
#' [fit_sdpd_procedure()],
#' [fit_sdpd_assemble()],
#' [parallelize_sdpd_data()]
#'
#' @export
fit_sdpd_model <- function(series,
                           model,
                           check = FALSE,
                           two_stage = FALSE,
                           na_covs = "pairwise.complete.obs",
                           parallelize = FALSE) {

  if (parallelize) {
    # Build a list of SDP-D series objects on parallelized groups.
    data_frame_data <- parallelize_sdpd_data(
      series_object = series,
      model = model
    )
  } else {
    data_frame_data <- list(series)
  }

  # Estimate the H-SDPD model parameters on each group of pixels.
  result <- data_frame_data |>
    purrr::map(
      fit_sdpd_procedure,
      model = model,
      check = check,
      two_stage = two_stage,
      na_covs = na_covs
    ) |>
    purrr::map(fit_sdpd_assemble) |>
    dplyr::bind_rows()

  # Output.
  result
}


