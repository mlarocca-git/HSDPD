#' Parallelize SDP-D Series Data by Group
#'
#' Splits an SDP-D series object into group-specific SDP-D series objects.
#'
#' The function groups spatial units according to the `COD` column in
#' `series_object$group` and builds one SDP-D series object per group. It is used
#' by [fit_sdpd_model()] when `parallelize = TRUE`.
#'
#' @param series_object SDP-D series object, typically created with
#'   [build_sdpd_series()]. It must contain `series`, `xx`, `px`, `lon`, `lat`,
#'   `group`, `ww_index`, `ww_values`, and `px_neighbors`.
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()].
#'
#' @return A list of group-specific SDP-D series objects. If there is only one
#'   group, the original series object is returned inside a single-element list.
#'
#' @details
#' Parallelization is based on the group labels stored in
#' `series_object$group[, "COD"]`.
#'
#' Each group-specific object is rebuilt with [build_sdpd_series()] using the
#' full series, covariates, coordinates, spatial weights, and neighbor structure.
#'
#' @examples
#' \dontrun{
#' grouped_series <- parallelize_sdpd_data(
#'   series_object = series,
#'   model = model
#' )
#'}
#' @seealso [build_sdpd_series()], [fit_sdpd_model()]
#'
#' @export
parallelize_sdpd_data <- function(series_object, model) {
  # Group pixels and data in series_object by group.
  n_groups <- length(table(series_object$group[, "COD"]))

  if (n_groups > 1) {
    px <- series_object$px

  group_data <- tibble::tibble(
    group_code = series_object$group[as.character(px), "COD"],
    group_label = series_object$group[as.character(px), "LABEL"],
    px = as.numeric(px)
    ) |>
  dplyr::group_by(.data$group_code, .data$group_label) |>
  tidyr::nest()

    # For each group of pixels, prepare data and quantities for H-SDPD model
    # estimation.
  data_frame_data <- purrr::map(
    group_data$data,
    build_sdpd_series,
    rr_y = series_object$series,
    rr_xx = series_object$xx,
    lon = series_object$lon,
    lat = series_object$lat,
    rr_groups = series_object$group,
    model = model,
    check = FALSE,
    ww_index = series_object$ww_index,
    ww_values = series_object$ww_values,
    px_neighbors = series_object$px_neighbors
)
  } else {
    cat("\nWarning: no parallelization has been performed because there is only one group.")

    data_frame_data <- list(series_object)
  }

  # Return list of objects.
  data_frame_data
}


