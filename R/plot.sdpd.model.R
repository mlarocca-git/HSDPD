#' Plot an SDP-D Model Fit
#'
#' Dispatches plotting of an SDP-D fitted model object.
#'
#' This function is a lightweight plotting wrapper. Currently, `which = 1`
#' dispatches to [plot_1_sdpd_model()].
#'
#' @param res_fit Fitted SDP-D model result object.
#' @param n_units Integer or `"all"`. Number of spatial units to plot. Defaults
#'   to `"all"`.
#' @param n_vars Integer or `"all"`. Number of variables to plot. Defaults to
#'   `"all"`.
#' @param which Integer scalar. Plot type selector. Currently only `1` is
#'   implemented. Defaults to `1`.
#' @param t_axis List. Time-axis options. Expected elements are `t_labels` and
#'   `t_points`.
#' @param x_limit Optional numeric vector of length 2. X-axis limits.
#' @param y_limit Optional numeric vector of length 2. Y-axis limits.
#' @param max_col Integer. Maximum number of plot columns. Defaults to `5`.
#' @param point_col Color specification for plotted points or lines. Defaults
#'   to `1`.
#'
#' @return The output of the selected plotting function, invisibly if that
#'   function returns invisibly.
#'
#' @details
#' The renamed API uses snake_case argument names. The old `t.axis` structure is
#' now `t_axis`, with elements `t_labels` and `t_points`.
#'
#' @examples
#' plot_sdpd_model(
#'   res_fit = fit,
#'   n_units = "all",
#'   n_vars = "all"
#' )
#'
#' @seealso [plot_1_sdpd_model()]
#'
#' @export
plot_sdpd_model <- function(res_fit,
                            n_units = "all",
                            n_vars = "all",
                            which = 1,
                            t_axis = list(t_labels = NULL, t_points = NULL),
                            x_limit = NULL,
                            y_limit = NULL,
                            max_col = 5,
                            point_col = 1) {
  if (which == 1) {
    return(plot_1_sdpd_model(
      res_fit = res_fit,
      n_units = n_units,
      n_vars = n_vars,
      t_axis = t_axis,
      x_limit = x_limit,
      y_limit = y_limit,
      max_col = max_col,
      point_col = point_col
    ))
  }

  stop("Only which = 1 is currently implemented.")
}

