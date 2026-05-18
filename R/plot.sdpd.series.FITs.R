#' Plot Fitted, Observed, and Residual SDP-D Series
#'
#' Plots fitted values, observed values, and residuals for the spatial location
#' closest to the supplied latitude and longitude after coordinate rounding.
#'
#' The function extracts fitted values and residuals from an SDP-D result object,
#' reconstructs the observed series as `fitted + resid`, and displays the three
#' series in a single `ggplot`.
#'
#' @param results List or data frame. SDP-D estimation results containing
#'   `fitted`, `resid`, `lat`, and `lon`.
#' @param latitude Numeric scalar. Latitude of the location to plot.
#' @param longitude Numeric scalar. Longitude of the location to plot.
#' @param n_digits Integer. Number of digits used to round coordinates before
#'   matching. Defaults to `1`.
#' @param main Optional character scalar. Plot title.
#' @param sub Optional character scalar. Plot subtitle.
#' @param xlab Optional character scalar. X-axis label.
#' @param ylab Optional character scalar. Y-axis label.
#'
#' @return A `ggplot` object, or a character message if the requested
#'   coordinates are not present in the result object.
#'
#' @details
#' Residuals are shifted vertically by an offset so they can be displayed on the
#' same plot as the fitted and observed series. A secondary y-axis is added for
#' the residual scale.
#'
#' Time labels are inferred from `colnames(results$fitted)`. Numeric time labels
#' are converted to dates using [as.Date()].
#'
#' @examples
#' plot_sdpd_series_fits(
#'   results = fit,
#'   latitude = 45.1,
#'   longitude = 9.2
#' )
#'
#' @seealso [ggplot2::ggplot()]
#'
#' @export
plot_sdpd_series_fits <- function(results,
                                  latitude,
                                  longitude,
                                  n_digits = 1,
                                  main = NULL,
                                  sub = NULL,
                                  xlab = NULL,
                                  ylab = NULL) {
  time_index <- dimnames(results$fitted)[[2]]

  if (inherits(time_index, "Date")) {
    time_range <- as.Date(time_index)
  } else if (all(!is.na(as.numeric(time_index)))) {
    time_range <- as.Date(as.numeric(time_index))
  } else {
    time_range <- time_index
  }

  if (is.null(sub)) {
    sub <- paste(
      "\n(based on data from ",
      range(time_index)[1],
      " to ",
      range(time_index)[2],
      ")",
      sep = ""
    )
  }

  if (is.null(main)) {
    main <- paste(
      "Estimated model series (longitude=",
      longitude,
      " - latitude=",
      latitude,
      ")",
      sep = ""
    )
  }

  selected_series <- results |>
    dplyr::mutate(latitude_rounded = round(lat, n_digits)) |>
    dplyr::mutate(longitude_rounded = round(lon, n_digits)) |>
    dplyr::filter(
      latitude_rounded == round(latitude, n_digits),
      longitude_rounded == round(longitude, n_digits)
    ) |>
    dplyr::select(fitted, resid)

  if (dim(selected_series)[1] == 0) {
    return("These coordinates are not present in the database.")
  }

  offset_axis <- range(selected_series$fitted, na.rm = TRUE)[1] -
    0.1 * diff(range(selected_series$fitted, na.rm = TRUE))

  plot_data <- data.frame(
    time = time_range,
    fitted = as.numeric(selected_series$fitted),
    resid = as.numeric(selected_series$resid)
  ) |>
    dplyr::mutate(observed = resid + fitted)

  result <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = time)) +
    ggplot2::geom_line(ggplot2::aes(y = fitted, color = "Fitted")) +
    ggplot2::geom_line(ggplot2::aes(y = observed, color = "Observed")) +
    ggplot2::geom_line(
      ggplot2::aes(y = resid + offset_axis, color = "Residuals")
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = mean(resid, na.rm = TRUE) + offset_axis)
    ) +
    ggplot2::guides(x = ggplot2::guide_axis(angle = 0)) +
    ggplot2::scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_y_continuous(
      name = ylab,
      sec.axis = ggplot2::sec_axis(
        transform = ~ . - offset_axis,
        name = "Secondary axis for residuals"
      )
    ) +
    ggplot2::labs(
      title = main,
      subtitle = sub,
      x = xlab,
      y = ylab,
      color = ""
    )

  result
}
