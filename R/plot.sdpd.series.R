#' Plot an SDP-D Series at a Spatial Location
#'
#' Plots a time series for the spatial location matching the supplied latitude
#' and longitude after coordinate rounding.
#'
#' @param rr_y Data frame. Series data containing `latitude`, `longitude`, and
#'   time-indexed series columns.
#' @param lat Numeric scalar. Latitude of the location to plot.
#' @param lon Numeric scalar. Longitude of the location to plot.
#' @param n_digits Integer. Number of digits used to round coordinates before
#'   matching. Defaults to `15`.
#' @param main Optional character scalar. Plot title.
#' @param sub Optional character scalar. Plot subtitle.
#' @param xlab Optional character scalar. X-axis label.
#' @param ylab Optional character scalar. Y-axis label.
#'
#' @return A `ggplot` object, or a character message if the requested
#'   coordinates are not present in the data.
#'
#' @details
#' The function expects `rr_y` to contain coordinate columns named `latitude`
#' and `longitude`. All other columns are interpreted as time-indexed series
#' values after the selected location has been filtered.
#'
#' Numeric time labels are converted to dates using [as.Date()].
#'
#' @examples
#' \dontrun{
#' plot_sdpd_series(
#'   rr_y = rr_y,
#'   lat = 45.1,
#'   lon = 9.2
#' )
#'}
#' @seealso [ggplot2::ggplot()]
#'
#' @export
plot_sdpd_series <- function(rr_y,
                             lat,
                             lon,
                             n_digits = 15,
                             main = NULL,
                             sub = NULL,
                             xlab = NULL,
                             ylab = NULL) {
  selected_series <- rr_y |>
    dplyr::mutate(latitude_rounded = round(.data$latitude, n_digits)) |>
    dplyr::mutate(longitude_rounded = round(.data$longitude, n_digits)) |>
    dplyr::filter(
      .data$latitude_rounded == round(lat, n_digits),
      .data$longitude_rounded == round(lon, n_digits)
    ) |>
    dplyr::select(
      -dplyr::any_of(c(
        "longitude_rounded",
        "latitude_rounded",
        "longitude",
        "latitude"
      ))
    )

  if (nrow(selected_series) == 0) {
    return("These coordinates are not present in the database.")
  }

  time_index <- dimnames(selected_series)[[2]]

  if (inherits(time_index, "Date")) {
    time_range <- as.Date(time_index)
  } else if (all(!is.na(as.numeric(time_index)))) {
    time_range <- as.Date(as.numeric(time_index))
  } else {
    time_range <- time_index
  }

  plot_data <- data.frame(
    time = time_range,
    series = as.numeric(selected_series)
  )

  result <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = .data$time)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$series, color = "Series")) +
    ggplot2::guides(x = ggplot2::guide_axis(angle = 0)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_y_continuous(name = ylab) +
    ggplot2::labs(
      title = main,
      subtitle = sub,
      x = xlab,
      y = ylab,
      color = ""
    )

  if (inherits(time_range, "Date")) {
    result <- result +
      ggplot2::scale_x_date(date_labels = "%Y", date_breaks = "1 year")
  }

  result
}


