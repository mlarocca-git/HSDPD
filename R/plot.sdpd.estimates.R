#' Plot Spatial SDP-D Estimates
#'
#' Creates spatial point plots for one or more estimated SDP-D quantities.
#'
#' The function extracts a matrix-like item from an SDP-D result object and
#' creates one `ggplot` object for each column. Each plot displays the selected
#' estimate at the corresponding longitude and latitude coordinates.
#'
#' @param results List or data frame. SDP-D estimation results containing
#'   longitude coordinates, latitude coordinates, group information, and the
#'   selected `item`.
#' @param item Character scalar. Name of the result component to plot. Defaults
#'   to `"coeff_hat"`.
#' @param main Optional character scalar. Text prepended to each plot title.
#' @param sub Optional character scalar. Text appended to each plot title.
#' @param limits Optional numeric vector of length 2. Limits passed to
#'   [ggplot2::scale_colour_gradient2()].
#' @param mid Numeric scalar or vector. Midpoint value for the diverging color
#'   scale. If scalar, it is recycled across all plotted columns. Defaults to
#'   `0`.
#' @param size_point Numeric scalar. Point size used in the plot. Defaults to
#'   `1`.
#'
#' @return A named list of `ggplot` objects, one for each plotted estimate.
#'
#' @details
#' The function assumes that `results$lon`, `results$lat`, and
#' `results$group$LABEL` are aligned with the rows of `results[[item]]`.
#'
#' For consistency with the renamed package API, the default `item` is
#' `"coeff_hat"`.
#'
#' @examples
#' plots <- plot_sdpd_estimates(
#'   results = fit,
#'   item = "coeff_hat",
#'   main = "Coefficient: "
#' )
#'
#' plots$lambda_0
#'
#' @seealso [ggplot2::ggplot()]
#'
#' @export
plot_sdpd_estimates <- function(results,
                                item = "coeff_hat",
                                main = NULL,
                                sub = NULL,
                                limits = NULL,
                                mid = 0,
                                size_point = 1) {
  plot_list <- list()

  item_matrix <- as.matrix(results[[item]])
  item_names <- dimnames(item_matrix)[[2]]

  if (length(item_names) == 0 && dim(item_matrix)[2] == 1) {
    item_names <- as.character(item)
  }

  if (length(mid) == 1) {
    mid <- rep(mid, length(item_names))
  }

  for (ii in seq_along(item_names)) {
    if (dim(item_matrix)[2] == 1) {
      values <- as.numeric(item_matrix)
    } else {
      values <- item_matrix[, ii]
    }

    plot_data <- data.frame(
      lon = results$lon,
      lat = results$lat,
      group = results$group$LABEL,
      value = values
    )

    plot_list[[ii]] <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = lon, y = lat, group = group, colour = value)
    ) +
      ggplot2::geom_point(size = size_point) +
      ggplot2::scale_colour_gradient2(
        high = "red",
        low = "blue",
        mid = "white",
        midpoint = mid[ii],
        limits = limits
      ) +
      ggplot2::labs(
        x = "",
        y = "",
        title = paste(main, item_names[ii], sub, sep = ""),
        colour = "value"
      )
  }

  names(plot_list) <- item_names

  plot_list
}
