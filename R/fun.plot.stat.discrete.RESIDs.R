#' Plot Spatial Residual Statistics
#'
#' Creates a spatial point plot of residual summary statistics.
#'
#' The function computes a user-defined statistic on model residuals and plots
#' the resulting value at each spatial location. It can also plot residual
#' diagnostic p-values as significant versus not significant.
#'
#' @param results List or data frame. SDP-D estimation results containing
#'   residuals, longitude coordinates, and latitude coordinates.
#' @param statistic Function. Summary function applied to residuals. Defaults to
#'   [mean()].
#' @param main Optional character scalar. Plot title.
#' @param significant_test Logical scalar. Whether values should be converted
#'   into significance classes using `alpha`. Defaults to `FALSE`.
#' @param by_adjusted Logical scalar. Whether p-values should be adjusted using
#'   the Benjamini-Yekutieli method. Defaults to `FALSE`.
#' @param alpha Optional numeric scalar. Significance level used when
#'   `significant_test = TRUE`.
#' @param mid_value Numeric scalar. Midpoint used for the diverging color scale.
#'   Defaults to `0`.
#' @param size_point Numeric scalar. Point size used in the plot. Defaults to
#'   `1`.
#' @param ... Additional arguments passed to `statistic`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' If `results$resid` is a data frame, the statistic is computed row-wise after
#' dropping the first residual column. If `results$resid` is a list, the
#' statistic is applied to each list element. Otherwise, the function assumes a
#' simple vector-like residual object.
#'
#' If `by_adjusted = TRUE`, values are adjusted using [stats::p.adjust()] with
#' `method = "BY"`.
#'
#' @examples
#' \dontrun{
#' plot_stat_discrete_resids(
#'   results = fit,
#'   statistic = mean,
#'   main = "Mean residuals"
#' )
#'
#' plot_stat_discrete_resids(
#'   results = fit,
#'   statistic = fun_lb_test,
#'   significant_test = TRUE,
#'   by_adjusted = TRUE,
#'   alpha = 0.05
#' )
#'}
#' @seealso [ggplot2::ggplot()], [stats::p.adjust()]
#'
#' @export
plot_stat_discrete_resids <- function(results,
                                      statistic = mean,
                                      main = NULL,
                                      significant_test = FALSE,
                                      by_adjusted = FALSE,
                                      alpha = NULL,
                                      mid_value = 0,
                                      size_point = 1,
                                      ...) {
  # Compute summary statistics on residuals.
  if (is.data.frame(results$resid)) {
    plot_data <- data.frame(
      lon = results$lon,
      lat = results$lat,
      value = apply(results$resid[, -1], 1, FUN = statistic, ...)
    )
  } else if (is.list(results$resid)) {
    plot_data <- data.frame(
      lon = results$lon,
      lat = results$lat,
      value = unlist(lapply(results$resid, FUN = statistic, ...))
    )
  } else {
    plot_data <- data.frame(
      lon = results$lon,
      lat = results$lat,
      value = statistic(results$resid, ...)
    )
  }

  # Adjust p-values if requested.
  if (by_adjusted) {
    plot_data$value <- stats::p.adjust(plot_data$value, method = "BY")
  }

  if (significant_test) {
    if (is.null(alpha)) {
      stop("The alpha argument must be supplied when significant_test = TRUE.")
    }

    plot_data$value <- as.factor(
      ifelse(plot_data$value < alpha, "Significant", "Not significant")
    )

    result <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$lon, y = .data$lat, colour = .data$value)
    ) +
      ggplot2::geom_point(size = size_point) +
      ggplot2::guides(fill = "none") +
      ggplot2::labs(
        title = main,
        x = "Longitude",
        y = "Latitude",
        colour = ""
      )
  } else {
    result <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$lon, y = .data$lat, colour = .data$value)
    ) +
      ggplot2::geom_point(size = size_point) +
      ggplot2::guides(fill = "none") +
      ggplot2::scale_colour_gradient2(
        high = "red",
        low = "blue",
        mid = "white",
        midpoint = mid_value
      ) +
      ggplot2::labs(
        title = main,
        x = "Longitude",
        y = "Latitude",
        colour = ""
      )
  }

  result
}

#' Plot SDP-D Residual Statistics
#'
#' User-facing wrapper for [plot_stat_discrete_resids()] using the current
#' package naming style.
#'
#' @param results List or data frame. SDP-D estimation results containing
#'   residuals, longitude coordinates, and latitude coordinates.
#' @param statistic Function. Summary function applied to residuals. Defaults to
#'   [mean()].
#' @param main Optional character scalar. Plot title.
#' @param significant_test Logical scalar. Whether values should be converted
#'   into significance classes using `alpha`. Defaults to `FALSE`.
#' @param by_adjusted Logical scalar. Whether p-values should be adjusted using
#'   the Benjamini-Yekutieli method. Defaults to `FALSE`.
#' @param alpha Optional numeric scalar. Significance level used when
#'   `significant_test = TRUE`.
#' @param mid_value Numeric scalar. Midpoint used for the diverging color scale.
#'   Defaults to `0`.
#' @param size_point Numeric scalar. Point size used in the plot. Defaults to
#'   `1`.
#' @param ... Additional arguments passed to `statistic`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_sdpd_residuals_stat(
#'   results = fit,
#'   statistic = mean,
#'   main = "Mean residuals",
#'   mid_value = 0
#' )
#' }
#' @seealso [plot_stat_discrete_resids()]
#'
#' @export
plot_sdpd_residuals_stat <- function(results,
                                     statistic = mean,
                                     main = NULL,
                                     significant_test = FALSE,
                                     by_adjusted = FALSE,
                                     alpha = NULL,
                                     mid_value = 0,
                                     size_point = 1,
                                     ...) {
  plot_stat_discrete_resids(
    results = results,
    statistic = statistic,
    main = main,
    significant_test = significant_test,
    by_adjusted = by_adjusted,
    alpha = alpha,
    mid_value = mid_value,
    size_point = size_point,
    ...
  )
}
