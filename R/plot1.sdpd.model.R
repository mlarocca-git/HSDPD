#' Plot Detailed SDP-D Model Diagnostics
#'
#' Plots diagnostic panels for selected spatial units and selected model
#' components of a fitted SDP-D model.
#'
#' For each selected spatial unit and model component, the function displays:
#' one scatter plot of the relevant signal against the observed series, and one
#' time-series plot comparing the observed series, component signal, residuals,
#' and optionally fixed effects.
#'
#' @param res_fit Fitted SDP-D model result object. Expected to contain `data`,
#'   `coeff_hat`, `resid`, and `model`.
#' @param n_units Integer, character vector, logical vector, or `"all"`.
#'   Spatial units to plot. Defaults to `"all"`.
#' @param n_vars Integer, character vector, logical vector, or `"all"`.
#'   Model components to plot. Defaults to `"all"`.
#' @param t_axis List. Time-axis options with elements `t_labels` and
#'   `t_points`.
#' @param x_limit Optional numeric vector of length 2. X-axis limits for the
#'   scatter plots.
#' @param y_limit Optional numeric vector of length 2. Y-axis limits for the
#'   scatter plots.
#' @param max_col Integer. Maximum number of component columns shown per unit.
#'   Defaults to `5`.
#' @param point_col Color specification for scatter-plot points. Defaults to
#'   `1`.
#'
#' @return Invisibly returns `NULL`, or a list with an `error` element if invalid
#'   input is detected.
#'
#' @details
#' The function assumes the renamed package API:
#' \itemize{
#'   \item `res_fit$coeff_hat`, not `res_fit$coeff.hat`;
#'   \item `res_fit$data$ww`, not `res_fit$data$W`;
#'   \item `res_fit$data$xx`, not `res_fit$data$X`;
#'   \item lambda columns named `lambda_0`, `lambda_1`, and `lambda_2`.
#' }
#'
#' The `fixed_effects` panel is excluded from the plotted panels because fixed
#' effects are represented as a horizontal dashed grey line in the time-series
#' panel.
#'
#' @examples
#' plot_1_sdpd_model(
#'   res_fit = fit,
#'   n_units = "all",
#'   n_vars = "all"
#' )
#'
#' @export
plot_1_sdpd_model <- function(res_fit,
                              n_units = "all",
                              n_vars = "all",
                              t_axis = list(t_labels = NULL, t_points = NULL),
                              x_limit = NULL,
                              y_limit = NULL,
                              max_col = 5,
                              point_col = 1) {
  # res_fit is expected to contain:
  # coeff_hat, data$series, data$ww, data$xx, fitted, resid, errors, warnings,
  # and model.

  data_series <- res_fit$data$series
  xx <- res_fit$data$xx

  nn <- dim(data_series)[2]
  pp <- dim(data_series)[1]

  alpha_hat <- res_fit$coeff_hat
  error <- NULL

  # Define panels to plot.
  panels <- dimnames(res_fit$coeff_hat)[[2]]
  panel_index <- rep(FALSE, length(panels))

  if (is.logical(n_vars)) {
    n_panels <- min(length(n_vars), length(panel_index))
    panel_index[seq_len(n_panels)] <- n_vars[seq_len(n_panels)]
  } else if (is.numeric(n_vars) && length(n_vars) == 1) {
    n_panels <- min(length(panel_index), as.integer(n_vars))
    panel_index[seq_len(n_panels)] <- TRUE
  } else if (is.numeric(n_vars)) {
    valid_index <- as.integer(n_vars) %in% seq_along(panel_index)
    n_vars <- as.integer(n_vars)[valid_index]

    n_panels <- min(length(n_vars), length(panel_index))
    panel_index[n_vars[seq_len(n_panels)]] <- TRUE
  } else if (n_vars[1] == "all") {
    n_vars <- dimnames(res_fit$coeff_hat)[[2]]
    panel_index <- panels %in% n_vars
  } else if (is.character(n_vars)) {
    panel_index <- panels %in% n_vars
  } else {
    return(list(
      error = "The value passed to n_vars does not follow the expected format."
    ))
  }

  panels <- panels[panel_index]

  if (length(which(panels == "fixed_effects")) > 0) {
    panels <- panels[-which(panels == "fixed_effects")]
  }

  n_col <- min(max_col, length(panels))

  if (n_units[1] == "all") {
    n_units <- seq_len(dim(res_fit$data$series)[1])
  }

  if (is.character(n_units) || is.numeric(n_units)) {
    indices <- seq_len(dim(res_fit$data$series)[1])
    names(indices) <- dimnames(res_fit$data$series)[[1]]
    n_units <- indices[n_units]
  } else {
    return(list(
      error = "The value passed to n_units does not follow the expected format."
    ))
  }

  if (res_fit$model$fixed_effects) {
    fixed_effects_text <- "fixed effect (dashed grey) & "
  } else {
    fixed_effects_text <- ""
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  # Main plotting loop.
  for (ii in n_units) {
    if (sum(is.na(res_fit$coeff_hat[ii, ])) > 0) {
      cat("\nThere are missing estimated coefficients for this series.")
      next
    }

    beta_count <- 0

    graphics::par(
      mfcol = c(2, n_col),
      mai = c(1, 0.8, 0.5, 0.5),
      las = 2
    )

    for (panel in panels[seq_len(n_col)]) {
      x_limit_panel <- x_limit
      y_limit_panel <- y_limit

      if (panel == "lambda_0") {
        panel_index_id <- 1
      } else if (panel == "lambda_1") {
        panel_index_id <- 2
      } else if (panel == "lambda_2") {
        panel_index_id <- 3
      } else {
        panel_index_id <- 4
        beta_count <- beta_count + 1

        if (is.matrix(res_fit$data$xx)) {
          title_expression <- substitute(
            bold(str0 * str1) * str2 * hat(beta)[i] == y,
            list(
              str0 = "i=",
              str1 = dimnames(data_series)[[1]][ii],
              str2 = ": ",
              y = round(alpha_hat[ii, panel], digits = 6)
            )
          )

          xx <- res_fit$data$xx
          x_name <- "X"
        } else if (is.array(res_fit$data$xx)) {
          title_expression <- substitute(
            bold(str0 * str1) * str2 * hat(beta)[i] == y,
            list(
              str0 = "i=",
              str1 = dimnames(data_series)[[1]][ii],
              str2 = ": ",
              y = round(alpha_hat[ii, panel], digits = 6)
            )
          )

          xx <- res_fit$data$xx[panel, , ]
          x_name <- panel
        }
      }

      if (panel_index_id == 4) {
        lag_value <- 0
        weighted_series <- xx

        y_label <- expression(y[i * t])
        x_label <- substitute(
          str0[str1 * str2],
          list(str0 = x_name, str1 = "i", str2 = "t")
        )

        plot_title <- paste(
          "Left axis: observed series (black) & ",
          x_name,
          "'s signal (green)\n Right axis: ",
          fixed_effects_text,
          "residuals (solid grey)",
          sep = ""
        )

        color_offset <- 1
      } else if (panel_index_id == 3) {
        lag_value <- 1
        weighted_series <- res_fit$data$ww %*% data_series

        y_label <- expression(y[i * t])
        x_label <- expression(bold(w)[i] * bold(y)[t - 1])

        plot_title <- paste(
          "Left axis: observed series (black) & spatial-dynamic signal (yellow)\n Right axis: ",
          fixed_effects_text,
          "residuals (solid grey)",
          sep = ""
        )

        title_expression <- substitute(
          bold(str0 * str1) * str2 * hat(lambda)[x * i] == y,
          list(
            str0 = "i=",
            str1 = dimnames(data_series)[[1]][ii],
            str2 = ": ",
            x = 2,
            y = round(alpha_hat[ii, panel], digits = 6)
          )
        )

        color_offset <- 5
      } else if (panel_index_id == 2) {
        lag_value <- 1
        weighted_series <- data_series

        y_label <- expression(y[i * t])
        x_label <- expression(y[i * (t - 1)])

        plot_title <- paste(
          "Left axis: observed series (black) & pure dynamic signal (blue)\n Right axis: ",
          fixed_effects_text,
          "residuals (solid grey)",
          sep = ""
        )

        title_expression <- substitute(
          bold(str0 * str1) * str2 * hat(lambda)[x * i] == y,
          list(
            str0 = "i=",
            str1 = dimnames(data_series)[[1]][ii],
            str2 = ": ",
            x = 1,
            y = round(alpha_hat[ii, panel], digits = 6)
          )
        )

        color_offset <- 2
      } else {
        lag_value <- 0
        weighted_series <- res_fit$data$ww %*% data_series

        y_label <- expression(y[i * t])
        x_label <- expression(bold(w)[i] * bold(y)[t])

        plot_title <- paste(
          "Left axis: observed series (black) & pure spatial signal (red)\n Right axis: ",
          fixed_effects_text,
          "residuals (solid grey)",
          sep = ""
        )

        title_expression <- substitute(
          bold(str0 * str1) * str2 * hat(lambda)[x * i] == y,
          list(
            str0 = "i=",
            str1 = dimnames(data_series)[[1]][ii],
            str2 = ": ",
            x = 0,
            y = round(alpha_hat[ii, panel], digits = 6)
          )
        )

        color_offset <- 0
      }

      # Plot spatial regression.
      if (is.null(x_limit_panel)) {
        x_limit_panel <- range(weighted_series[ii, ], na.rm = TRUE)
      }

      if (is.null(y_limit_panel)) {
        y_limit_panel <- range(data_series[ii, ], na.rm = TRUE)
      }

      graphics::plot(
        weighted_series[ii, seq_len(nn - lag_value)],
        data_series[ii, (1 + lag_value):nn],
        ylim = y_limit_panel,
        xlim = x_limit_panel,
        xlab = "",
        ylab = "",
        cex.lab = 1.5,
        cex.axis = 0.9,
        col = point_col
      )

      graphics::title(
        xlab = x_label,
        ylab = y_label,
        cex.lab = 1.5,
        cex.main = 1.5,
        main = title_expression
      )

      graphics::lines(
        x_limit_panel,
        x_limit_panel * alpha_hat[ii, panel],
        col = 2 + color_offset,
        lwd = 1,
        lty = 1
      )

      # Plot spatial and residual estimated series.
      transform_axis <- function(x, y_limit_1, y_limit_2) {
        y_limit_1[1] +
          (x - y_limit_2[1]) *
            diff(y_limit_1) /
            diff(y_limit_2)
      }

      y_limit_observed <- range(
        data_series[ii, ],
        alpha_hat[ii, panel] * weighted_series[ii, seq_len(nn - lag_value)],
        na.rm = TRUE
      )

      if (res_fit$model$fixed_effects) {
        y_limit_resid <- range(
          res_fit$resid[ii, ],
          res_fit$coeff_hat[ii, "fixed_effects"],
          na.rm = TRUE
        )
      } else {
        y_limit_resid <- range(res_fit$resid[ii, ], na.rm = TRUE)
      }

      if (abs(y_limit_resid[2] - y_limit_observed[2]) <
          diff(y_limit_observed) * 0.1 &&
          abs(y_limit_resid[1] - y_limit_observed[1]) <
          diff(y_limit_observed) * 0.1) {
        y_limit_observed <- range(y_limit_observed, y_limit_resid)
      } else if (y_limit_resid[2] <
                 y_limit_observed[1] - diff(y_limit_observed) * 0.1) {
        y_limit_observed[1] <- y_limit_observed[1] -
          diff(y_limit_observed) * 0.1
      } else if (y_limit_resid[1] >
                 y_limit_observed[2] + diff(y_limit_observed) * 0.1) {
        y_limit_observed[2] <- y_limit_observed[2] +
          diff(y_limit_observed) * 0.1
      }

      stats::ts.plot(
        data_series[ii, (1 + lag_value):nn],
        ylim = y_limit_observed,
        xlab = "",
        ylab = "",
        gpars = list(axes = FALSE, cex.main = 0.9),
        type = "n",
        main = paste(plot_title)
      )

      graphics::lines(
        seq_len(nn),
        transform_axis(
          res_fit$resid[ii, ],
          y_limit_1 = y_limit_resid,
          y_limit_2 = y_limit_resid
        ),
        col = "grey"
      )

      graphics::abline(
        h = transform_axis(
          0,
          y_limit_1 = y_limit_resid,
          y_limit_2 = y_limit_resid
        ),
        col = "grey"
      )

      if (res_fit$model$fixed_effects) {
        graphics::abline(
          h = transform_axis(
            res_fit$coeff_hat[ii, "fixed_effects"],
            y_limit_1 = y_limit_resid,
            y_limit_2 = y_limit_resid
          ),
          col = "grey",
          lty = 2
        )
      }

      graphics::lines(
        seq(1 + lag_value, nn),
        data_series[ii, (1 + lag_value):nn]
      )

      graphics::lines(
        seq(1 + lag_value, nn),
        alpha_hat[ii, panel] * weighted_series[ii, seq_len(nn - lag_value)],
        col = 2 + color_offset,
        lwd = 1
      )

      graphics::axis(2)

      if (!is.null(t_axis$t_points)) {
        if (length(t_axis$t_points) != length(t_axis$t_labels)) {
          error <- c(
            error,
            "The lengths of t_axis$t_points and t_axis$t_labels are different."
          )

          t_axis$t_points <- NULL
        }
      }

      if (is.null(t_axis$t_points)) {
        t_axis$t_points <- seq(1, nn, by = max(1, nn %/% 15))
        t_axis$t_labels <- dimnames(res_fit$data$series)[[2]][t_axis$t_points]
      }

      graphics::axis(
        1,
        at = t_axis$t_points,
        labels = t_axis$t_labels,
        cex.axis = 0.8
      )

      graphics::box()
    }

    cat("\nunit=", ii)

    if (ii != n_units[length(n_units)]) {
      cat("\nClick on the plot window for the next panel.")
      graphics::locator(1)
    } else {
      cat("\nPanels are finished.")
    }
  }

  if (length(error) > 0) {
    return(list(error = error))
  }

  invisible(NULL)
}


