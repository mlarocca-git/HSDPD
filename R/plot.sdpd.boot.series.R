#' Plot a Bootstrap SDP-D Series
#'
#' Plots one bootstrapped SDP-D time series against the corresponding observed
#' series for diagnostic control.
#'
#' The function randomly selects one spatial unit from the first bootstrap group
#' and saves a JPEG plot comparing the observed series and one bootstrapped
#' series. The plot title reports selected model diagnostics, lambda estimates,
#' neighborhood size, eigenvalue information, and bootstrap deviations.
#'
#' @param bb Integer. Current bootstrap iteration.
#' @param pp Integer. Number of spatial units.
#' @param nn Integer. Number of time observations.
#' @param eigen Numeric vector. Eigenvalue diagnostics.
#' @param yy_star_1 Numeric matrix. Bootstrapped time series.
#' @param boot_options List. Bootstrap options. Expected elements include
#'   `y_limits`, `label_index`, and `folder`.
#' @param res_fit List. Fitted SDP-D model result containing data, coefficients,
#'   and model information.
#' @param group_index Vector. Group index used to select the plotted spatial
#'   unit.
#'
#' @return Invisibly returns the output file path of the generated JPEG plot.
#'
#' @details
#' If `boot_options$y_limits` is `NULL`, the y-axis limits are inferred from
#' `res_fit$data$series`.
#'
#' The function assumes that `res_fit$data` contains `series`, `ww`,
#' `px_neighbors`, and optionally `xx`. It also assumes that `res_fit$coeff_hat`
#' follows the renamed coefficient convention, with columns such as `lambda_0`,
#' `lambda_1`, `lambda_2`, and `fixed_effects`.
#'
#' @examples
#' plot_sdpd_boot_series(
#'   bb = 500,
#'   pp = pp,
#'   nn = nn,
#'   eigen = eigen_values,
#'   yy_star_1 = yy_star_1,
#'   boot_options = boot_options,
#'   res_fit = res_fit,
#'   group_index = group_index
#' )
#'
#' @export
plot_sdpd_boot_series <- function(bb,
                                  pp,
                                  nn,
                                  eigen,
                                  yy_star_1,
                                  boot_options,
                                  res_fit,
                                  group_index) {
  if (is.null(boot_options$y_limits)) {
    boot_options$y_limits <- range(res_fit$data$series, na.rm = TRUE)
  }

  cat(
    "\nBootstrap iterations:",
    bb - 300 + 1,
    "-",
    bb,
    "... one of the bootstrapped time series is plotted for control."
  )

  location_index <- sample(
    seq_len(pp)[group_index == 1],
    size = 1
  )

  title_text <- "True series (red) vs bootstrap series (black)"
  title_suffix <- NULL

  title_text <- paste(
    title_text,
    " - location=",
    dimnames(res_fit$data$series)[[1]][location_index],
    sep = ""
  )

  if (!is.null(boot_options$label_index)) {
    title_text <- paste(
      title_text,
      " - group=",
      boot_options$label_index[location_index],
      sep = ""
    )
  }

  title_text <- paste(
    title_text,
    "\n core.points=",
    sum(abs(res_fit$data$ww[location_index, ]) > 0),
    ", n.points=",
    sum(
      res_fit$data$px_neighbors$index[location_index, ] > 0,
      na.rm = TRUE
    ),
    sep = ""
  )

  title_text <- paste(title_text, "  (", sep = "")

  lambda_sum <- 0

  if ("lambda_0" %in% dimnames(res_fit$coeff_hat)[[2]]) {
    title_text <- paste(
      title_text,
      "l0=",
      round(res_fit$coeff_hat[location_index, "lambda_0"], digits = 2),
      sep = ""
    )

    lambda_sum <- lambda_sum + res_fit$coeff_hat[location_index, "lambda_0"]
  }

  if ("lambda_1" %in% dimnames(res_fit$coeff_hat)[[2]]) {
    title_text <- paste(
      title_text,
      ", l1=",
      round(res_fit$coeff_hat[location_index, "lambda_1"], digits = 2),
      sep = ""
    )

    lambda_sum <- lambda_sum + res_fit$coeff_hat[location_index, "lambda_1"]
  }

  if ("lambda_2" %in% dimnames(res_fit$coeff_hat)[[2]]) {
    title_text <- paste(
      title_text,
      ", l2=",
      round(res_fit$coeff_hat[location_index, "lambda_2"], digits = 2),
      sep = ""
    )

    lambda_sum <- lambda_sum + res_fit$coeff_hat[location_index, "lambda_2"]
  }

  title_text <- paste(
    title_text,
    title_suffix,
    ") eigen=",
    round(eigen[1], digits = 2),
    sep = ""
  )

  if ("lambda_1" %in% dimnames(res_fit$coeff_hat)[[2]]) {
    title_text <- paste(
      title_text,
      ", #{|l1|>1}=",
      sum(ifelse(abs(res_fit$coeff_hat[, "lambda_1"]) > 1, 1, 0)),
      sep = ""
    )
  }

  title_text <- paste(
    title_text,
    ", #{|l0+l1+l2|>1}=",
    as.integer(abs(lambda_sum) > 1),
    sep = ""
  )

  title_text <- paste(
    title_text,
    ", DEVS=(",
    round(mean(yy_star_1[location_index, ] - res_fit$data$series[location_index, ]), 2),
    ";",
    round(stats::sd(yy_star_1[location_index, ] - res_fit$data$series[location_index, ]), 2),
    ")",
    sep = ""
  )

  output_file <- paste(
    boot_options$folder,
    boot_options$label_index[location_index],
    dimnames(res_fit$data$series)[[1]][location_index],
    "bootplot_",
    bb,
    ".jpeg",
    sep = ""
  )

  subtitle_text <- paste(
    res_fit$model$name_endogenous,
    "~H-SDPD  (based on data from ",
    dimnames(res_fit$data$series)[[2]][1],
    " to ",
    dimnames(res_fit$data$series)[[2]][2],
    ")",
    sep = ""
  )

  grDevices::jpeg(
    file = output_file,
    width = 1100,
    height = 600
  )

  stats::ts.plot(
    yy_star_1[location_index, ],
    col = group_index[location_index],
    ylim = boot_options$y_limits,
    ylab = "values"
  )

  graphics::title(
    main = title_text,
    sub = subtitle_text,
    cex.main = 1.5
  )

  graphics::lines(
    seq_len(nn),
    res_fit$data$series[location_index, ],
    col = "red"
  )

  if ("fixed_effects" %in% dimnames(res_fit$coeff_hat)[[2]]) {
    graphics::abline(
      h = res_fit$coeff_hat[location_index, "fixed_effects"]
    )
  }

  kk <- res_fit$data$kk

  if (!is.null(res_fit$data$xx) && !is.null(kk)) {
    if (kk == 1) {
      graphics::lines(
        seq_len(nn),
        res_fit$data$xx[location_index, ],
        col = "blue",
        lty = 2
      )
    } else if (kk > 1) {
      graphics::lines(
        seq_len(nn),
        res_fit$data$xx[1, location_index, ],
        col = "blue",
        lty = 2
      )
    }
  }

  grDevices::dev.off()

  invisible(output_file)
}

