#' Compute Fitted Values and Residuals for an SDP-D Series
#'
#' Computes fitted values and residuals for an SDP-D model, or simulates/update
#' a series using supplied residuals.
#'
#' The function builds the active lambda, beta, fixed-effect, and time-effect
#' components of the SDP-D model. It then computes fitted values in the standard
#' case, or iteratively updates the series when residuals are supplied.
#'
#' @param data_series Numeric matrix. Endogenous data series with spatial units
#'   in rows and time observations in columns.
#' @param ww Numeric matrix. Spatial weight matrix.
#' @param x_centered Optional numeric matrix or three-dimensional array.
#'   Mean-centered regressors. Required when the model has covariates.
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()].
#' @param coeff_hat Numeric matrix. Estimated coefficients.
#' @param time_effects Optional numeric vector of time effects.
#' @param px_sim Optional character or numeric vector identifying the spatial
#'   units to update during simulation.
#' @param resids Optional numeric matrix of residuals. If supplied, the function
#'   computes an iterative simulated/update path.
#' @param markovian Logical scalar. Whether to use the Markovian autoregressive
#'   update when `resids` are supplied. Defaults to `TRUE`.
#' @param num_steps Integer. Number of iterative update steps when `resids` are
#'   supplied. Defaults to `10`.
#'
#' @return A list containing:
#' \describe{
#'   \item{series}{Possibly updated endogenous series.}
#'   \item{fitted}{Fitted values.}
#'   \item{resid}{Residuals.}
#' }
#'
#' @details
#' If `resids` is `NULL`, the function computes standard fitted values and
#' residuals from the supplied data series.
#'
#' If `resids` is supplied and `markovian = TRUE`, the function uses a Markovian
#' autoregressive update based on the inverted contemporaneous spatial component.
#'
#' If `resids` is supplied and `markovian = FALSE`, the function uses a
#' fixed-point autoregressive update.
#'
#' @examples
#' \dontrun{
#' fitted_series <- fit_sdpd_series(
#'   data_series = data_series,
#'   ww = ww,
#'   x_centered = x_centered,
#'   model = model,
#'   coeff_hat = coeff_hat
#' )
#'
#' simulated_series <- fit_sdpd_series(
#'   data_series = data_series,
#'   ww = ww,
#'   x_centered = x_centered,
#'   model = model,
#'   coeff_hat = coeff_hat,
#'   resids = resids,
#'   markovian = TRUE
#' )
#'}
#' @seealso [build_sdpd_model()]
#'
#' @export
fit_sdpd_series <- function(data_series,
                            ww,
                            x_centered = NULL,
                            model,
                            coeff_hat,
                            time_effects = NULL,
                            px_sim = NULL,
                            resids = NULL,
                            markovian = TRUE,
                            num_steps = 10) {
  beta_names <- names(model$beta_coeffs)[model$beta_coeffs]

  px <- dimnames(data_series)[[1]]
  time_index <- dimnames(data_series)[[2]]

  data <- list(
    pp = length(px),
    nn = length(time_index),
    kk = length(beta_names)
  )

  # Build model components.
  lambda_0_matrix <- matrix(0, nrow = data$pp, ncol = data$pp)
  lambda_1_matrix <- matrix(0, nrow = data$pp, ncol = data$pp)
  lambda_2_matrix <- matrix(0, nrow = data$pp, ncol = data$pp)

  identity_matrix <- diag(rep(1, data$pp))

  if (model$lambda_coeffs["lambda_0"]) {
    lambda_0_matrix <- diag(coeff_hat[, "lambda_0"])
  }

  if (model$lambda_coeffs["lambda_1"]) {
    lambda_1_matrix <- diag(coeff_hat[, "lambda_1"])
  }

  if (model$lambda_coeffs["lambda_2"]) {
    lambda_2_matrix <- diag(coeff_hat[, "lambda_2"])
  }

  beta_component <- matrix(0, nrow = data$pp, ncol = data$nn)
  fixed_component <- matrix(0, nrow = data$pp, ncol = data$nn)

  if (is.null(model$time_effects) || !model$time_effects || is.null(time_effects)) {
    time_effects <- numeric(data$nn)
  }

  time_component <- matrix(
    rep(time_effects, data$pp),
    byrow = TRUE,
    nrow = data$pp
  )

  if (model$fixed_effects) {
    fixed_component <- coeff_hat[, "fixed_effects"] %*% t(rep(1, data$nn))
  }

  if (data$kk == 1) {
    beta_component <- diag(coeff_hat[, beta_names[1]]) %*% x_centered
  } else if (data$kk > 1) {
    for (jj in seq_len(data$kk)) {
      beta_component <- beta_component +
        diag(coeff_hat[, beta_names[jj]]) %*% x_centered[jj, , ]
    }
  }

  if (is.null(px_sim)) {
    px_sim <- as.character(px)
  } else {
    px_sim <- as.character(px_sim)
  }

  if (is.null(dimnames(data_series))) {
    dimnames(data_series)[[1]] <- seq_len(dim(data_series)[1])
  }

  if (!is.null(resids) && is.null(dimnames(resids))) {
    dimnames(resids)[[1]] <- dimnames(data_series)[[1]]
  }

  # Estimation results.
  if (!is.null(resids) && markovian) {
    # Markovian autoregressive case.
    fitted_values <- new_resid <- as.matrix(resids)

    inverse_matrix <- solve(identity_matrix - lambda_0_matrix %*% ww)

    transition_matrix <- inverse_matrix %*%
      (lambda_1_matrix + lambda_2_matrix %*% ww)

    innovation_component <- inverse_matrix %*%
      (beta_component + fixed_component + new_resid)

    for (step in seq_len(num_steps)) {
      if (step > 1) {
        data_series[, 1] <- fitted_values[, 2] + innovation_component[, 1]
      }

      for (tt in 2:data$nn) {
        fitted_values[, tt] <- transition_matrix %*% data_series[, tt - 1]

        data_series[px_sim, tt] <- fitted_values[px_sim, tt] +
          innovation_component[px_sim, tt]
      }
    }

    fitted_values <- data_series - new_resid
  } else if (!is.null(resids)) {
    # Fixed-point autoregressive case.
    fitted_values <- new_resid <- as.matrix(resids)

    contemporaneous_matrix <- lambda_0_matrix %*% ww
    lagged_matrix <- lambda_1_matrix + lambda_2_matrix %*% ww

    innovation_component <- beta_component +
      fixed_component +
      time_component +
      new_resid

    dimnames(innovation_component) <- dimnames(data_series)

    for (step in seq_len(num_steps)) {
      if (step > 1) {
        data_series[, 1] <- fitted_values[, 2] + innovation_component[, 1]
      }

      for (tt in 2:data$nn) {
        fitted_values[, tt] <- contemporaneous_matrix %*% data_series[, tt] +
          lagged_matrix %*% data_series[, tt - 1]

        data_series[px_sim, tt] <- fitted_values[px_sim, tt] +
          innovation_component[px_sim, tt]
      }
    }

    fitted_values <- data_series - new_resid
  } else {
    # Standard fitted-value computation.
    fitted_values <- as.matrix(data_series)

    fitted_values[, 2:data$nn] <-
      lambda_0_matrix %*% ww %*% data_series[, 2:data$nn] +
      (lambda_1_matrix + lambda_2_matrix %*% ww) %*%
        data_series[, -data$nn] +
      beta_component[, -1] +
      fixed_component[, -1] +
      time_component[, -1]

    fitted_values[, 1] <- rep(NA, data$pp)

    new_resid <- data_series - fitted_values
  }

  # Return.
  list(
    series = data_series,
    fitted = fitted_values,
    resid = new_resid
  )
}

