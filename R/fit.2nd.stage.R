.fit_second_stage_from_design <- function(design, coeff_hat) {
  fit_second_stage(
    data_series = design$series,
    x_centered = design$x,
    ww = design$ww,
    model = design$model,
    coeff_hat = coeff_hat
  )
}

#' Fit the Second-Stage SDP-D Estimator
#'
#' Fits the second-stage estimator for the beta coefficients and, when included
#' in the model, the fixed effects.
#'
#' The function uses the dynamic spatial model components and an existing matrix
#' of coefficient estimates to update the beta coefficients through
#' location-specific least-squares regressions.
#'
#' @param data_series Numeric matrix. Endogenous data series with spatial units
#'   in rows and time observations in columns.
#' @param x_centered Numeric matrix or array. Mean-centered regressors. If the
#'   model has more than one covariate, this is expected to be a
#'   three-dimensional array indexed by covariate, spatial unit, and time.
#' @param ww Numeric matrix. Spatial weight matrix.
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()].
#' @param coeff_hat Numeric matrix. Current coefficient estimates to be updated.
#' @param data Optional list. Data summary object. Recomputed internally from
#'   `data_series` and the active beta coefficients.
#'
#' @return A numeric matrix of updated coefficient estimates.
#'
#' @details
#' The function constructs the lambda components of the reduced-form model and
#' computes the transformed response used in the second-stage regression.
#'
#' If `model$fixed_effects = TRUE`, the local regression includes an intercept
#' term and updates both beta coefficients and `fixed_effects`.
#'
#' @examples
#' \dontrun{
#' coeff_hat <- fit_second_stage(
#'   data_series = data_series,
#'   x_centered = x_centered,
#'   ww = ww,
#'   model = model,
#'   coeff_hat = coeff_hat
#' )
#'}
#' @seealso [build_sdpd_model()]
#'
#' @export
fit_second_stage <- function(data_series,
                             x_centered,
                             ww,
                             model,
                             coeff_hat,
                             data = NULL) {

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

  if (model$lambda_coeffs["lambda_0"]) {
    lambda_0_matrix <- diag(coeff_hat[, "lambda_0"])
  }

  if (model$lambda_coeffs["lambda_1"]) {
    lambda_1_matrix <- diag(coeff_hat[, "lambda_1"])
  }

  if (model$lambda_coeffs["lambda_2"]) {
    lambda_2_matrix <- diag(coeff_hat[, "lambda_2"])
  }

  # Second-stage estimator for beta coefficients.
  y_second_stage <- data_series[, 2:data$nn] -
    lambda_0_matrix %*% ww %*% data_series[, 2:data$nn] -
    lambda_1_matrix %*% data_series[, 1:(data$nn - 1)] -
    lambda_2_matrix %*% ww %*% data_series[, 1:(data$nn - 1)]

  if (model$fixed_effects) {
    if (data$kk > 1) {
      for (ii in seq_len(data$pp)) {
        x_i <- cbind(
          t(x_centered[beta_names, ii, -1]),
          rep(1, data$nn - 1)
        )

        coeff_hat[ii, c(beta_names, "fixed_effects")] <-
          solve(t(x_i) %*% x_i) %*% t(x_i) %*% y_second_stage[ii, ]
      }
    } else {
      for (ii in seq_len(data$pp)) {
        x_i <- cbind(
          t(x_centered[ii, -1]),
          rep(1, data$nn - 1)
        )

        coeff_hat[ii, c(beta_names, "fixed_effects")] <-
          solve(t(x_i) %*% x_i) %*% t(x_i) %*% y_second_stage[ii, ]
      }
    }
  } else if (data$kk > 1) {
    for (ii in seq_len(data$pp)) {
      x_i <- t(x_centered[beta_names, ii, -1])

      coeff_hat[ii, beta_names] <-
        solve(t(x_i) %*% x_i) %*% t(x_i) %*% y_second_stage[ii, ]
    }
  } else if (data$kk == 1) {
    for (ii in seq_len(data$pp)) {
      x_i <- t(x_centered[ii, -1])

      coeff_hat[ii, beta_names] <-
        solve(t(x_i) %*% x_i) %*% t(x_i) %*% y_second_stage[ii, ]
    }
  }

  coeff_hat
}
