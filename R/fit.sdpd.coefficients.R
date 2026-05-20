.fit_sdpd_coefficients_from_design <- function(design, covs) {
  fit_sdpd_coefficients(
    ww = design$ww,
    covs = covs,
    mu = design$mu,
    model = design$model
  )
}

#' Estimate SDP-D Coefficients
#'
#' Estimates the first-stage SDP-D model coefficients using spatial covariance
#' components and, when included in the model, computes fixed effects.
#'
#' This is a low-level estimation helper. Most users should call
#' [fit_sdpd_model()] on a series object instead of calling this function
#' directly.
#'
#' The function estimates active lambda coefficients and beta coefficients for
#' each spatial unit. If fixed effects are included in the model, they are
#' computed from the reduced-form coefficient matrix and the spatial-unit means.
#'
#' @param ww Numeric matrix. Spatial weight matrix.
#' @param covs List. Covariance components used for coefficient estimation.
#'   Expected elements are `index`, `cov12`, `cov11`, and, when covariates are
#'   included, `cov_x`.
#' @param mu Numeric vector. Spatial-unit means used to compute fixed effects.
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()].
#'
#' @return A list containing:
#' \describe{
#'   \item{coeff_hat}{Numeric matrix of estimated coefficients.}
#' }
#'
#' @details
#' Coefficients are estimated separately for each spatial unit through local
#' least-squares systems based on the covariance components.
#'
#' If the local normal-equation matrix is not invertible, the corresponding
#' coefficient estimates are set to `NA`.
#'
#' @examples
#' \dontrun{
#' coeffs <- fit_sdpd_coefficients(
#'   ww = ww,
#'   covs = covs,
#'   mu = mu,
#'   model = model
#' )
#'}
#' @seealso [build_sdpd_model()]
#'
#' @export
fit_sdpd_coefficients <- function(ww, covs, mu, model) {
  lambda_names <- names(model$lambda_coeffs)[model$lambda_coeffs]
  beta_names <- names(model$beta_coeffs)[model$beta_coeffs]
  fixed_effects_name <- "fixed_effects"[model$fixed_effects]

  px <- dimnames(ww)[[1]]

  data <- list(
    ww = ww,
    pp = length(px),
    kk = length(beta_names)
  )

  # Variable definition.
  coeff_hat <- matrix(
    0,
    nrow = data$pp,
    ncol = sum(model$fixed_effects) +
      sum(model$lambda_coeffs) +
      sum(model$beta_coeffs)
  )

  dimnames(coeff_hat)[[2]] <- c(
    lambda_names,
    beta_names,
    fixed_effects_name
  )

  dimnames(coeff_hat)[[1]] <- px

  unit_vector <- numeric(data$pp)
  names(unit_vector) <- px

  # Estimation of coefficients.
  is_invertible <- function(matrix_obj) {
    inherits(try(solve(matrix_obj), silent = TRUE), "matrix")
  }

  if (data$kk == 0) {
    coeff_names <- lambda_names

    for (ii in px) {
      weight_i <- ww[ii, px]

      unit_vector[px] <- 0
      unit_vector[ii] <- 1

      indices <- as.character(stats::na.exclude(covs$index[ii, ]))

      y_i <- t(covs$cov12[px, indices]) %*% unit_vector

      x_i <- cbind(
        t(covs$cov12[px, indices]) %*% weight_i,
        t(covs$cov11[px, indices]) %*% unit_vector,
        t(covs$cov11[px, indices]) %*% weight_i
      )[, model$lambda_coeffs]

      if (!is_invertible(t(x_i) %*% x_i)) {
        coeff_hat[ii, coeff_names] <- NA
        next
      }

      coeff_hat[ii, coeff_names] <-
        solve(t(x_i) %*% x_i) %*% t(x_i) %*% y_i
    }
  } else if (data$kk == 1) {
    coeff_names <- c(lambda_names, beta_names)

    for (ii in px) {
      weight_i <- ww[ii, px]

      unit_vector[px] <- 0
      unit_vector[ii] <- 1

      indices <- as.character(stats::na.exclude(covs$index[ii, ]))

      y_i <- t(covs$cov12[px, indices]) %*% unit_vector

      x_i <- cbind(
        t(covs$cov12[px, indices]) %*% weight_i,
        t(covs$cov11[px, indices]) %*% unit_vector,
        t(covs$cov11[px, indices]) %*% weight_i,
        t(covs$cov_x[px, indices]) %*% unit_vector
      )[, c(model$lambda_coeffs, TRUE)]

      if (!is_invertible(t(x_i) %*% x_i)) {
        coeff_hat[ii, coeff_names] <- NA
        next
      }

      coeff_hat[ii, coeff_names] <-
        solve(t(x_i) %*% x_i) %*% t(x_i) %*% y_i
    }
  } else if (data$kk > 1) {
    coeff_names <- c(lambda_names, beta_names)

    for (ii in px) {
      weight_i <- ww[ii, px]

      unit_vector[px] <- 0
      unit_vector[ii] <- 1

      indices <- as.character(stats::na.exclude(covs$index[ii, ]))

      y_i <- t(covs$cov12[px, indices]) %*% unit_vector

      x_i <- cbind(
        t(covs$cov12[px, indices]) %*% weight_i,
        t(covs$cov11[px, indices]) %*% unit_vector,
        t(covs$cov11[px, indices]) %*% weight_i
      )[, model$lambda_coeffs]

      x_i <- cbind(
        x_i,
        t(covs$cov_x[beta_names, ii, indices])
      )

      if (!is_invertible(t(x_i) %*% x_i)) {
        coeff_hat[ii, coeff_names] <- NA
        next
      }

      coeff_hat[ii, coeff_names] <-
        solve(t(x_i) %*% x_i) %*% t(x_i) %*% y_i
    }
  }

  if (model$fixed_effects) {
    lambda_0_matrix <- matrix(0, ncol = data$pp, nrow = data$pp)
    lambda_1_matrix <- matrix(0, ncol = data$pp, nrow = data$pp)
    lambda_2_matrix <- matrix(0, ncol = data$pp, nrow = data$pp)

    if (model$lambda_coeffs[1]) {
      lambda_0_matrix <- diag(coeff_hat[px, "lambda_0"]) %*% ww
    }

    if (model$lambda_coeffs[2]) {
      lambda_1_matrix <- diag(coeff_hat[px, "lambda_1"])
    }

    if (model$lambda_coeffs[3]) {
      lambda_2_matrix <- diag(coeff_hat[px, "lambda_2"]) %*% ww
    }

    b_matrix <- diag(rep(1, data$pp)) -
      lambda_0_matrix -
      lambda_1_matrix -
      lambda_2_matrix

    coeff_hat[px, fixed_effects_name] <- b_matrix %*% mu[px]
  }

  # Estimation results.
  list(coeff_hat = coeff_hat)
}

