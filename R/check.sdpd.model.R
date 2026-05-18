#' Check an SDP-D Model
#'
#' Checks the validity of an SDP-D model or fitted SDP-D model by evaluating
#' stationarity-related diagnostics.
#'
#' The function computes eigenvalue-based diagnostics for the reduced-form model
#' matrix using the coefficient matrix and the spatial weight matrix.
#'
#' @param res_fit Optional fitted SDP-D model object. If supplied, coefficients
#'   are taken from `res_fit$coeff_hat` and the model structure is taken from
#'   `res_fit$model`.
#' @param model Optional SDP-D model object, typically created with
#'   [build_sdpd_model()]. Ignored if `res_fit` is supplied.
#' @param ww_index Optional matrix. Spatial-neighbor index matrix. If omitted,
#'   it is taken from `model$ww_index`.
#' @param ww_values Optional matrix. Spatial-weight value matrix. If omitted,
#'   it is taken from `model$ww_values`.
#'
#' @return A list containing:
#' \describe{
#'   \item{diagnostics}{A data frame with eigenvalue-based diagnostics:
#'   `max_mod_eigen_1`, `max_mod_eigen_2`, and `max_mod_eigen_a`.}
#' }
#'
#' If the spatial matrix components are missing, a list with an `error` element
#' is returned.
#'
#' @details
#' If `res_fit` is supplied, the function evaluates the model inherited from
#' `res_fit`; any separately supplied `model` object is ignored.
#'
#' The function builds the spatial matrix using [build_spatial_matrix()] and
#' then computes diagnostics based on the lambda components available in the
#' coefficient matrix.
#'
#' @examples
#' diagnostics <- check_sdpd_model(model = model)
#'
#' diagnostics <- check_sdpd_model(
#'   model = model,
#'   ww_index = model$ww_index,
#'   ww_values = model$ww_values
#' )
#'
#' @seealso [build_sdpd_model()], [build_spatial_matrix()]
#'
#' @export
check_sdpd_model <- function(res_fit = NULL,
                             model = NULL,
                             ww_index = NULL,
                             ww_values = NULL) {

  ## This function checks the validity of a model or fitted model by evaluating
  ## stationarity conditions.

  if (!is.null(res_fit)) {
    coeffs <- res_fit$coeff_hat

    if (!is.null(model)) {
      cat("\nWarning: the model inherited from res_fit has been evaluated; the supplied model object has been ignored.")
    }

    model <- res_fit$model
  } else {
    coeffs <- NULL
  }

  if (!is.null(model)) {
    if (is.null(coeffs)) {
      coeffs <- model$coeffs
    }

    px <- dimnames(coeffs)[[1]]
  }

  # Check validity of the spatial matrix.
  if (is.null(ww_index) || is.null(ww_values)) {
    ww_index <- model$ww_index[px, ]
    ww_values <- model$ww_values[px, ]
  }

  if (is.null(ww_index) || is.null(ww_values)) {
    return(list(error = "The spatial matrix components are missing."))
  }

  ww <- build_spatial_matrix(
    ww_index = ww_index,
    ww_values = ww_values
  )

  if (!is.null(coeffs)) {
    pp <- dim(coeffs)[1]

    # Compute eigenvalues of the reduced-form matrix components.
    lambda_0_matrix <- matrix(0, ncol = pp, nrow = pp)
    lambda_1_matrix <- matrix(0, ncol = pp, nrow = pp)
    lambda_2_matrix <- matrix(0, ncol = pp, nrow = pp)

    eigen_values_1 <- rep(0, pp)
    eigen_values_2 <- rep(0, pp)
    eigen_values_a <- rep(0, pp)

    if ("lambda_0" %in% dimnames(coeffs)[[2]]) {
      lambda_0_matrix <- diag(coeffs[, "lambda_0"]) %*% ww
      matrix_1 <- solve(diag(rep(1, pp)) - lambda_0_matrix)

      if (sum(is.na(matrix_1)) == 0) {
        eigen_values_1 <- eigen(matrix_1)$values
      } else {
        eigen_values_1 <- rep(NA, pp)
      }
    } else {
      matrix_1 <- diag(eigen_values_1)
    }

    if ("lambda_1" %in% dimnames(coeffs)[[2]]) {
      lambda_1_matrix <- diag(coeffs[, "lambda_1"])
    }

    if ("lambda_2" %in% dimnames(coeffs)[[2]]) {
      lambda_2_matrix <- diag(coeffs[, "lambda_2"]) %*% ww
    }

    matrix_2 <- lambda_1_matrix + lambda_2_matrix

    if (sum(is.na(matrix_2)) == 0) {
      eigen_values_2 <- eigen(matrix_2)$values
    } else {
      eigen_values_2 <- rep(NA, pp)
    }

    if (sum(is.na(matrix_1 %*% matrix_2)) == 0) {
      eigen_values_a <- eigen(matrix_1 %*% matrix_2)$values
    } else {
      eigen_values_a <- rep(NA, pp)
    }

    diagnostics <- cbind(
      max_mod_eigen_1 = Mod(eigen_values_1),
      max_mod_eigen_2 = Mod(eigen_values_2),
      max_mod_eigen_a = Mod(eigen_values_a)
    )
  } else {
    cat("\nNo coefficients to evaluate.")
    diagnostics <- NULL
  }

  # Output.
  list(diagnostics = data.frame(diagnostics))
}

