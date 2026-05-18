#' Generate an SDP-D Series
#'
#' Simulates a multivariate or spatio-temporal time series from an SDP-D model.
#'
#' The function generates the endogenous series, optional exogenous regressors,
#' and innovation errors using the coefficient and spatial-structure components
#' stored in an SDP-D model object.
#'
#' @param nn Optional integer. Number of time observations to simulate. Required
#'   when `rr_y` is `NULL`.
#' @param rr_y Optional numeric matrix. Initial endogenous series. If supplied
#'   and `nn` is `NULL`, the number of time observations is inferred from
#'   `ncol(rr_y)`.
#' @param rr_xx Optional numeric matrix, data frame, or three-dimensional array.
#'   Exogenous regressor series. If `NULL` and the model has covariates,
#'   regressors are simulated.
#' @param model SDP-D model object, typically created with
#'   [build_sdpd_model()]. It must contain coefficient and spatial-structure
#'   components.
#' @param markovian Logical scalar. Whether to use the Markovian autoregressive
#'   update in [fit_sdpd_series()].
#' @param num_steps Integer. Number of iterative update steps used by
#'   [fit_sdpd_series()].
#' @param covariates_sim_model List. ARIMA model specification passed to
#'   [stats::arima.sim()] when covariates are simulated.
#'
#' @return A list containing:
#' \describe{
#'   \item{rr_y}{Simulated endogenous series.}
#'   \item{rr_xx}{Observed or simulated exogenous regressors.}
#'   \item{eps}{Innovation-error matrix.}
#' }
#'
#' If required inputs are missing or incompatible, a list with an `error` element
#' is returned.
#'
#' @details
#' If `rr_y` is not supplied, the function initializes it as a zero matrix with
#' dimensions inferred from `model$coeffs` and `nn`.
#'
#' If the model contains active beta coefficients and `rr_xx` is not supplied,
#' covariates are simulated as stationary ARIMA processes using
#' [stats::arima.sim()]. A covariate named `"trend"` is generated as a normalized
#' deterministic time trend.
#'
#' @examples
#' \dontrun{
#' simulated_data <- generate_sdpd_series(
#'   nn = 50,
#'   rr_y = NULL,
#'   rr_xx = NULL,
#'   model = model,
#'   markovian = TRUE,
#'   num_steps = 10,
#'   covariates_sim_model = list(ar = c(0.8), sd = 1)
#' )
#'}
#' @seealso [build_sdpd_model()], [check_sdpd_series()],
#'   [build_spatial_matrix()], [fit_sdpd_series()], [stats::arima.sim()]
#'
#' @export
generate_sdpd_series <- function(nn,
                                 rr_y,
                                 rr_xx,
                                 model,
                                 markovian,
                                 num_steps,
                                 covariates_sim_model) {

  # This function simulates a multivariate or spatio-temporal time series from
  # an SDP-D model.

  if (is.null(rr_y) && (is.null(nn) || is.null(model$coeffs))) {
    return(list(
      error = paste(
        "When rr_y is not supplied, simulating the data requires both the",
        "series length and the coefficient matrix."
      )
    ))
  }

  if (is.null(nn)) {
    nn <- dim(rr_y)[2]
  }

  if (is.null(rr_y)) {
    pp <- dim(model$coeffs)[1]

    rr_y <- matrix(0, nrow = pp, ncol = nn)
    dimnames(rr_y)[[1]] <- dimnames(model$coeffs)[[1]]
    dimnames(rr_y)[[2]] <- seq_len(nn)
  } else {
    pp <- dim(rr_y)[1]
  }

  px <- dimnames(rr_y)[[1]]

  kk <- sum(model$beta_coeffs)
  beta_names <- names(model$beta_coeffs)[model$beta_coeffs]

  data <- list(
    pp = pp,
    nn = nn,
    kk = kk
  )

  # Generate errors.
  if (is.null(model$eps)) {
    if (is.null(model$sigma_eps)) {
      model$sigma_eps <- matrix(rep(1, data$pp), ncol = 1)
    } else {
      model$sigma_eps <- as.matrix(model$sigma_eps, ncol = 1)
    }

    series_eps <- t(apply(
      model$sigma_eps,
      1,
      FUN = function(xx, nn) stats::rnorm(nn, sd = xx),
      nn = data$nn
    ))
  } else if (dim(model$eps)[1] != data$pp || dim(model$eps)[2] != data$nn) {
    return(list(
      error = "The dimension of the error matrix in model$eps is not compatible with the data matrix."
    ))
  } else {
    series_eps <- model$eps
  }

  if (is.null(dimnames(series_eps)[[1]])) {
    dimnames(series_eps)[[1]] <- dimnames(rr_y)[[1]]
  }

  # Generate exogenous regressors as stationary ARIMA processes.
  if (data$kk > 0) {
    if (is.null(rr_xx)) {
      if (data$kk == 1) {
        if (beta_names[1] == "trend") {
          xx <- matrix(
            rep(seq_len(data$nn) / data$nn, data$pp),
            byrow = TRUE,
            nrow = data$pp
          )
        } else {
          xx <- matrix(
            stats::arima.sim(
              n = data$nn * data$pp,
              model = covariates_sim_model
            ),
            ncol = data$nn,
            nrow = data$pp,
            byrow = TRUE
          )
        }

        dimnames(xx)[[1]] <- dimnames(rr_y)[[1]]
        dimnames(xx)[[2]] <- dimnames(rr_y)[[2]]
      } else if (data$kk > 1) {
        xx <- array(0, dim = c(data$kk, data$pp, data$nn))

        for (jj in seq_len(data$kk)) {
          if (beta_names[jj] == "trend") {
            xx[jj, , ] <- matrix(
              rep(seq_len(data$nn) / data$nn, data$pp),
              byrow = TRUE,
              nrow = data$pp
            )
          } else {
            xx[jj, , ] <- matrix(
              stats::arima.sim(
                n = data$nn * data$pp,
                model = covariates_sim_model
              ),
              ncol = data$nn,
              nrow = data$pp,
              byrow = TRUE
            )
          }
        }

        dimnames(xx)[[1]] <- beta_names
        dimnames(xx)[[2]] <- dimnames(rr_y)[[1]]
        dimnames(xx)[[3]] <- dimnames(rr_y)[[2]]
      }
    } else if ((is.matrix(rr_xx) || is.data.frame(rr_xx)) &&
               data$kk == 1 &&
               dim(rr_xx)[1] == data$pp &&
               dim(rr_xx)[2] == data$nn &&
               sum(is.na(rr_xx)) == 0) {
      xx <- as.matrix(rr_xx)
    } else if (is.array(rr_xx) &&
               data$kk > 1 &&
               dim(rr_xx)[1] == data$kk &&
               dim(rr_xx)[2] == data$pp &&
               dim(rr_xx)[3] == data$nn &&
               sum(is.na(rr_xx)) == 0) {
      xx <- rr_xx
    } else {
      return(list(error = "Something is wrong with rr_xx."))
    }
  } else {
    xx <- NULL
  }

  # Check validity of data.
  data <- check_sdpd_series(
    series = rr_y,
    xx = xx,
    model = model
  )

  if (length(data$errors) > 0) {
    return(list(
      errors = data$errors,
      warnings = data$warnings
    ))
  }

  # Build spatial weights.
  ww <- build_spatial_matrix(
    ww_index = data$ww_index,
    ww_values = data$ww_values
  )

  rr_y <- fit_sdpd_series(
    data_series = data$series,
    ww = ww,
    x_centered = data$xx,
    num_steps = num_steps,
    model = model,
    markovian = markovian,
    coeff_hat = model$coeffs,
    resids = series_eps
  )$series

  # Results.
  list(
    rr_y = rr_y,
    rr_xx = xx,
    eps = series_eps
  )
}




