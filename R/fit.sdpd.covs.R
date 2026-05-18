#' Compute SDP-D Covariance Components
#'
#' Computes covariance components used by the SDP-D coefficient-estimation
#' routine.
#'
#' The function builds covariance matrices between the observed series, boundary
#' or neighbor-augmented series, and optional regressors. These covariance
#' components are then used by [fit_sdpd_coefficients()].
#'
#' @param series Numeric matrix. Endogenous data series with spatial units in
#'   rows and time observations in columns.
#' @param x Optional numeric matrix or three-dimensional array. Regressor data.
#'   If `kk = 1`, `x` is expected to be a matrix. If `kk > 1`, `x` is expected
#'   to be a three-dimensional array indexed by covariate, spatial unit, and
#'   time.
#' @param kk Integer. Number of covariates. Defaults to `0`.
#' @param px_neighbors Optional list. Pixel or proximity-neighbor information.
#'   If supplied, expected elements are `series_boundary` and `index`.
#' @param nn Integer. Number of time observations.
#' @param pp Integer. Number of spatial units.
#' @param na_covs Character scalar. Missing-value handling method passed to
#'   [stats::cov()]. Defaults to `"pairwise.complete.obs"`.
#'
#' @return A list containing:
#' \describe{
#'   \item{index}{Neighbor index matrix.}
#'   \item{cov11}{Covariance matrix between `series` and the boundary-augmented
#'   series.}
#'   \item{cov12}{Lagged covariance matrix between `series` and the
#'   boundary-augmented series.}
#'   \item{cov_x}{Covariance component for regressors, or `NULL` when no
#'   covariates are used.}
#' }
#'
#' @details
#' If `px_neighbors` is `NULL`, the function uses the original `series` as the
#' boundary series and constructs a full index matrix from the row names of
#' `series`.
#'
#' If `px_neighbors` is supplied, the function appends
#' `px_neighbors$series_boundary` to the original series and uses
#' `px_neighbors$index` as the neighbor index matrix.
#'
#' @examples
#' covs <- fit_sdpd_covs(
#'   series = series,
#'   kk = 0,
#'   nn = ncol(series),
#'   pp = nrow(series)
#' )
#'
#' covs <- fit_sdpd_covs(
#'   series = series,
#'   x = x,
#'   kk = 1,
#'   px_neighbors = px_neighbors,
#'   nn = ncol(series),
#'   pp = nrow(series)
#' )
#'
#' @seealso [fit_sdpd_coefficients()], [stats::cov()]
#'
#' @export
fit_sdpd_covs <- function(series,
                          x = NULL,
                          kk = 0,
                          px_neighbors = NULL,
                          nn,
                          pp,
                          na_covs = "pairwise.complete.obs") {
  if (is.null(px_neighbors)) {
    series_boundary <- series

    index <- matrix(
      rep(dimnames(series)[[1]], pp),
      byrow = TRUE,
      nrow = pp
    )

    dimnames(index)[[1]] <- dimnames(series)[[1]]
  } else {
    series_boundary <- rbind(series, px_neighbors$series_boundary)
    index <- px_neighbors$index
  }

  cov12 <- stats::cov(
    t(series)[2:nn, ],
    t(series_boundary)[-nn, ],
    use = na_covs
  )

  cov11 <- stats::cov(
    t(series),
    t(series_boundary),
    use = na_covs
  )

  if (kk == 0 || is.null(kk)) {
    cov_x <- NULL
  } else if (kk == 1) {
    cov_x <- stats::cov(
      t(x)[2:nn, ],
      t(series_boundary)[-nn, ],
      use = na_covs
    )
  } else if (kk > 1) {
    cov_x <- array(
      0,
      dim = c(dim(x)[1], dim(x)[2], dim(series_boundary)[1])
    )

    dimnames(cov_x)[[1]] <- dimnames(x)[[1]]
    dimnames(cov_x)[[2]] <- dimnames(x)[[2]]
    dimnames(cov_x)[[3]] <- dimnames(series_boundary)[[1]]

    for (jj in seq_len(kk)) {
      cov_x[jj, , ] <- stats::cov(
        t(x[jj, , 2:nn]),
        t(series_boundary)[-nn, ],
        use = na_covs
      )
    }
  }

  list(
    index = index,
    cov11 = cov11,
    cov12 = cov12,
    cov_x = cov_x
  )
}


