#' Assemble SDP-D Estimation Results
#'
#' Assembles the output of a single spatial-time SDP-D estimation into a
#' standardized result object.
#'
#' The function extracts fitted values, residuals, coefficient estimates,
#' residual diagnostics, eigenvalue diagnostics, spatial coordinates, group
#' information, and mean-equation components from an estimation object.
#'
#' @param obj_stime List. Single spatial-time estimation object containing
#'   fitted values, residuals, coefficient estimates, diagnostics, spatial
#'   metadata, and mean-equation components.
#'
#' @return A list containing the assembled estimation results, or `NULL` if
#'   `obj_stime` is `NULL` or contains an estimation error.
#'
#' The returned list contains:
#' \describe{
#'   \item{px}{Named spatial-unit index.}
#'   \item{lon}{Longitude coordinates.}
#'   \item{lat}{Latitude coordinates.}
#'   \item{group}{Group information.}
#'   \item{fitted}{Fitted values.}
#'   \item{resid}{Residuals.}
#'   \item{coeff_hat}{Estimated coefficients.}
#'   \item{mean_resid}{Mean residual by spatial unit.}
#'   \item{sd_resid}{Residual standard deviation by spatial unit.}
#'   \item{pvalue_lb}{Benjamini-Yekutieli adjusted Ljung-Box p-values.}
#'   \item{pvalue_jb}{Benjamini-Yekutieli adjusted Jarque-Bera p-values.}
#'   \item{max_eigen_a}{Maximum eigenvalue diagnostic.}
#'   \item{mu_means}{Mean-equation fixed-effect means.}
#'   \item{mu_tmeans}{Time-weighted mean-equation fixed-effect means.}
#' }
#'
#' @details
#' Residual diagnostics are computed row-wise. Ljung-Box and Jarque-Bera
#' p-values are adjusted using the Benjamini-Yekutieli method.
#'
#' If `obj_stime$error` is not `NULL`, the function reports the estimation error
#' and returns `NULL`.
#'
#' @examples
#' \dontrun{
#' assembled_result <- fit_sdpd_assemble(obj_stime)
#'}
#' @seealso [p.adjust()]
#'
#' @export
fit_sdpd_assemble <- function(obj_stime) {
  if (is.null(obj_stime)) {
    cat("\nA NULL object has been received.")
    result <- NULL
  } else if (is.null(obj_stime$error)) {
    indices <- as.character(obj_stime$px)

    px <- obj_stime$px
    names(px) <- indices

    mean_resid <- apply(obj_stime$resid, 1, mean, na.rm = TRUE)
    sd_resid <- apply(obj_stime$resid, 1, stats::sd, na.rm = TRUE)

    pvalue_lb <- apply(obj_stime$resid, 1, fun_lb_test)
    pvalue_lb <- stats::p.adjust(pvalue_lb, method = "BY")

    pvalue_jb <- apply(obj_stime$resid, 1, fun_jb_test)
    pvalue_jb <- stats::p.adjust(pvalue_jb, method = "BY")

    eigen_a <- round(
      rep(obj_stime$diagnostics, dim(obj_stime$resid)[1]),
      2
    )

    names(mean_resid) <-
      names(sd_resid) <-
      names(pvalue_lb) <-
      names(pvalue_jb) <-
      names(eigen_a) <-
      dimnames(obj_stime$resid)[[1]]

    result <- list(
      px = px,
      lon = obj_stime$lon,
      lat = obj_stime$lat,
      group = data.frame(obj_stime$group)[indices, ],
      fitted = data.frame(obj_stime$fitted, check.names = FALSE)[indices, ],
      resid = data.frame(obj_stime$resid, check.names = FALSE)[indices, ],
      coeff_hat = data.frame(obj_stime$coeff_hat, check.names = FALSE)[indices, ],
      mean_resid = mean_resid[indices],
      sd_resid = sd_resid[indices],
      pvalue_lb = pvalue_lb[indices],
      pvalue_jb = pvalue_jb[indices],
      max_eigen_a = eigen_a[indices],
      mu_means = data.frame(
        obj_stime$mean_equation$mus,
        check.names = FALSE
      )[indices, ],
      mu_tmeans = data.frame(
        obj_stime$mean_equation$mus_t_weighted,
        check.names = FALSE
      )[indices, ]
    )
  } else {
    cat(
      "\nThere was an error in the estimation:",
      obj_stime$group[1, ],
      obj_stime$error
    )

    result <- NULL
  }

  result
}


