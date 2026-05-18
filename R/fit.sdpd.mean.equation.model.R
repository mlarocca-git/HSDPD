#' Compute the SDP-D Mean-Equation Model
#'
#' Computes mean-equation summaries from fitted SDP-D values.
#'
#' The function computes spatial-unit means over the full period, spatially
#' lagged means using the spatial weight matrix, and optionally weighted spatial
#' means. If time weights are supplied, it also computes summaries by time group
#' or by numeric time weights.
#'
#' @param result List. SDP-D estimation result object containing at least a
#'   `fitted` matrix.
#' @param ww Numeric matrix. Spatial weight matrix.
#' @param time_weights Optional factor or numeric vector. If a factor, means are
#'   computed separately for each time group. If numeric, weighted means over
#'   time are computed.
#' @param index_weights Optional numeric vector. Spatial-unit weights used to
#'   compute weighted neighborhood means.
#'
#' @return A list containing:
#' \describe{
#'   \item{mus}{Data frame with full-period mean-equation summaries.}
#'   \item{mus_t_weighted}{Optional data frame or list with time-grouped or
#'   time-weighted mean-equation summaries.}
#' }
#'
#' @details
#' The `mus` component always contains:
#' \itemize{
#'   \item `mu_i`: mean fitted value for each spatial unit;
#'   \item `mu_around_i`: spatially lagged mean fitted value.
#' }
#'
#' If `index_weights` is supplied, `mus` also contains
#' `mu_j_weighted_around_i`, the weighted neighborhood mean.
#'
#' If `time_weights` is a factor, the function computes separate means for each
#' factor level. If `time_weights` is numeric, the function computes
#' time-weighted means.
#'
#' @examples
#' \dontrun{
#' mean_equation <- fit_sdpd_mean_equation_model(
#'   result = result,
#'   ww = ww
#' )
#'
#' mean_equation <- fit_sdpd_mean_equation_model(
#'   result = result,
#'   ww = ww,
#'   time_weights = time_weights,
#'   index_weights = index_weights
#' )
#'}
#'
#' @export
fit_sdpd_mean_equation_model <- function(result,
                                         ww,
                                         time_weights = NULL,
                                         index_weights = NULL) {

  # Compute means over the full period.
  mu_i <- apply(result$fitted, 1, mean, na.rm = TRUE)
  mu_around_i <- as.vector(ww %*% mu_i)

  mus <- cbind(
    mu_i = mu_i,
    mu_around_i = mu_around_i
  )

  if (!is.null(index_weights)) {
    weights_around_i <- as.vector(ww %*% index_weights)

    mu_j_weighted_around_i <- as.vector(
      ww %*% (mu_i * index_weights)
    ) / weights_around_i

    mus <- cbind(
      mus,
      mu_j_weighted_around_i = mu_j_weighted_around_i
    )
  }

  # Compute means over time groups, for example seasons, or by time weighting.
  mu_group_i <- NULL
  mu_group_around_i <- NULL
  mu_j_weighted_group_around_i <- NULL
  mus_t_weighted <- NULL

  if (!is.null(time_weights)) {
    if (is.factor(time_weights)) {
      groups <- levels(time_weights)

      for (group in groups) {
        mu_temp <- apply(
          result$fitted[, time_weights == group],
          1,
          mean,
          na.rm = TRUE
        )

        mu_group_i <- cbind(mu_group_i, mu_temp)
        mu_group_around_i <- cbind(
          mu_group_around_i,
          as.vector(ww %*% mu_temp)
        )

        if (!is.null(index_weights)) {
          mu_j_weighted_group_around_i <- cbind(
            mu_j_weighted_group_around_i,
            as.vector(ww %*% (mu_temp * index_weights)) / weights_around_i
          )
        }
      }

      dimnames(mu_group_i)[[2]] <- groups
      dimnames(mu_group_around_i)[[2]] <- groups

      if (!is.null(index_weights)) {
        dimnames(mu_j_weighted_group_around_i)[[2]] <- groups
      }

      mus_t_weighted <- data.frame(
        mu_i = data.frame(mu_group_i),
        mu_around_i = data.frame(mu_group_around_i),
        mu_j_weighted_around_i = data.frame(mu_j_weighted_group_around_i)
      )
    } else if (is.numeric(time_weights)) {
      total_time_weight <- sum(time_weights)

      mus_t_weighted <- list()

      mus_t_weighted$mu_i <- apply(
        result$fitted,
        1,
        FUN = function(x, w) {
          sum(x * w, na.rm = TRUE) / total_time_weight
        },
        w = time_weights
      )

      mus_t_weighted$mu_around_i <- as.vector(
        ww %*% mus_t_weighted$mu_i
      )

      if (!is.null(index_weights)) {
        mus_t_weighted$mu_j_weighted_around_i <- as.vector(
          ww %*% (mus_t_weighted$mu_i * index_weights)
        ) / weights_around_i
      }
    }
  }

  # Return results.
  list(
    mus = data.frame(mus),
    mus_t_weighted = mus_t_weighted
  )
}



