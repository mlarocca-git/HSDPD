#' Generate an SDP-D Model
#'
#' Generates SDP-D model components for simulation.
#'
#' The function generates a spatial weight structure, proximity-neighbor
#' structure, coefficient matrix, innovation standard deviations, and optional
#' index weights. It is mainly used internally by [build_sdpd_model()] when
#' `sim = TRUE`.
#'
#' @param pp Integer. Number of spatial units.
#' @param model SDP-D model object or partial SDP-D model specification,
#'   typically created by [build_sdpd_model()].
#' @param ww_index Optional matrix. Spatial-neighbor index matrix. If `NULL`,
#'   it is generated according to `sim_options$ww_type`.
#' @param ww_values Optional matrix. Spatial-weight value matrix. If `NULL`,
#'   it is generated and row-normalized.
#' @param px_neighbors Optional list. Pixel or proximity-neighbor structure. If
#'   `NULL`, it is generated from `ww_index`.
#' @param sim_options List. Simulation options controlling spatial weights,
#'   neighbors, coefficient ranges, constant-by-group parameters, and innovation
#'   standard deviations.
#'
#' @return A list containing:
#' \describe{
#'   \item{coeffs}{Generated coefficient matrix.}
#'   \item{ww_values}{Generated or supplied spatial-weight values.}
#'   \item{ww_index}{Generated or supplied spatial-neighbor index matrix.}
#'   \item{px_neighbors}{Generated or supplied proximity-neighbor structure.}
#'   \item{sigma_eps}{Innovation standard deviations.}
#'   \item{index_weights}{Index weights used for group-constant parameters.}
#' }
#'
#' @details
#' If `ww_index` or `ww_values` is missing, the function generates a spatial
#' structure according to `sim_options$ww_type`. Supported values are:
#' \itemize{
#'   \item `1` or `"rook"`;
#'   \item `2` or `"queen"`;
#'   \item `3` or `"corr"`.
#' }
#'
#' The generated spatial matrix is row-normalized using the L1 norm.
#'
#' Coefficients are generated from uniform distributions using the ranges
#' specified in `sim_options`. If a coefficient is marked as constant, one value
#' is generated within each `index_weights` group.
#'
#' @examples
#' \dontrun{
#' generated_model <- generate_sdpd_model(
#'   pp = 100,
#'   model = model,
#'   ww_index = NULL,
#'   ww_values = NULL,
#'   px_neighbors = NULL,
#'   sim_options = sim_options
#' )
#'}
#' @seealso [build_sdpd_model()]
#'
#' @export
generate_sdpd_model <- function(pp,
                                model,
                                ww_index,
                                ww_values,
                                px_neighbors,
                                sim_options) {
  # Generate the spatial weight matrix.
  if (is.null(ww_index) || is.null(ww_values)) {
    has_valid_ww <- FALSE

    while (!has_valid_ww) {
      if (sim_options$ww_type == 1 || sim_options$ww_type == "rook") {
        ww <- ww_index <- matrix(0, ncol = 5, nrow = pp)

        ww_index[, 1] <- seq_len(pp)
        ww_index[, 2] <- ifelse(seq_len(pp) + 1 > pp, 0, seq_len(pp) + 1)
        ww_index[, 3] <- ifelse(seq_len(pp) + 2 > pp, 0, seq_len(pp) + 2)
        ww_index[, 4] <- ifelse(seq_len(pp) - 1 < 1, 0, seq_len(pp) - 1)
        ww_index[, 5] <- ifelse(seq_len(pp) - 2 < 1, 0, seq_len(pp) - 2)

        ww[, 1] <- 0
        ww[, 2:5] <- 1

        has_valid_ww <- TRUE
      } else if (sim_options$ww_type == 2 || sim_options$ww_type == "queen") {
        ww <- ww_index <- matrix(0, ncol = 9, nrow = pp)

        ww_index[, 1] <- seq_len(pp)
        ww_index[, 2] <- ifelse(seq_len(pp) + 1 > pp, 0, seq_len(pp) + 1)
        ww_index[, 3] <- ifelse(seq_len(pp) + 2 > pp, 0, seq_len(pp) + 2)
        ww_index[, 4] <- ifelse(seq_len(pp) - 1 < 1, 0, seq_len(pp) - 1)
        ww_index[, 5] <- ifelse(seq_len(pp) - 2 < 1, 0, seq_len(pp) - 2)
        ww_index[, 6] <- ifelse(seq_len(pp) + 3 > pp, 0, seq_len(pp) + 3)
        ww_index[, 7] <- ifelse(seq_len(pp) + 4 > pp, 0, seq_len(pp) + 4)
        ww_index[, 8] <- ifelse(seq_len(pp) - 3 < 1, 0, seq_len(pp) - 3)
        ww_index[, 9] <- ifelse(seq_len(pp) - 4 < 1, 0, seq_len(pp) - 4)

        ww[, 1] <- 0
        ww[, 2:9] <- 1

        has_valid_ww <- TRUE
      } else if (sim_options$ww_type == 3 || sim_options$ww_type == "corr") {
        ww <- matrix(stats::rnorm(pp * pp, sd = 1), ncol = pp)
        ww <- ww %*% t(ww)
        diag(ww) <- 0

        has_valid_ww <- det(ww) != 0

        ww_index <- matrix(
          rep(seq_len(pp), pp),
          nrow = pp,
          byrow = TRUE
        )
      }
    }

    dimnames(ww_index)[[1]] <- seq_len(pp)
    dimnames(ww)[[1]] <- seq_len(pp)
  } else {
    ww <- ww_values
  }

  # Row-normalize the spatial matrix by the L1 norm.
  ww_values <- t(apply(
    ww,
    1,
    FUN = function(x) {
      x / sum(abs(x))
    }
  ))

  # Generate matrix with neighbor points.
  if (is.null(px_neighbors)) {
    n_neighbor_cols <- (2 * sim_options$px_neighbors + 1)^2 - 1

    neighbors <- matrix(
      NA,
      nrow = pp,
      ncol = n_neighbor_cols + 1
    )

    neighbors[, seq_len(min(n_neighbor_cols, dim(ww_index)[2]))] <-
      ww_index[, seq_len(min(n_neighbor_cols, dim(ww_index)[2]))]

    identify_neighbors <- function(x, ww_i, n_col) {
      result <- rep(NA, n_col)

      neighbors_old <- x[abs(x) > 0][-1]
      neighbors_old <- stats::na.exclude(neighbors_old)

      candidate_neighbors <- as.vector(ww_i[as.character(neighbors_old), ])
      neighbors_new <- candidate_neighbors[abs(candidate_neighbors) > 0]
      neighbors_new <- unique(neighbors_new)
      neighbors_new <- neighbors_new[
        !(neighbors_new %in% c(x[1], neighbors_old))
      ]

      neighbors_new <- c(neighbors_old, neighbors_new)

      if (length(neighbors_new) >= n_col) {
        result[seq_len(n_col)] <- neighbors_new[seq_len(n_col)]
      } else if (length(neighbors_new) > 0) {
        result[seq_along(neighbors_new)] <- neighbors_new
      }

      c(x[1], result)
    }

    for (zz in seq_len(100)) {
      old_neighbors <- neighbors

      neighbors <- t(apply(
        old_neighbors,
        1,
        FUN = identify_neighbors,
        ww_i = ww_index,
        n_col = n_neighbor_cols
      ))
    }

    neighbors <- neighbors[, -1]
    dimnames(neighbors)[[1]] <- dimnames(ww_index)[[1]]

    px_neighbors <- list(
      index = neighbors,
      series_boundary = NULL
    )
  }

  # Generate lambda and beta parameters.
  lambda_names <- names(model$lambda_coeffs)[model$lambda_coeffs]
  beta_names <- names(model$beta_coeffs)[model$beta_coeffs]
  fixed_effects_name <- "fixed_effects"[model$fixed_effects]

  parameters <- matrix(
    0,
    nrow = pp,
    ncol = sum(model$lambda_coeffs) +
      sum(model$beta_coeffs) +
      sum(model$fixed_effects)
  )

  dimnames(parameters)[[2]] <- c(
    lambda_names,
    beta_names,
    fixed_effects_name
  )

  # Fixed-effect values.
  if (model$fixed_effects) {
    if (sim_options$fixed_effects$constant[1]) {
      parameters[, "fixed_effects"] <- rep(
        stats::runif(
          1,
          sim_options$fixed_effects$min[1],
          sim_options$fixed_effects$max[1]
        ),
        pp
      )
    } else {
      parameters[, "fixed_effects"] <- stats::runif(
        pp,
        sim_options$fixed_effects$min[1],
        sim_options$fixed_effects$max[1]
      )
    }
  }

  # Set default index weights.
  if (is.null(sim_options$index_weights)) {
    sim_options$index_weights <- rep(1, pp)
  }

  # Beta coefficients for exogenous regressors.
  kk <- length(beta_names)

  if (kk > 0) {
    if (is.null(sim_options$betas$constant)) {
      sim_options$betas$constant <- rep(FALSE, kk)
    }

    if (length(sim_options$betas$constant) == 1 ||
        length(sim_options$betas$constant) != kk) {
      sim_options$betas$constant <- rep(sim_options$betas$constant[1], kk)
    }

    if (length(sim_options$betas$min) == 1 ||
        length(sim_options$betas$min) != kk) {
      sim_options$betas$min <- rep(sim_options$betas$min[1], kk)
    }

    if (length(sim_options$betas$max) == 1 ||
        length(sim_options$betas$max) != kk) {
      sim_options$betas$max <- rep(sim_options$betas$max[1], kk)
    }

    for (ss in seq_len(kk)) {
      if (sim_options$betas$constant[ss]) {
        for (jj in seq_len(max(sim_options$index_weights))) {
          parameters[
            sim_options$index_weights == jj,
            beta_names[ss]
          ] <- rep(
            round(
              stats::runif(
                1,
                min = sim_options$betas$min[ss],
                max = sim_options$betas$max[ss]
              ),
              digits = 4
            ),
            sum(sim_options$index_weights == jj)
          )
        }
      } else {
        parameters[, beta_names[ss]] <- round(
          stats::runif(
            pp,
            min = sim_options$betas$min[ss],
            max = sim_options$betas$max[ss]
          ),
          digits = 4
        )
      }
    }
  }

  # Generate lambda coefficients.
  if (model$lambda_coeffs["lambda_0"] &&
      is.null(sim_options$lambda_0$constant)) {
    sim_options$lambda_0$constant <- FALSE
  }

  if (model$lambda_coeffs["lambda_1"] &&
      is.null(sim_options$lambda_1$constant)) {
    sim_options$lambda_1$constant <- FALSE
  }

  if (model$lambda_coeffs["lambda_2"] &&
      is.null(sim_options$lambda_2$constant)) {
    sim_options$lambda_2$constant <- FALSE
  }

  if (model$lambda_coeffs["lambda_0"]) {
    if (sim_options$lambda_0$constant) {
      for (jj in seq_len(max(sim_options$index_weights))) {
        parameters[
          sim_options$index_weights == jj,
          "lambda_0"
        ] <- rep(
          round(
            stats::runif(
              1,
              min = sim_options$lambda_0$min,
              max = sim_options$lambda_0$max
            ),
            digits = 4
          ),
          sum(sim_options$index_weights == jj)
        )
      }
    } else {
      parameters[, "lambda_0"] <- round(
        stats::runif(
          pp,
          min = sim_options$lambda_0$min,
          max = sim_options$lambda_0$max
        ),
        digits = 4
      )
    }
  }

  if (model$lambda_coeffs["lambda_1"]) {
    if (sim_options$lambda_1$constant) {
      for (jj in seq_len(max(sim_options$index_weights))) {
        parameters[
          sim_options$index_weights == jj,
          "lambda_1"
        ] <- rep(
          round(
            stats::runif(
              1,
              min = sim_options$lambda_1$min,
              max = sim_options$lambda_1$max
            ),
            digits = 4
          ),
          sum(sim_options$index_weights == jj)
        )
      }
    } else {
      parameters[, "lambda_1"] <- round(
        stats::runif(
          pp,
          min = sim_options$lambda_1$min,
          max = sim_options$lambda_1$max
        ),
        digits = 4
      )
    }
  }

  if (model$lambda_coeffs["lambda_2"]) {
    if (sim_options$lambda_2$constant) {
      for (jj in seq_len(max(sim_options$index_weights))) {
        parameters[
          sim_options$index_weights == jj,
          "lambda_2"
        ] <- rep(
          round(
            stats::runif(
              1,
              min = sim_options$lambda_2$min,
              max = sim_options$lambda_2$max
            ),
            digits = 4
          ),
          sum(sim_options$index_weights == jj)
        )
      }
    } else {
      parameters[, "lambda_2"] <- round(
        stats::runif(
          pp,
          min = sim_options$lambda_2$min,
          max = sim_options$lambda_2$max
        ),
        digits = 4
      )
    }
  }

  # Generate sigma_eps values.
  sigma_eps <- numeric(pp)

  if (is.null(sim_options$sigma_eps$constant)) {
    sim_options$sigma_eps$constant <- TRUE
  }

  if (length(sim_options$sigma_eps$constant) != 1) {
    sim_options$sigma_eps$constant <- sim_options$sigma_eps$constant[1]
  }

  if (sim_options$sigma_eps$constant) {
    for (jj in seq_len(max(sim_options$index_weights))) {
      sigma_eps[sim_options$index_weights == jj] <- rep(
        round(
          stats::runif(
            1,
            min = sim_options$sigma_eps$min[1],
            max = sim_options$sigma_eps$max[1]
          ),
          digits = 4
        ),
        sum(sim_options$index_weights == jj)
      )
    }
  } else {
    sigma_eps <- round(
      stats::runif(
        pp,
        min = sim_options$sigma_eps$min[1],
        max = sim_options$sigma_eps$max[1]
      ),
      digits = 4
    )
  }

  # Return model parameters.
  dimnames(parameters)[[1]] <- seq_len(pp)
  names(sigma_eps) <- seq_len(pp)
  names(sim_options$index_weights) <- seq_len(pp)

  list(
    coeffs = parameters,
    ww_values = ww_values,
    ww_index = ww_index,
    px_neighbors = px_neighbors,
    sigma_eps = sigma_eps,
    index_weights = sim_options$index_weights
  )
}
