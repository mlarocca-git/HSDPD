.test_sdpd_design_from_checked_data <- function(data, model_obj) {
  .new_sdpd_design(
    series = data$series,
    x = data$xx,
    unit_index = rownames(data$series),
    time_index = colnames(data$series),
    pp = data$pp,
    nn = data$nn,
    kk = data$kk,
    ww = data$ww,
    px_neighbors = data$px_neighbors,
    mu = data$mu,
    model = model_obj
  )
}

#' Test an SDP-D Model by Residual Bootstrap
#'
#' Performs bootstrap-based hypothesis tests on fitted SDP-D model coefficients.
#'
#' The function approximates the sampling distribution of estimated SDP-D
#' coefficients using residual bootstrap. It supports several null hypotheses,
#' including zero coefficients, constant coefficients, grouped coefficients, no
#' spatial effects, no autoregressive effects, no covariate effects, and selected
#' coefficient constraints.
#'
#' @param res_fit Fitted SDP-D model result object.
#' @param px Vector. Spatial-unit identifiers used in the fitted model.
#' @param model Optional SDP-D model object. If `NULL`, `res_fit$model` is used.
#' @param n_boot Integer. Number of bootstrap replications. Defaults to `399`.
#' @param h0 Character or integer. Null hypothesis to test. Supported values are
#'   `"zero"`, `"constant"`, `"grouped"`, `"nospatial"`,
#'   `"noautoregressive"`, `"no_x"`, and `"constrained"`.
#' @param coeff_h0 Numeric scalar or matrix. Coefficient value under the null
#'   when `h0 = "zero"`. Defaults to `0`.
#' @param group_index Vector. Group index used when `h0 = "grouped"`. Defaults
#'   to one common group.
#' @param boot_options List. Bootstrap options. Expected elements are
#'   `markovian`, `resid`, `sigma_resid`, `boot_plot`, `folder`, `y_limits`,
#'   and `label_index`.
#'
#' @return A list containing:
#' \describe{
#'   \item{pvalue}{Bootstrap-test p-values.}
#'   \item{n_boot}{Number of bootstrap replications.}
#'   \item{h0}{Null hypothesis tested.}
#'   \item{diagnostics_model}{Model diagnostics.}
#'   \item{coeff_hat}{Estimated coefficients from the fitted model.}
#'   \item{diagnostics_coeff_boot}{Bootstrap coefficient diagnostics.}
#'   \item{diagnostics_sdevs_tsboot}{Bootstrap time-series deviation diagnostics.}
#'   \item{warnings}{Model warnings, if available.}
#' }
#'
#' If validation fails, a list with an `error` or `errors` element is returned.
#'
#' @details
#' The bootstrap can use either fitted residuals, `boot_options$resid =
#' "fitted"`, or newly simulated normal residuals, `boot_options$resid =
#' "normal"`.
#'
#' The test uses a normal basic bootstrap approximation:
#' \deqn{
#'   p = 2 \Phi \left(-\left|\frac{\hat\theta - \theta_0}{s^*}\right|\right)
#' }
#' where `s^*` is the bootstrap standard deviation.
#'
#' @examples
#' \dontrun{
#' test_result <- test_sdpd_model(
#'   res_fit = fit,
#'   px = fit$px,
#'   h0 = "zero",
#'   n_boot = 399
#' )
#'}
#' @seealso
#' [check_sdpd_model()],
#' [fit_sdpd_series()],
#' [fit_sdpd_covs()],
#' [fit_sdpd_coefficients()],
#' [plot_sdpd_boot_series()]
#'
#' @export
test_sdpd_model <- function(res_fit,
                            px,
                            model = NULL,
                            n_boot = 399,
                            h0 = c(
                              "zero",
                              "constant",
                              "grouped",
                              "nospatial",
                              "noautoregressive",
                              "no_x",
                              "constrained"
                            )[1],
                            coeff_h0 = 0,
                            group_index = rep(1, dim(res_fit$fitted)[1]),
                            boot_options = list(
                              markovian = FALSE,
                              resid = c("fitted", "normal")[1],
                              sigma_resid = 1,
                              boot_plot = FALSE,
                              folder = "",
                              y_limits = NULL,
                              label_index = NULL
                            )) {
  # h0 defines the null hypothesis.

  if (is.null(model)) {
    model_obj <- res_fit$model
  } else {
    model_obj <- model
  }

  # Check stationarity conditions.
  diagnostics_model <- check_sdpd_model(res_fit = res_fit)

  if (!is.null(diagnostics_model$errors) &&
      length(diagnostics_model$errors) > 0) {
    cat("\nThere are errors:\n", diagnostics_model$errors)

    return(list(
      errors = diagnostics_model$errors,
      warnings = diagnostics_model$warnings
    ))
  }

  if (!is.null(diagnostics_model$warnings) &&
      length(diagnostics_model$warnings) > 0) {
    cat("\nThere are warnings:\n", diagnostics_model$warnings)
  }

  data <- res_fit$data
  data$pp <- dim(data$series)[1]
  data$nn <- dim(data$series)[2]
  data$kk <- sum(model_obj$beta_coeffs)
  data$mu <- if (!is.null(data$mu)) data$mu else apply(data$series, 1, mean, na.rm = TRUE)

  # Sample distribution approximation based on residual bootstrap.
  if (n_boot <= 20) {
    return(list(error = "The number of bootstrap replications is insufficient."))
  }

  if (boot_options$resid == "fitted") {
    residuals <- res_fit$resid[, -1, drop = FALSE]
  } else if (boot_options$resid == "normal") {
    if (length(boot_options$sigma_resid) == 1) {
      sigma_vector <- rep(boot_options$sigma_resid, data$pp)
    } else if (length(boot_options$sigma_resid) >= data$pp) {
      sigma_vector <- boot_options$sigma_resid[seq_len(data$pp)]
    } else {
      sigma_vector <- rep(1, data$pp)
    }

    residuals <- t(apply(
      as.matrix(sigma_vector),
      1,
      FUN = function(xx, nn) stats::rnorm(n = nn, sd = xx),
      nn = data$nn
    ))
  } else {
    return(list(error = "The value of boot_options$resid must be either 'fitted' or 'normal'."))
  }

  residual_index <- seq_len(dim(residuals)[2])

  boot_index <- matrix(
    sample(
      residual_index,
      size = data$nn * n_boot,
      replace = TRUE
    ),
    nrow = n_boot
  )

  boot_rep <- array(
    0,
    dim = c(n_boot, data$pp, dim(res_fit$coeff_hat)[2])
  )

  dimnames(boot_rep)[[3]] <- dimnames(res_fit$coeff_hat)[[2]]

  devs_tsboot <- array(
    0,
    dim = c(n_boot, 2, data$pp)
  )

  design <- .test_sdpd_design_from_checked_data(data, model_obj)

  # Bootstrap iterations.
  for (bb in seq_len(n_boot)) {
    residual_boot <- residuals[, boot_index[bb, ], drop = FALSE]

    yy_star_1 <- res_fit$data$series
    boot_design <- design

    for (step in seq_len(5)) {
      boot_design$series <- yy_star_1
      yy_star_1 <- .fit_sdpd_series_from_design(
        design = boot_design,
        coeff_hat = res_fit$coeff_hat,
        resids = residual_boot,
        markovian = boot_options$markovian
      )$series
    }

    boot_design$series <- yy_star_1
    covs <- .fit_sdpd_covs_from_design(boot_design)

    boot_rep[bb, , ] <- .fit_sdpd_coefficients_from_design(
      design = boot_design,
      covs = covs
    )$coeff_hat

    devs_tsboot[bb, , ] <- apply(
      yy_star_1 - res_fit$data$series,
      1,
      FUN = function(x) {
        c(mean = mean(x), sd = stats::sd(x))
      }
    )

    if (bb %% 300 == 0 && isTRUE(boot_options$boot_plot)) {
      plot_sdpd_boot_series(
        bb = bb,
        pp = data$pp,
        nn = data$nn,
        eigen = diagnostics_model$diagnostics[1, "max_mod_eigen_a"],
        yy_star_1 = yy_star_1,
        boot_options = boot_options,
        res_fit = res_fit,
        group_index = group_index
      )
    }
  }

  # Collect results.
  results <- coeff_constr <- array(
    0,
    dim = c(data$pp, dim(boot_rep)[3])
  )

  dimnames(results)[[1]] <-
    dimnames(coeff_constr)[[1]] <-
    dimnames(boot_rep)[[2]] <-
    dimnames(res_fit$coeff_hat)[[1]]

  dimnames(results)[[2]] <-
    dimnames(coeff_constr)[[2]] <-
    dimnames(res_fit$coeff_hat)[[2]]

  # Define test statistics for each null hypothesis type.
  if (h0 == "zero" || h0 == 1) {
    boot_constr <- boot_rep
    boot_constr[, , ] <- 0

    coeff_constr[,] <- coeff_h0
  } else if (h0 == "constant" || h0 == 2) {
    boot_constr <- apply(
      boot_rep,
      c(1, 3),
      FUN = function(x) rep(mean(x, na.rm = TRUE), length(x))
    )

    boot_constr <- apply(
      boot_constr,
      c(1, 3),
      FUN = function(x) x
    )

    coeff_constr[,] <- apply(
      res_fit$coeff_hat,
      2,
      FUN = function(x) rep(mean(x, na.rm = TRUE), length(x))
    )
  } else if (h0 == "grouped" || h0 == 3) {
    group_constrain <- function(x, index) {
      x_constrained <- x

      for (jj in seq_len(max(index))) {
        x_constrained[index == jj] <- rep(
          mean(x[index == jj], na.rm = TRUE),
          sum(index == jj)
        )
      }

      x_constrained
    }

    boot_constr <- apply(
      boot_rep,
      c(1, 3),
      FUN = group_constrain,
      index = group_index
    )

    boot_constr <- apply(
      boot_constr,
      c(1, 3),
      FUN = function(x) x
    )

    coeff_constr[,] <- apply(
      res_fit$coeff_hat,
      2,
      FUN = group_constrain,
      index = group_index
    )
  } else if (h0 == "nospatial" || h0 == 4) {
    boot_constr <- boot_rep

    spatial_coefficients <- intersect(
      c("lambda_0", "lambda_2"),
      dimnames(boot_rep)[[3]]
    )

    boot_constr[, , spatial_coefficients] <- 0

    coeff_constr[,] <- res_fit$coeff_hat
    coeff_constr[, spatial_coefficients] <- 0
  } else if (h0 == "noautoregressive" || h0 == 5) {
    boot_constr <- boot_rep

    if ("lambda_1" %in% dimnames(boot_rep)[[3]]) {
      boot_constr[, , "lambda_1"] <- 0
      coeff_constr[,] <- res_fit$coeff_hat
      coeff_constr[, "lambda_1"] <- 0
    } else {
      coeff_constr[,] <- res_fit$coeff_hat
    }
  } else if (h0 == "no_x" || h0 == "noX" || h0 == 6) {
    beta_names <- names(model_obj$beta_coeffs)[model_obj$beta_coeffs]

    boot_constr <- boot_rep

    beta_names <- intersect(beta_names, dimnames(boot_rep)[[3]])

    boot_constr[, , beta_names] <- 0

    coeff_constr[,] <- res_fit$coeff_hat
    coeff_constr[, beta_names] <- 0
  } else if (h0 == "constrained" || h0 == "constrain" || h0 == 7) {
    boot_constr <- boot_rep

    if (all(c("lambda_1", "lambda_2") %in% dimnames(boot_rep)[[3]])) {
      lambda_1_values <- boot_rep[, , "lambda_1"]
      lambda_2_values <- boot_rep[, , "lambda_2"]

      boot_constr[, , "lambda_1"] <- -lambda_2_values
      boot_constr[, , "lambda_2"] <- -lambda_1_values
    }

    coeff_constr[,] <- 0
  } else {
    return(list(error = paste("\nThe null hypothesis", h0, "is not available.")))
  }

  # Summary of bootstrap distribution.
  statistic <- res_fit$coeff_hat - coeff_constr
  boot_rep <- boot_rep - boot_constr

  summarize_boot <- function(x) {
    n_na_x <- sum(is.na(x))
    sd_x <- stats::sd(x, na.rm = TRUE)
    mean_x <- mean(x, na.rm = TRUE)

    if (n_na_x == 0) {
      p_value <- fun_jb_test(x)
    } else {
      p_value <- NA_real_
    }

    c(mean_x, sd_x, p_value, n_na_x)
  }

  boot <- apply(boot_rep, c(2, 3), summarize_boot)

  dimnames(boot)[[1]] <- c(
    "bias_boot",
    "sd_boot",
    "pv_jb_normaltest_boot",
    "nas_boot"
  )

  devs_tsboot <- t(apply(devs_tsboot, c(2, 3), mean))

  dimnames(devs_tsboot)[[2]] <- c("mean", "sd")

  # Test using normal basic bootstrap.
  pvalue <- 2 * stats::pnorm(
    abs(statistic / boot["sd_boot", , ]),
    lower.tail = FALSE
  )

  # Results.
  dimnames(pvalue)[[1]] <- dimnames(devs_tsboot)[[1]] <- dimnames(boot)[[2]]
  dimnames(pvalue)[[2]] <- dimnames(boot)[[3]]

  boot["bias_boot", , ] <- boot["bias_boot", , ] - res_fit$coeff_hat

  list(
    pvalue = pvalue,
    n_boot = n_boot,
    h0 = h0,
    diagnostics_model = diagnostics_model$diagnostics,
    coeff_hat = res_fit$coeff_hat,
    diagnostics_coeff_boot = boot,
    diagnostics_sdevs_tsboot = devs_tsboot,
    warnings = diagnostics_model$warnings
  )
}
