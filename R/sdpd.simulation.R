.as_simulation_matrix <- function(x, row_names = NULL) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.vector(x) && !is.list(x)) {
    x_names <- names(x)
    x <- matrix(x, nrow = 1)

    if (!is.null(x_names)) {
      colnames(x) <- x_names
    }
  } else {
    x <- as.matrix(x)
  }

  storage.mode(x) <- "numeric"

  has_default_rownames <- identical(rownames(x), as.character(seq_len(nrow(x))))

  if ((is.null(rownames(x)) || has_default_rownames) &&
      !is.null(row_names) &&
      nrow(x) == length(row_names)) {
    rownames(x) <- row_names
  }

  x
}

.fit_result_column <- function(fit_results, names) {
  for (name in names) {
    if (!is.null(fit_results[[name]])) {
      return(fit_results[[name]])
    }
  }

  NULL
}

.fit_result_coefficients <- function(fit_results, row_names) {
  coefficient_estimates <- .fit_result_column(fit_results, c("coeff_hat", "coeff.hat"))

  if (!is.null(coefficient_estimates)) {
    return(.as_simulation_matrix(coefficient_estimates, row_names))
  }

  coefficient_columns <- grep(
    "^(coeff_hat|coeff\\.hat)\\.",
    names(fit_results),
    value = TRUE
  )

  if (length(coefficient_columns) == 0) {
    return(NULL)
  }

  coefficient_estimates <- as.matrix(fit_results[, coefficient_columns, drop = FALSE])
  colnames(coefficient_estimates) <- sub(
    "^(coeff_hat|coeff\\.hat)\\.",
    "",
    coefficient_columns
  )

  .as_simulation_matrix(coefficient_estimates, row_names)
}

.align_simulation_matrix <- function(source, template) {
  aligned <- matrix(
    NA_real_,
    nrow = nrow(template),
    ncol = ncol(template),
    dimnames = dimnames(template)
  )

  if (is.null(source)) {
    return(aligned)
  }

  source_rows <- intersect(rownames(template), rownames(source))
  source_cols <- intersect(colnames(template), colnames(source))

  if (length(source_rows) > 0 && length(source_cols) > 0) {
    aligned[source_rows, source_cols] <- source[source_rows, source_cols, drop = FALSE]
  }

  aligned
}

.residual_statistic_vector <- function(fit_results, name, row_names) {
  current_names <- switch(
    name,
    mean = c("mean_resid", "mean.resid", "mean"),
    sd = c("sd_resid", "sd.resid", "sd"),
    pvalue.LB = c("pvalue_lb", "pvalue.LB", "pvalue_lb_adjusted"),
    pvalue.JB = c("pvalue_jb", "pvalue.JB", "pvalue_jb_adjusted"),
    max.eigenA = c("max_eigen_a", "max.eigenA", "max_mod_eigen_a"),
    max_eigen_a = c("max_eigen_a", "max.eigenA", "max_mod_eigen_a"),
    c(name)
  )

  values <- .fit_result_column(fit_results, current_names)

  if (is.null(values)) {
    return(rep(NA_real_, length(row_names)))
  }

  values <- as.numeric(values)
  names(values) <- names(.fit_result_column(fit_results, current_names))

  if (is.null(names(values))) {
    names(values) <- row_names[seq_along(values)]
  }

  values[row_names]
}

#' Increment SDP-D Simulation Partial Sums
#'
#' Updates cumulative simulation summaries from one fitted SDP-D result.
#'
#' @param partial_sum List. Partial simulation accumulator containing
#'   `estimation.bias`, `estimation.mse`, `resid.sim`, `contatore1`, and
#'   `contatore2`.
#' @param fit_results Data frame or list. Fitted SDP-D results from
#'   [fit_sdpd_model()] or a compatible object containing coefficient estimates
#'   and residual diagnostics.
#' @param simulated_coefficients Matrix or data frame. True coefficients used in
#'   the simulation.
#'
#' @return The updated `partial_sum` list.
#'
#' @details
#' The function preserves the legacy accumulator field names used by older
#' simulation scripts. Bias contributions are computed as estimated coefficient
#' minus simulated coefficient, and MSE contributions are the squared
#' differences. `contatore1` and `contatore2` are incremented only for
#' non-missing contributions.
#'
#' Residual diagnostics are read from current names such as `mean_resid`,
#' `sd_resid`, `pvalue_lb`, `pvalue_jb`, and `max_eigen_a`, while also accepting
#' legacy names such as `pvalue.LB`, `pvalue.JB`, and `max.eigenA`.
#'
#' @examples
#' partial_sum <- list(
#'   estimation.bias = data.frame(lambda_1 = c("1" = 0)),
#'   estimation.mse = data.frame(lambda_1 = c("1" = 0)),
#'   resid.sim = data.frame(mean = c("1" = 0)),
#'   contatore1 = data.frame(lambda_1 = c("1" = 0)),
#'   contatore2 = data.frame(mean = c("1" = 0))
#' )
#' fit_results <- data.frame(
#'   coeff_hat = data.frame(lambda_1 = c("1" = 0.3)),
#'   mean_resid = c("1" = 0.2)
#' )
#' simulated_coefficients <- data.frame(lambda_1 = c("1" = 0.1))
#'
#' increment_partial_sum(partial_sum, fit_results, simulated_coefficients)
#'
#' @export
increment_partial_sum <- function(partial_sum,
                                  fit_results,
                                  simulated_coefficients) {
  required_fields <- c(
    "estimation.bias",
    "estimation.mse",
    "resid.sim",
    "contatore1",
    "contatore2"
  )

  missing_fields <- required_fields[!required_fields %in% names(partial_sum)]

  if (length(missing_fields) > 0) {
    stop(
      "partial_sum is missing required fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  bias_sum <- .as_simulation_matrix(partial_sum$estimation.bias)
  mse_sum <- .as_simulation_matrix(partial_sum$estimation.mse, rownames(bias_sum))
  coefficient_counter <- .as_simulation_matrix(partial_sum$contatore1, rownames(bias_sum))
  residual_sum <- .as_simulation_matrix(partial_sum$resid.sim)
  residual_counter <- .as_simulation_matrix(partial_sum$contatore2, rownames(residual_sum))

  coefficient_estimates <- .fit_result_coefficients(fit_results, rownames(bias_sum))
  true_coefficients <- .as_simulation_matrix(
    simulated_coefficients,
    rownames(bias_sum)
  )

  coefficient_estimates <- .align_simulation_matrix(coefficient_estimates, bias_sum)
  true_coefficients <- .align_simulation_matrix(true_coefficients, bias_sum)

  coefficient_difference <- coefficient_estimates - true_coefficients
  valid_coefficients <- !is.na(coefficient_difference)

  bias_sum[valid_coefficients] <- bias_sum[valid_coefficients] +
    coefficient_difference[valid_coefficients]
  mse_sum[valid_coefficients] <- mse_sum[valid_coefficients] +
    coefficient_difference[valid_coefficients]^2
  coefficient_counter[valid_coefficients] <- coefficient_counter[valid_coefficients] + 1

  for (statistic_name in colnames(residual_sum)) {
    statistic_values <- .residual_statistic_vector(
      fit_results = fit_results,
      name = statistic_name,
      row_names = rownames(residual_sum)
    )
    valid_values <- !is.na(statistic_values)

    residual_sum[valid_values, statistic_name] <-
      residual_sum[valid_values, statistic_name] + statistic_values[valid_values]
    residual_counter[valid_values, statistic_name] <-
      residual_counter[valid_values, statistic_name] + 1
  }

  partial_sum$estimation.bias <- data.frame(bias_sum, check.names = FALSE)
  partial_sum$estimation.mse <- data.frame(mse_sum, check.names = FALSE)
  partial_sum$resid.sim <- data.frame(residual_sum, check.names = FALSE)
  partial_sum$contatore1 <- data.frame(coefficient_counter, check.names = FALSE)
  partial_sum$contatore2 <- data.frame(residual_counter, check.names = FALSE)

  partial_sum
}

#' Summarise SDP-D Simulation Results
#'
#' Builds compact column summaries from averaged SDP-D simulation outputs.
#'
#' @param results List or data frame. Simulation results containing one or more
#'   of `estimation.bias`, `estimation.mse`, and `resid.sim`.
#' @param na_rm Logical scalar. Whether missing values should be removed before
#'   computing column means. Defaults to `TRUE`.
#'
#' @return A data frame with `component`, `statistic`, and `value` columns.
#'
#' @export
summarise_sdpd_simulation <- function(results, na_rm = TRUE) {
  components <- c("estimation.bias", "estimation.mse", "resid.sim")

  summaries <- lapply(
    components[components %in% names(results)],
    function(component) {
      values <- .as_simulation_matrix(results[[component]])
      data.frame(
        component = component,
        statistic = colnames(values),
        value = colMeans(values, na.rm = na_rm),
        row.names = NULL,
        check.names = FALSE
      )
    }
  )

  if (length(summaries) == 0) {
    return(data.frame(
      component = character(),
      statistic = character(),
      value = numeric()
    ))
  }

  do.call(rbind, summaries)
}
