#' Build an SDP-D Model Object
#'
#' Builds an `sdpd_model` object, which is the main input used by subsequent
#' estimation functions in the package.
#'
#' The function defines the active spatial-dynamic components, endogenous and
#' exogenous variables, fixed effects, time effects, spatial weight structures,
#' and optional simulation settings. If `sim = TRUE`, the model components are
#' randomly generated using the parameters supplied through `sim_options`.
#'
#' If `check = TRUE`, and the generated or supplied model contains coefficient
#' and spatial-weight information, the function also runs model diagnostics using
#' `check_sdpd_model()`.
#'
#' @param endogenous Character scalar. Name of the endogenous variable.
#'   Defaults to `"var_y"`.
#' @param lambda_0 Logical scalar. Whether the first lambda component is active.
#'   Defaults to `TRUE`.
#' @param lambda_1 Logical scalar. Whether the second lambda component is active.
#'   Defaults to `TRUE`.
#' @param lambda_2 Logical scalar. Whether the third lambda component is active.
#'   Defaults to `TRUE`.
#' @param covariates Character or numeric vector. Names or indices of exogenous
#'   variables to include in the model. Use `0` for no covariates.
#'   Defaults to `0`.
#' @param fixed_effects Logical scalar. Whether fixed effects are included.
#'   Defaults to `TRUE`.
#' @param time_effects Logical scalar. Whether time effects are included.
#'   Defaults to `FALSE`.
#' @param coeffs Optional numeric matrix or data frame. Coefficient matrix used
#'   to infer available variable names and model components. If supplied, it must
#'   have column names.
#' @param check Logical scalar. Whether model diagnostics should be computed
#'   when enough information is available. Defaults to `TRUE`.
#' @param ww_index Optional object defining the index structure of the spatial
#'   weight matrix.
#' @param ww_values Optional object defining the values of the spatial weight
#'   matrix.
#' @param px_neighbors Optional object defining pixel or proximity neighbors.
#' @param rr_groups Optional matrix or data frame defining regional groups. If
#'   `NULL` and `sim = TRUE`, groups are generated automatically.
#' @param sim Logical scalar. Whether the model should be simulated. Defaults to
#'   `FALSE`.
#' @param pp Optional integer. Cross-sectional dimension of the model. Required
#'   when `sim = TRUE`.
#' @param d_i Integer. Approximate group-size divisor used when groups are
#'   generated automatically. Defaults to `100`.
#' @param sim_options List. Simulation options used when `sim = TRUE`.
#'   Expected elements include `ww_type`, `px_neighbors`, `index_group`,
#'   `fixed_effects`, `lambda_0`, `lambda_1`, `lambda_2`, `betas`, and
#'   `sigma_eps`.
#'
#' @return A list of classless model components representing an SDP-D model.
#'   The returned object may contain the following elements:
#'   \describe{
#'     \item{lambda_coeffs}{Logical vector identifying active lambda components.}
#'     \item{name_endogenous}{Name of the endogenous variable.}
#'     \item{beta_coeffs}{Logical vector identifying active beta coefficients.}
#'     \item{fixed_effects}{Logical value indicating whether fixed effects are included.}
#'     \item{time_effects}{Logical value indicating whether time effects are included.}
#'     \item{coeffs}{Coefficient matrix, when generated or supplied.}
#'     \item{ww_index}{Spatial-weight index structure, when available.}
#'     \item{ww_values}{Spatial-weight values, when available.}
#'     \item{px_neighbors}{Pixel or proximity-neighbor information, when available.}
#'     \item{sigma_eps}{Error-term variance parameter, when available.}
#'     \item{groups}{Regional group structure, when available.}
#'     \item{diagnostics}{Model diagnostics, when `check = TRUE` and diagnostics can be computed.}
#'   }
#'
#' @details
#' At least one of `lambda_0`, `lambda_1`, or `lambda_2` must be active.
#'
#' If `coeffs` is supplied, the function uses its column names to verify whether
#' requested covariates and fixed effects are available. Covariates not present
#' in `coeffs` are removed from the model with a warning.
#'
#' When `sim = TRUE`, `pp` must be supplied. If `rr_groups` is not supplied, the
#' function automatically builds a group structure using `pp` and `d_i`.
#'
#' @examples
#' model <- build_sdpd_model(
#'   endogenous = "var_y",
#'   covariates = c("x_1", "x_2"),
#'   fixed_effects = TRUE,
#'   time_effects = FALSE,
#'   sim = FALSE
#' )
#'
#' simulated_model <- build_sdpd_model(
#'   endogenous = "var_y",
#'   covariates = c("x_1", "x_2"),
#'   sim = TRUE,
#'   pp = 100
#' )
#'
#' @seealso [generate_sdpd_model()], [check_sdpd_model()]
#'
#' @export
build_sdpd_model <- function(endogenous = "var_y",
                             lambda_0 = TRUE,
                             lambda_1 = TRUE,
                             lambda_2 = TRUE,
                             covariates = 0,
                             fixed_effects = TRUE,
                             time_effects = FALSE,
                             coeffs = NULL,
                             check = TRUE,
                             ww_index = NULL,
                             ww_values = NULL,
                             px_neighbors = NULL,
                             rr_groups = NULL,
                             sim = FALSE,
                             pp = NULL,
                             d_i = 100,
                             sim_options = list(
                               ww_type = "queen",
                               px_neighbors = 8,
                               index_group = rep(1, pp),
                               fixed_effects = list(constant = FALSE, min = 10, max = 20),
                               lambda_0 = list(constant = FALSE, min = -0.5, max = 0.5),
                               lambda_1 = list(constant = FALSE, min = -0.5, max = 0.5),
                               lambda_2 = list(constant = FALSE, min = -0.5, max = 0.5),
                               betas = list(
                                 constant = FALSE,
                                 min = c(0, 0, 0, 1),
                                 max = c(1, 1, 1, 4)
                               ),
                               sigma_eps = list(constant = TRUE, min = 1, max = 1)
                             )) {

  ## This function builds the sdpd_model object, which is the main input
  ## for the subsequent estimation functions.
  ##
  ## If sim = TRUE, all model components are randomly generated using
  ## the parameters passed through sim_options.
  ##
  ## If check = TRUE, the structure and stationarity of the reduced-form
  ## model are also checked afterwards.

  sdpd_model <- list()

  # The following conditions check that the inputs can be interpreted
  # as Boolean values.
  if (lambda_0) {
    lambda_0 <- TRUE
  }

  if (lambda_1) {
    lambda_1 <- TRUE
  }

  if (lambda_2) {
    lambda_2 <- TRUE
  }

  # Define lambda components.
  if (lambda_0 + lambda_1 + lambda_2 == 0) {
    return(list(error = "At least one lambda component must be active"))
  }

  sdpd_model$lambda_coeffs <- c(lambda_0, lambda_1, lambda_2)

  lambda_names <- c("lambda_0", "lambda_1", "lambda_2")
  names(sdpd_model$lambda_coeffs) <- lambda_names
  lambda_names <- lambda_names[sdpd_model$lambda_coeffs]

  # Initialize variable names possibly provided through coeffs.
  if (is.null(coeffs)) {
    coeff_names <- NULL
  } else if (is.data.frame(coeffs) || is.matrix(coeffs)) {
    coeffs <- as.matrix(coeffs)
    coeff_names <- dimnames(coeffs)[[2]]

    if (is.null(coeff_names)) {
      return(list(error = "The coefficient matrix must have defined column names"))
    }
  } else {
    return(list(error = "The coeffs parameter has an invalid format: it must be a numeric matrix or dataframe"))
  }

  # Check endogenous variable name/index.
  if (is.vector(endogenous) && is.character(endogenous)) {
    if (length(endogenous) > 1) {
      return(list(error = "The endogenous variable must be unique; multiple values are not allowed"))
    }

    if (sum(endogenous %in% covariates) > 0) {
      return(list(error = "The endogenous variable cannot be included among the covariates"))
    }

    sdpd_model$name_endogenous <- endogenous
  } else {
    return(list(error = "The endogenous variable has an invalid format: it must be a character value"))
  }

  # Check exogenous variable names/indices.
  if (is.vector(covariates) && (is.numeric(covariates) || is.character(covariates))) {
    if (covariates[1] == 0) {
      sdpd_model$beta_coeffs <- NULL
      beta_names <- NULL
    } else if (sum(covariates %in% coeff_names) > 0) {
      beta_names <- covariates[covariates %in% coeff_names]

      if (length(beta_names) < length(covariates)) {
        cat("\nWarning: some covariates were removed from the model because they are not present in the coeffs matrix\n")
      }

      sdpd_model$beta_coeffs <- rep(TRUE, length(beta_names))
      names(sdpd_model$beta_coeffs) <- beta_names
    } else if (is.null(coeff_names)) {
      sdpd_model$beta_coeffs <- rep(TRUE, length(covariates))
      beta_names <- names(sdpd_model$beta_coeffs) <- covariates
    }
  } else {
    return(list(error = "The exogenous variables parameter has an invalid format: it must be a numeric or character vector"))
  }

  # Define the fixed-effects component.
  if (!is.null(coeff_names) && fixed_effects) {
    if ("fixed_effects" %in% coeff_names) {
      sdpd_model$fixed_effects <- TRUE
      fixed_effects_name <- "fixed_effects"
    } else {
      sdpd_model$fixed_effects <- FALSE
      cat("\nWarning: fixed effects were removed from the model because they are not present in the coeffs matrix.\n")
      fixed_effects_name <- NULL
    }
  } else {
    sdpd_model$fixed_effects <- fixed_effects
    fixed_effects_name <- "fixed_effects"[fixed_effects]
  }

  # Define the time-effects component.
  if (is.logical(time_effects) && length(time_effects) == 1) {
    sdpd_model$time_effects <- time_effects
  } else {
    return(list(error = "The time_effects parameter has an invalid format: it must be a Boolean value"))
  }

  # Simulated model case.
  if (sim) {
    if (is.null(pp)) {
      return(list(error = "When coefficients are simulated, the pp dimension must be provided"))
    }

    if (is.null(rr_groups)) {
      n_groups <- pp %/% d_i

      groups <- rep(n_groups, pp)
      groups[seq_along(rep(seq(1, n_groups), each = pp %/% n_groups))] <-
        rep(seq(1, n_groups), each = pp %/% n_groups)

      group_labels <- paste("group", groups, sep = "_")
      groups <- cbind(cod = groups, label = group_labels)
      dimnames(groups)[[1]] <- seq(1, pp)
    } else {
      groups <- rr_groups
    }

    new_model <- generate_sdpd_model(
      pp = pp,
      model = sdpd_model,
      ww_index = ww_index,
      ww_values = ww_values,
      px_neighbors = px_neighbors,
      sim_options = sim_options
    )

    sdpd_model$coeffs <- new_model$coeffs
    sdpd_model$ww_index <- new_model$ww_index
    sdpd_model$ww_values <- new_model$ww_values
    sdpd_model$px_neighbors <- new_model$px_neighbors
    sdpd_model$sigma_eps <- new_model$sigma_eps
    sdpd_model$groups <- groups
  }

  # Diagnostics.
  if (!is.null(sdpd_model$coeffs) &&
      !is.null(sdpd_model$ww_index) &&
      !is.null(sdpd_model$ww_values) &&
      check) {
    sdpd_model$diagnostics <- check_sdpd_model(model = sdpd_model)$diagnostics[1, ]
  }

  # Output.
  sdpd_model
}
