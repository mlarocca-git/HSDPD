build.sdpd.model <- function(...) {
  .Deprecated("build_sdpd_model")
  build_sdpd_model(...)
}

#' Build an SDP-D Series Object
#'
#' Deprecated wrapper for [build_sdpd_series()].
#'
#' @param df.obj Optional data object used by parallelized workflows.
#' @param px Optional vector of spatial-unit identifiers.
#' @param rry Optional endogenous data series.
#' @param rrxx Optional exogenous regressors.
#' @param lon Optional longitude vector or column identifier.
#' @param lat Optional latitude vector or column identifier.
#' @param rrgroups Optional group matrix or data frame.
#' @param label_groups Optional labels for raster group codes.
#' @param model SDP-D model object.
#' @param check Logical scalar. Whether to check the generated series.
#' @param ww.index Optional spatial-neighbor index matrix.
#' @param ww.values Optional spatial-weight value matrix.
#' @param px.neighbors Optional pixel-neighbor structure.
#' @param SIM Logical scalar. Whether to simulate the data.
#' @param nn Optional integer. Number of time observations for simulation.
#' @param vec.options List of vectorization and simulation options.
#'
#' @return See [build_sdpd_series()].
#' @seealso [build_sdpd_series()]
#' @export
build.sdpd.series <- function(df.obj = NULL,
                              px = NULL,
                              rry = NULL,
                              rrxx = NULL,
                              lon = NULL,
                              lat = NULL,
                              rrgroups = NULL,
                              label_groups = NULL,
                              model,
                              check = TRUE,
                              ww.index = NULL,
                              ww.values = NULL,
                              px.neighbors = NULL,
                              SIM = FALSE,
                              nn = NULL,
                              vec.options = list(
                                px.core = 1,
                                px.neighbors = 6,
                                na.rm = TRUE,
                                NAcovs = "pairwise.complete.obs",
                                covariates.sim.model = list(ar = c(0.8), sd = 1),
                                markovian = TRUE,
                                num.steps = 10
                              )) {
  .Deprecated("build_sdpd_series")

  if (!is.null(df.obj) && is.null(px)) {
    px <- df.obj$px
  }

  # Map old vec.options names to the new snake_case API.
  if (!is.null(vec.options$px.core) && is.null(vec.options$px_core)) {
    vec.options$px_core <- vec.options$px.core
  }

  if (!is.null(vec.options$px.neighbors) && is.null(vec.options$px_neighbors)) {
    vec.options$px_neighbors <- vec.options$px.neighbors
  }

  if (!is.null(vec.options$na.rm) && is.null(vec.options$na_rm)) {
    vec.options$na_rm <- vec.options$na.rm
  }

  if (!is.null(vec.options$NAcovs) && is.null(vec.options$na_covs)) {
    vec.options$na_covs <- vec.options$NAcovs
  }

  if (!is.null(vec.options$covariates.sim.model) &&
      is.null(vec.options$covariates_sim_model)) {
    vec.options$covariates_sim_model <- vec.options$covariates.sim.model
  }

  if (!is.null(vec.options$num.steps) && is.null(vec.options$num_steps)) {
    vec.options$num_steps <- vec.options$num.steps
  }

  build_sdpd_series(
    df_obj = df.obj,
    px = px,
    rr_y = rry,
    rr_xx = rrxx,
    lon = lon,
    lat = lat,
    rr_groups = rrgroups,
    label_groups = label_groups,
    model = model,
    check = check,
    ww_index = ww.index,
    ww_values = ww.values,
    px_neighbors = px.neighbors,
    sim = SIM,
    nn = nn,
    vec_options = vec.options
  )
}


#' Build a Spatial Weight Matrix
#'
#' Deprecated wrapper for [build_spatial_matrix()].
#'
#' @param ww.index Matrix. Spatial-neighbor index matrix.
#' @param ww.values Matrix. Spatial-weight value matrix.
#'
#' @return See [build_spatial_matrix()].
#' @seealso [build_spatial_matrix()]
#' @export
build.spatial.matrix <- function(ww.index, ww.values) {
  .Deprecated("build_spatial_matrix")

  build_spatial_matrix(
    ww_index = ww.index,
    ww_values = ww.values
  )
}

#' Check an SDP-D Model
#'
#' Deprecated wrapper for [check_sdpd_model()].
#'
#' @param res.fit Optional fitted SDP-D model object.
#' @param model Optional SDP-D model object.
#' @param ww.index Optional spatial-neighbor index matrix.
#' @param ww.values Optional spatial-weight value matrix.
#'
#' @return See [check_sdpd_model()].
#' @seealso [check_sdpd_model()]
#' @export
check.sdpd.model <- function(res.fit = NULL,
                             model = NULL,
                             ww.index = NULL,
                             ww.values = NULL) {
  .Deprecated("check_sdpd_model")

  check_sdpd_model(
    res_fit = res.fit,
    model = model,
    ww_index = ww.index,
    ww_values = ww.values
  )
}


#' Check an SDP-D Series Object
#'
#' Deprecated wrapper for [check_sdpd_series()].
#'
#' @param series Matrix, data frame, or list. Endogenous series or SDP-D series
#'   object.
#' @param XX Optional matrix, data frame, or three-dimensional array of
#'   regressors.
#' @param model Optional SDP-D model object.
#' @param ww.index Optional spatial-neighbor index matrix.
#' @param ww.values Optional spatial-weight value matrix.
#' @param px.neighbors Optional pixel or proximity-neighbor object.
#' @param px Optional spatial-unit or pixel structure.
#' @param lat Optional latitude vector.
#' @param lon Optional longitude vector.
#' @param group Optional group object.
#' @param index.weights Optional spatial-unit weights.
#' @param time.weights Optional time weights.
#'
#' @return See [check_sdpd_series()].
#' @seealso [check_sdpd_series()]
#' @export
check.sdpd.series <- function(series,
                              XX = NULL,
                              model = NULL,
                              ww.index = NULL,
                              ww.values = NULL,
                              px.neighbors = NULL,
                              px = NULL,
                              lat = NULL,
                              lon = NULL,
                              group = NULL,
                              index.weights = NULL,
                              time.weights = NULL) {
  .Deprecated("check_sdpd_series")

  check_sdpd_series(
    series = series,
    xx = XX,
    model = model,
    ww_index = ww.index,
    ww_values = ww.values,
    px_neighbors = px.neighbors,
    px = px,
    lat = lat,
    lon = lon,
    group = group,
    index_weights = index.weights,
    time_weights = time.weights
  )
}

#' Fit the Second-Stage SDP-D Estimator
#'
#' Deprecated wrapper for [fit_second_stage()].
#'
#' @param dseries Numeric matrix. Endogenous data series.
#' @param X.centr Numeric matrix or array of centered regressors.
#' @param W Numeric spatial weight matrix.
#' @param model SDP-D model object.
#' @param coeff.hat Numeric matrix of coefficient estimates.
#' @param data Optional data summary object.
#'
#' @return See [fit_second_stage()].
#' @seealso [fit_second_stage()]
#' @export
fit.2nd.stage <- function(dseries,
                          X.centr,
                          W,
                          model,
                          coeff.hat,
                          data = NULL) {
  .Deprecated("fit_second_stage")

  fit_second_stage(
    data_series = dseries,
    x_centered = X.centr,
    ww = W,
    model = model,
    coeff_hat = coeff.hat,
    data = data
  )
}

#' Assemble SDP-D Estimation Results
#'
#' Deprecated wrapper for [fit_sdpd_assemble()].
#'
#' @param obj.stime List. Single spatial-time estimation object.
#'
#' @return See [fit_sdpd_assemble()].
#' @seealso [fit_sdpd_assemble()]
#' @export
fit.sdpd.assemble <- function(obj.stime) {
  .Deprecated("fit_sdpd_assemble")

  fit_sdpd_assemble(obj_stime = obj.stime)
}

#' Estimate SDP-D Coefficients
#'
#' Deprecated wrapper for [fit_sdpd_coefficients()].
#'
#' @param W Numeric matrix. Spatial weight matrix.
#' @param COVs List. Covariance components.
#' @param mu Numeric vector. Spatial-unit means.
#' @param model SDP-D model object.
#'
#' @return See [fit_sdpd_coefficients()].
#' @seealso [fit_sdpd_coefficients()]
#' @export
fit.sdpd.coefficients <- function(W, COVs, mu, model) {
  .Deprecated("fit_sdpd_coefficients")

  fit_sdpd_coefficients(
    ww = W,
    covs = COVs,
    mu = mu,
    model = model
  )
}

#' Compute SDP-D Covariance Components
#'
#' Deprecated wrapper for [fit_sdpd_covs()].
#'
#' @param series Numeric matrix. Endogenous data series.
#' @param X Optional numeric matrix or three-dimensional array of regressors.
#' @param kk Integer. Number of covariates.
#' @param px.neighbors Optional pixel or proximity-neighbor information.
#' @param nn Integer. Number of time observations.
#' @param pp Integer. Number of spatial units.
#' @param NAcovs Character scalar. Missing-value handling method passed to
#'   [stats::cov()].
#'
#' @return See [fit_sdpd_covs()].
#' @seealso [fit_sdpd_covs()]
#' @export
fit.sdpd.covs <- function(series,
                          X = NULL,
                          kk = 0,
                          px.neighbors = NULL,
                          nn,
                          pp,
                          NAcovs = "pairwise.complete.obs") {
  .Deprecated("fit_sdpd_covs")

  fit_sdpd_covs(
    series = series,
    x = X,
    kk = kk,
    px_neighbors = px.neighbors,
    nn = nn,
    pp = pp,
    na_covs = NAcovs
  )
}

#' Compute the SDP-D Mean-Equation Model
#'
#' Deprecated wrapper for [fit_sdpd_mean_equation_model()].
#'
#' @param res List. SDP-D estimation result object.
#' @param WW Numeric matrix. Spatial weight matrix.
#' @param time.weights Optional factor or numeric vector of time weights.
#' @param index.weights Optional numeric vector of spatial-unit weights.
#'
#' @return See [fit_sdpd_mean_equation_model()].
#' @seealso [fit_sdpd_mean_equation_model()]
#' @export
fit.sdpd.mean.equation.model <- function(res,
                                         WW,
                                         time.weights = NULL,
                                         index.weights = NULL) {
  .Deprecated("fit_sdpd_mean_equation_model")

  fit_sdpd_mean_equation_model(
    result = res,
    ww = WW,
    time_weights = time.weights,
    index_weights = index.weights
  )
}

#' Fit an SDP-D Model
#'
#' Deprecated wrapper for [fit_sdpd_model()].
#'
#' @param series SDP-D series object.
#' @param model SDP-D model object.
#' @param check Logical scalar. Whether model checks should be performed.
#' @param two.stage Logical scalar. Whether the second-stage estimator should be
#'   used.
#' @param NAcovs Character scalar. Missing-value handling method passed to
#'   covariance computations.
#' @param parallelize Logical scalar. Whether the series should be split into
#'   parallelized spatial groups before estimation.
#'
#' @return See [fit_sdpd_model()].
#' @seealso [fit_sdpd_model()]
#' @export
fit.sdpd.model <- function(series,
                           model,
                           check = FALSE,
                           two.stage = FALSE,
                           NAcovs = "pairwise.complete.obs",
                           parallelize = FALSE) {
  .Deprecated("fit_sdpd_model")

  fit_sdpd_model(
    series = series,
    model = model,
    check = check,
    two_stage = two.stage,
    na_covs = NAcovs,
    parallelize = parallelize
  )
}

#' Run the SDP-D Estimation Procedure
#'
#' Deprecated wrapper for [fit_sdpd_procedure()].
#'
#' @param series SDP-D series object.
#' @param model SDP-D model object.
#' @param check Logical scalar. Whether stationarity diagnostics should be
#'   computed.
#' @param two.stage Logical scalar. Whether the second-stage estimator should be
#'   applied.
#' @param NAcovs Character scalar. Missing-value handling method passed to
#'   covariance computations.
#'
#' @return See [fit_sdpd_procedure()].
#' @seealso [fit_sdpd_procedure()]
#' @export
fit.sdpd.procedure <- function(series,
                               model,
                               check = FALSE,
                               two.stage = FALSE,
                               NAcovs = "pairwise.complete.obs") {
  .Deprecated("fit_sdpd_procedure")

  fit_sdpd_procedure(
    series = series,
    model = model,
    check = check,
    two_stage = two.stage,
    na_covs = NAcovs
  )
}



#' Compute Fitted Values and Residuals for an SDP-D Series
#'
#' Deprecated wrapper for [fit_sdpd_series()].
#'
#' @param dseries Numeric matrix. Endogenous data series.
#' @param W Numeric matrix. Spatial weight matrix.
#' @param X.centr Optional numeric matrix or three-dimensional array of centered
#'   regressors.
#' @param model SDP-D model object.
#' @param coeff.hat Numeric matrix of estimated coefficients.
#' @param time_effects Optional numeric vector of time effects.
#' @param px.sim Optional vector identifying spatial units to update during
#'   simulation.
#' @param resids Optional numeric matrix of residuals.
#' @param markovian Logical scalar. Whether to use the Markovian update.
#' @param num.steps Integer. Number of iterative update steps.
#'
#' @return See [fit_sdpd_series()].
#' @seealso [fit_sdpd_series()]
#' @export
fit.sdpd.series <- function(dseries,
                            W,
                            X.centr = NULL,
                            model,
                            coeff.hat,
                            time_effects = NULL,
                            px.sim = NULL,
                            resids = NULL,
                            markovian = TRUE,
                            num.steps = 10) {
  .Deprecated("fit_sdpd_series")

  fit_sdpd_series(
    data_series = dseries,
    ww = W,
    x_centered = X.centr,
    model = model,
    coeff_hat = coeff.hat,
    time_effects = time_effects,
    px_sim = px.sim,
    resids = resids,
    markovian = markovian,
    num_steps = num.steps
  )
}


#' Compute a Jarque-Bera Normality-Test P-Value
#'
#' Deprecated wrapper for [fun_jb_test()].
#'
#' @param x Numeric vector. Input values to test for normality.
#' @param ... Additional arguments passed to [fun_jb_test()].
#'
#' @return See [fun_jb_test()].
#' @seealso [fun_jb_test()]
#' @export
fun.JBtest <- function(x, ...) {
  .Deprecated("fun_jb_test")

  fun_jb_test(x = x, ...)
}


#' Compute a Ljung-Box Test P-Value
#'
#' Deprecated wrapper for [fun_lb_test()].
#'
#' @param x Numeric vector. Input values to test for serial autocorrelation.
#' @param ... Additional arguments passed to [fun_lb_test()].
#'
#' @return See [fun_lb_test()].
#' @seealso [fun_lb_test()]
#' @export
fun.LBtest <- function(x, ...) {
  .Deprecated("fun_lb_test")

  fun_lb_test(x = x, ...)
}


#' Plot Spatial Residual Statistics
#'
#' Deprecated wrapper for [plot_stat_discrete_resids()].
#'
#' @param df.results List or data frame. SDP-D estimation results.
#' @param statistic Function. Summary function applied to residuals.
#' @param main Optional character scalar. Plot title.
#' @param significant.test Logical scalar. Whether values should be converted
#'   into significance classes.
#' @param BYadjusted Logical scalar. Whether p-values should be adjusted using
#'   the Benjamini-Yekutieli method.
#' @param alpha Optional numeric scalar. Significance level.
#' @param mid_value Numeric scalar. Midpoint used for the diverging color scale.
#' @param size.point Numeric scalar. Point size used in the plot.
#' @param ... Additional arguments passed to `statistic`.
#'
#' @return See [plot_stat_discrete_resids()].
#' @seealso [plot_stat_discrete_resids()]
#' @export
fun.plot.stat.discrete.RESIDs <- function(df.results,
                                          statistic = mean,
                                          main = NULL,
                                          significant.test = FALSE,
                                          BYadjusted = FALSE,
                                          alpha = NULL,
                                          mid_value = 0,
                                          size.point = 1,
                                          ...) {
  .Deprecated("plot_stat_discrete_resids")

  plot_stat_discrete_resids(
    results = df.results,
    statistic = statistic,
    main = main,
    significant_test = significant.test,
    by_adjusted = BYadjusted,
    alpha = alpha,
    mid_value = mid_value,
    size_point = size.point,
    ...
  )
}

#' Generate an SDP-D Series
#'
#' Deprecated wrapper for [generate_sdpd_series()].
#'
#' @param nn Optional integer. Number of time observations to simulate.
#' @param rry Optional numeric matrix. Initial endogenous series.
#' @param rrxx Optional numeric matrix, data frame, or three-dimensional array
#'   of exogenous regressors.
#' @param model SDP-D model object.
#' @param markovian Logical scalar. Whether to use the Markovian update.
#' @param num.steps Integer. Number of iterative update steps.
#' @param covariates.sim.model List. ARIMA model specification passed to
#'   [stats::arima.sim()].
#'
#' @return See [generate_sdpd_series()].
#' @seealso [generate_sdpd_series()]
#' @export
generate.sdpd.series <- function(nn,
                                 rry,
                                 rrxx,
                                 model,
                                 markovian,
                                 num.steps,
                                 covariates.sim.model) {
  .Deprecated("generate_sdpd_series")

  generate_sdpd_series(
    nn = nn,
    rr_y = rry,
    rr_xx = rrxx,
    model = model,
    markovian = markovian,
    num_steps = num.steps,
    covariates_sim_model = covariates.sim.model
  )
}



#' Generate an SDP-D Model
#'
#' Deprecated wrapper for [generate_sdpd_model()].
#'
#' @param pp Integer. Number of spatial units.
#' @param model SDP-D model object.
#' @param ww.index Optional spatial-neighbor index matrix.
#' @param ww.values Optional spatial-weight value matrix.
#' @param px.neighbors Optional proximity-neighbor structure.
#' @param sim.options List. Simulation options.
#'
#' @return See [generate_sdpd_model()].
#' @seealso [generate_sdpd_model()]
#' @export
generate.sdpd.model <- function(pp,
                                model,
                                ww.index,
                                ww.values,
                                px.neighbors,
                                sim.options) {
  .Deprecated("generate_sdpd_model")

  # Map old option names to the new snake_case API.
  if (!is.null(sim.options$ww.type)) {
    sim.options$ww_type <- sim.options$ww.type
  }

  if (!is.null(sim.options$px.neighbors)) {
    sim.options$px_neighbors <- sim.options$px.neighbors
  }

  if (!is.null(sim.options$fixed.effects)) {
    sim.options$fixed_effects <- sim.options$fixed.effects
  }

  if (!is.null(sim.options$lambda0)) {
    sim.options$lambda_0 <- sim.options$lambda0
  }

  if (!is.null(sim.options$lambda1)) {
    sim.options$lambda_1 <- sim.options$lambda1
  }

  if (!is.null(sim.options$lambda2)) {
    sim.options$lambda_2 <- sim.options$lambda2
  }

  if (!is.null(sim.options$sigma.eps)) {
    sim.options$sigma_eps <- sim.options$sigma.eps
  }

  if (!is.null(sim.options$index.weights)) {
    sim.options$index_weights <- sim.options$index.weights
  }

  generate_sdpd_model(
    pp = pp,
    model = model,
    ww_index = ww.index,
    ww_values = ww.values,
    px_neighbors = px.neighbors,
    sim_options = sim.options
  )
}


#' Parallelize SDP-D Series Data by Group
#'
#' Deprecated wrapper for [parallelize_sdpd_data()].
#'
#' @param obj.series SDP-D series object.
#' @param model SDP-D model object.
#'
#' @return See [parallelize_sdpd_data()].
#' @seealso [parallelize_sdpd_data()]
#' @export
parallelize.sdpd.data <- function(obj.series, model) {
  .Deprecated("parallelize_sdpd_data")

  # Map old series fields to the new snake_case API when needed.
  if (!is.null(obj.series$X) && is.null(obj.series$xx)) {
    obj.series$xx <- obj.series$X
  }

  if (!is.null(obj.series$ww.index) && is.null(obj.series$ww_index)) {
    obj.series$ww_index <- obj.series$ww.index
  }

  if (!is.null(obj.series$ww.values) && is.null(obj.series$ww_values)) {
    obj.series$ww_values <- obj.series$ww.values
  }

  if (!is.null(obj.series$px.neighbors) && is.null(obj.series$px_neighbors)) {
    obj.series$px_neighbors <- obj.series$px.neighbors
  }

  if (!is.null(obj.series$index.weights) && is.null(obj.series$index_weights)) {
    obj.series$index_weights <- obj.series$index.weights
  }

  if (!is.null(obj.series$time.weights) && is.null(obj.series$time_weights)) {
    obj.series$time_weights <- obj.series$time.weights
  }

  parallelize_sdpd_data(
    series_object = obj.series,
    model = model
  )
}

#' Plot a Bootstrap SDP-D Series
#'
#' Deprecated wrapper for [plot_sdpd_boot_series()].
#'
#' @param bb Integer. Current bootstrap iteration.
#' @param pp Integer. Number of spatial units.
#' @param nn Integer. Number of time observations.
#' @param eigen Numeric vector. Eigenvalue diagnostics.
#' @param yy.star1 Numeric matrix. Bootstrapped time series.
#' @param opts.boot List. Bootstrap options.
#' @param res.fit List. Fitted SDP-D model result.
#' @param group.index Vector. Group index.
#'
#' @return Invisibly returns the output file path of the generated JPEG plot.
#' @seealso [plot_sdpd_boot_series()]
#' @export
plot.sdpd.boot.series <- function(bb,
                                  pp,
                                  nn,
                                  eigen,
                                  yy.star1,
                                  opts.boot,
                                  res.fit,
                                  group.index) {
  .Deprecated("plot_sdpd_boot_series")

  # Map old bootstrap-option names to the new snake_case API.
  if (!is.null(opts.boot$ylimiti) && is.null(opts.boot$y_limits)) {
    opts.boot$y_limits <- opts.boot$ylimiti
  }

  if (!is.null(opts.boot$label.index) && is.null(opts.boot$label_index)) {
    opts.boot$label_index <- opts.boot$label.index
  }

  if (!is.null(opts.boot$cartel) && is.null(opts.boot$folder)) {
    opts.boot$folder <- opts.boot$cartel
  }

  # Map old fitted-result field names to the new snake_case API.
  if (!is.null(res.fit$coeff.hat) && is.null(res.fit$coeff_hat)) {
    res.fit$coeff_hat <- res.fit$coeff.hat
  }

  if (!is.null(res.fit$data$W) && is.null(res.fit$data$ww)) {
    res.fit$data$ww <- res.fit$data$W
  }

  if (!is.null(res.fit$data$px.neighbors) &&
      is.null(res.fit$data$px_neighbors)) {
    res.fit$data$px_neighbors <- res.fit$data$px.neighbors
  }

  if (!is.null(res.fit$data$X) && is.null(res.fit$data$xx)) {
    res.fit$data$xx <- res.fit$data$X
  }

  plot_sdpd_boot_series(
    bb = bb,
    pp = pp,
    nn = nn,
    eigen = eigen,
    yy_star_1 = yy.star1,
    boot_options = opts.boot,
    res_fit = res.fit,
    group_index = group.index
  )
}


#' Plot Spatial SDP-D Estimates
#'
#' Deprecated wrapper for [plot_sdpd_estimates()].
#'
#' @param obj.results List or data frame. SDP-D estimation results.
#' @param item Character scalar. Name of the result component to plot.
#' @param main Optional character scalar. Text prepended to each plot title.
#' @param sub Optional character scalar. Text appended to each plot title.
#' @param limits Optional numeric vector of length 2.
#' @param mid Numeric scalar or vector. Midpoint value for the diverging color
#'   scale.
#' @param size.point Numeric scalar. Point size used in the plot.
#'
#' @return See [plot_sdpd_estimates()].
#' @seealso [plot_sdpd_estimates()]
#' @export
plot.sdpd.estimates <- function(obj.results,
                                item = "coeff.hat",
                                main = NULL,
                                sub = NULL,
                                limits = NULL,
                                mid = 0,
                                size.point = 1) {
  .Deprecated("plot_sdpd_estimates")

  # Map old result-field names to the new snake_case API when needed.
  if (item == "coeff.hat") {
    item <- "coeff_hat"

    if (!is.null(obj.results$coeff.hat) &&
        is.null(obj.results$coeff_hat)) {
      obj.results$coeff_hat <- obj.results$coeff.hat
    }
  }

  plot_sdpd_estimates(
    results = obj.results,
    item = item,
    main = main,
    sub = sub,
    limits = limits,
    mid = mid,
    size_point = size.point
  )
}

#' Plot an SDP-D Model Fit
#'
#' Deprecated wrapper for [plot_sdpd_model()].
#'
#' @param res.fit Fitted SDP-D model result object.
#' @param n.units Integer or `"all"`. Number of spatial units to plot.
#' @param n.vars Integer or `"all"`. Number of variables to plot.
#' @param which Integer scalar. Plot type selector.
#' @param t.axis List. Time-axis options with elements `t.labels` and
#'   `t.points`.
#' @param xlimit Optional numeric vector of length 2. X-axis limits.
#' @param ylimit Optional numeric vector of length 2. Y-axis limits.
#' @param max.col Integer. Maximum number of plot columns.
#' @param col.punti Color specification for plotted points or lines.
#'
#' @return See [plot_sdpd_model()].
#' @seealso [plot_sdpd_model()]
#' @export
plot.sdpd.model <- function(res.fit,
                            n.units = "all",
                            n.vars = "all",
                            which = c(1, 2)[1],
                            t.axis = list(t.labels = NULL, t.points = NULL),
                            xlimit = NULL,
                            ylimit = NULL,
                            max.col = 5,
                            col.punti = 1) {
  .Deprecated("plot_sdpd_model")

  t_axis <- list(
    t_labels = t.axis$t.labels,
    t_points = t.axis$t.points
  )

  plot_sdpd_model(
    res_fit = res.fit,
    n_units = n.units,
    n_vars = n.vars,
    which = which,
    t_axis = t_axis,
    x_limit = xlimit,
    y_limit = ylimit,
    max_col = max.col,
    point_col = col.punti
  )
}


#' Plot Fitted, Observed, and Residual SDP-D Series
#'
#' Deprecated wrapper for [plot_sdpd_series_fits()].
#'
#' @param df.results List or data frame. SDP-D estimation results.
#' @param latitude Numeric scalar. Latitude of the location to plot.
#' @param longitude Numeric scalar. Longitude of the location to plot.
#' @param n.digits Integer. Number of digits used to round coordinates.
#' @param main Optional character scalar. Plot title.
#' @param sub Optional character scalar. Plot subtitle.
#' @param xlab Optional character scalar. X-axis label.
#' @param ylab Optional character scalar. Y-axis label.
#'
#' @return See [plot_sdpd_series_fits()].
#' @seealso [plot_sdpd_series_fits()]
#' @export
plot.sdpd.series.FITs <- function(df.results,
                                  latitude,
                                  longitude,
                                  n.digits = 1,
                                  main = NULL,
                                  sub = NULL,
                                  xlab = NULL,
                                  ylab = NULL) {
  .Deprecated("plot_sdpd_series_fits")

  plot_sdpd_series_fits(
    results = df.results,
    latitude = latitude,
    longitude = longitude,
    n_digits = n.digits,
    main = main,
    sub = sub,
    xlab = xlab,
    ylab = ylab
  )
}


#' Plot an SDP-D Series at a Spatial Location
#'
#' Deprecated wrapper for [plot_sdpd_series()].
#'
#' @param rry Data frame. Series data containing `latitude`, `longitude`, and
#'   time-indexed series columns.
#' @param lat Numeric scalar. Latitude of the location to plot.
#' @param lon Numeric scalar. Longitude of the location to plot.
#' @param n.digits Integer. Number of digits used to round coordinates.
#' @param main Optional character scalar. Plot title.
#' @param sub Optional character scalar. Plot subtitle.
#' @param xlab Optional character scalar. X-axis label.
#' @param ylab Optional character scalar. Y-axis label.
#'
#' @return See [plot_sdpd_series()].
#' @seealso [plot_sdpd_series()]
#' @export
plot.sdpd.series <- function(rry,
                             lat,
                             lon,
                             n.digits = 15,
                             main = NULL,
                             sub = NULL,
                             xlab = NULL,
                             ylab = NULL) {
  .Deprecated("plot_sdpd_series")

  plot_sdpd_series(
    rr_y = rry,
    lat = lat,
    lon = lon,
    n_digits = n.digits,
    main = main,
    sub = sub,
    xlab = xlab,
    ylab = ylab
  )
}

#' Plot Detailed SDP-D Model Diagnostics
#'
#' Deprecated wrapper for [plot_1_sdpd_model()].
#'
#' @param res.fit Fitted SDP-D model result object.
#' @param n.units Integer, character vector, logical vector, or `"all"`.
#' @param n.vars Integer, character vector, logical vector, or `"all"`.
#' @param t.axis List. Time-axis options with elements `t.labels` and
#'   `t.points`.
#' @param xlimit Optional numeric vector of length 2. X-axis limits.
#' @param ylimit Optional numeric vector of length 2. Y-axis limits.
#' @param max.col Integer. Maximum number of component columns shown per unit.
#' @param col.punti Color specification for scatter-plot points.
#'
#' @return See [plot_1_sdpd_model()].
#' @seealso [plot_1_sdpd_model()]
#' @export
plot1.sdpd.model <- function(res.fit,
                             n.units = "all",
                             n.vars = "all",
                             t.axis = list(t.labels = NULL, t.points = NULL),
                             xlimit = NULL,
                             ylimit = NULL,
                             max.col = 5,
                             col.punti = 1) {
  .Deprecated("plot_1_sdpd_model")

  # Map old fitted-result field names to the new snake_case API.
  if (!is.null(res.fit$coeff.hat) && is.null(res.fit$coeff_hat)) {
    res.fit$coeff_hat <- res.fit$coeff.hat
  }

  if (!is.null(res.fit$data$W) && is.null(res.fit$data$ww)) {
    res.fit$data$ww <- res.fit$data$W
  }

  if (!is.null(res.fit$data$X) && is.null(res.fit$data$xx)) {
    res.fit$data$xx <- res.fit$data$X
  }

  t_axis <- list(
    t_labels = t.axis$t.labels,
    t_points = t.axis$t.points
  )

  plot_1_sdpd_model(
    res_fit = res.fit,
    n_units = n.units,
    n_vars = n.vars,
    t_axis = t_axis,
    x_limit = xlimit,
    y_limit = ylimit,
    max_col = max.col,
    point_col = col.punti
  )
}


#' Read SDP-D Data from a Data Frame or Matrix
#'
#' Deprecated wrapper for [read_data_from_dataframe()].
#'
#' @param px Optional vector. Spatial-unit identifiers to include.
#' @param latit Optional numeric vector or column index identifying latitudes.
#' @param longit Optional numeric vector or column index identifying longitudes.
#' @param rry Data frame or matrix. Endogenous series.
#' @param rrXX Optional list or array of exogenous regressors.
#' @param rrgroups Optional data frame or matrix with columns `COD` and `LABEL`.
#' @param model SDP-D model object.
#' @param ww.index Matrix. Spatial-neighbor index matrix.
#' @param ww.values Matrix. Spatial-weight value matrix.
#' @param px.neighbors Pixel or proximity-neighbor structure.
#'
#' @return See [read_data_from_dataframe()].
#' @seealso [read_data_from_dataframe()]
#' @export
read.data.from.dataframe <- function(px,
                                     latit = NULL,
                                     longit = NULL,
                                     rry,
                                     rrXX = NULL,
                                     rrgroups = NULL,
                                     model,
                                     ww.index,
                                     ww.values,
                                     px.neighbors) {
  .Deprecated("read_data_from_dataframe")

  if (!is.null(px.neighbors$seriesBoundary) &&
      is.null(px.neighbors$series_boundary)) {
    px.neighbors$series_boundary <- px.neighbors$seriesBoundary
  }

  read_data_from_dataframe(
    px = px,
    lat = latit,
    lon = longit,
    rr_y = rry,
    rr_xx = rrXX,
    rr_groups = rrgroups,
    model = model,
    ww_index = ww.index,
    ww_values = ww.values,
    px_neighbors = px.neighbors
  )
}

#' Read SDP-D Data from Raster Objects
#'
#' Deprecated wrapper for [read_data_from_raster()].
#'
#' @param px Optional vector. Raster cell identifiers to include.
#' @param latit Optional latitude vector.
#' @param longit Optional longitude vector.
#' @param rry `SpatRaster`. Endogenous raster time series.
#' @param rrXX Optional `SpatRaster` or list of `SpatRaster` objects.
#' @param rrgroups Optional group raster.
#' @param label_groups Optional named vector of group labels.
#' @param model SDP-D model object.
#' @param vec.options List. Vectorization options.
#' @param type.w Character scalar. Spatial-weight type.
#'
#' @return See [read_data_from_raster()].
#' @seealso [read_data_from_raster()]
#' @export
read.data.from.raster <- function(px = NULL,
                                  latit = NULL,
                                  longit = NULL,
                                  rry,
                                  rrXX = NULL,
                                  rrgroups = NULL,
                                  label_groups = NULL,
                                  model,
                                  vec.options,
                                  type.w = "distance") {
  .Deprecated("read_data_from_raster")

  if (!is.null(vec.options$px.core) && is.null(vec.options$px_core)) {
    vec.options$px_core <- vec.options$px.core
  }

  if (!is.null(vec.options$px.neighbors) && is.null(vec.options$px_neighbors)) {
    vec.options$px_neighbors <- vec.options$px.neighbors
  }

  if (!is.null(vec.options$na.rm) && is.null(vec.options$na_rm)) {
    vec.options$na_rm <- vec.options$na.rm
  }

  read_data_from_raster(
    px = px,
    lat = latit,
    lon = longit,
    rr_y = rry,
    rr_xx = rrXX,
    rr_groups = rrgroups,
    label_groups = label_groups,
    model = model,
    vec_options = vec.options,
    type_w = type.w
  )
}

#' Test an SDP-D Model by Residual Bootstrap
#'
#' Deprecated wrapper for [test_sdpd_model()].
#'
#' @param res.fit Fitted SDP-D model result object.
#' @param px Vector. Spatial-unit identifiers.
#' @param model Optional SDP-D model object.
#' @param n.boot Integer. Number of bootstrap replications.
#' @param H0 Character or integer. Null hypothesis to test.
#' @param coeff.H0 Numeric scalar or matrix. Coefficient value under the null.
#' @param group.index Vector. Group index.
#' @param opts.boot List. Bootstrap options.
#'
#' @return See [test_sdpd_model()].
#' @seealso [test_sdpd_model()]
#' @export
test.sdpd.model <- function(res.fit,
                            px,
                            model = NULL,
                            n.boot = 399,
                            H0 = c(
                              "zero",
                              "constant",
                              "grouped",
                              "nospatial",
                              "noautoregressive",
                              "noX",
                              "constrained"
                            )[1],
                            coeff.H0 = 0,
                            group.index = rep(1, dim(res.fit$fitted)[1]),
                            opts.boot = list(
                              markovian = FALSE,
                              resid = c("fitted", "normal")[1],
                              sigma.resid = 1,
                              boot.plot = FALSE,
                              cartel = "",
                              ylimiti = NULL,
                              label.index = NULL
                            )) {
  .Deprecated("test_sdpd_model")

  # Map old result field names to the new snake_case API.
  if (!is.null(res.fit$coeff.hat) && is.null(res.fit$coeff_hat)) {
    res.fit$coeff_hat <- res.fit$coeff.hat
  }

  if (!is.null(res.fit$data$W) && is.null(res.fit$data$ww)) {
    res.fit$data$ww <- res.fit$data$W
  }

  if (!is.null(res.fit$data$X) && is.null(res.fit$data$xx)) {
    res.fit$data$xx <- res.fit$data$X
  }

  if (!is.null(res.fit$data$px.neighbors) &&
      is.null(res.fit$data$px_neighbors)) {
    res.fit$data$px_neighbors <- res.fit$data$px.neighbors
  }

  # Map old bootstrap-option names to the new snake_case API.
  if (!is.null(opts.boot$sigma.resid) && is.null(opts.boot$sigma_resid)) {
    opts.boot$sigma_resid <- opts.boot$sigma.resid
  }

  if (!is.null(opts.boot$boot.plot) && is.null(opts.boot$boot_plot)) {
    opts.boot$boot_plot <- opts.boot$boot.plot
  }

  if (!is.null(opts.boot$cartel) && is.null(opts.boot$folder)) {
    opts.boot$folder <- opts.boot$cartel
  }

  if (!is.null(opts.boot$ylimiti) && is.null(opts.boot$y_limits)) {
    opts.boot$y_limits <- opts.boot$ylimiti
  }

  if (!is.null(opts.boot$label.index) && is.null(opts.boot$label_index)) {
    opts.boot$label_index <- opts.boot$label.index
  }

  test_sdpd_model(
    res_fit = res.fit,
    px = px,
    model = model,
    n_boot = n.boot,
    h0 = H0,
    coeff_h0 = coeff.H0,
    group_index = group.index,
    boot_options = opts.boot
  )
}






