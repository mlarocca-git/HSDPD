fit.sdpd.model <- function(series=NULL, W=NULL, X=NULL, data.back=FALSE, 
                           model,
                           index.weights=NULL, time.weights=NULL,
                           two.stage=FALSE,
                           NAcovs="pairwise.complete.obs") 
{
  ## dseries is a matrix of dimension (nn,pp) where pp is the number of univariate time dseries and data$nn the number of time observations
  ## W is a spatial weight matrix of dimension (pp,pp)
  ## X is an array of dim=c(kk,nn,pp) which includes the data for kk exogenous regressors. If kk=1 then X is a matrix of dim=(nn,pp)

  if(is.null(series)){
    cat("\n Errore: df.data is NULL")
    return(NULL)
  }

  ## checking validity of data...
  data <- check.sdpd.series(series=series, WW=W, XX=X, model=model, index.weights=index.weights, time.weights=time.weights)

  if(length(data$errors)>0){
    cat("\n There are some errors with the series....\n", data$errors)
    return(list(errors=data$errors, warnings=data$warnings))
  }

  COVs <- fit.sdpd.covs(series=data$series, X=data$XX, px.neighbors=data$neighbours, kk=data$kk, nn=data$nn, pp=data$pp, NAcovs=NAcovs)
  
  ## estimating model parameters...
  fit <- fit.sdpd.coefficients(W=data$WW, COVs=COVs, mu=data$mu, model=data$model)
  if(sum(is.na(fit$coeff.hat))>0){
    data$warnings <- c(data$warnings, "There are some NA in estimated coefficients.")
    cat("\n There are some NA in estimated coefficients.\n")
  }

  ## second stage estimation...
  if(two.stage)
    fit$coeff.hat <- fit.2nd.stage(dseries=data$series, W=data$WW, X.centr=data$X, model=data$model, coeff.hat=fit$coeff.hat)

  ## estimating fitted values...
  res <- fit.sdpd.series(dseries=data$series, W=data$WW, X.centr=data$X, model=data$model, coeff.hat=fit$coeff.hat, time_effects=data$time_effects)

  ## estimating the mean equation meta-model
  mus <- fit.sdpd.mean.equation.model(res=res, WW=data$WW, time.weights=data$time.weights, index.weights=data$index.weights)

  ## model diagnostics
  diagnostics <- check.sdpd.model(res.fit=res, coeffs=fit$coeff.hat, WW=data$WW, model=model, mus=mus$mus$mu.i)
  if(length(diagnostics$errors)>0){
    cat("\n There are some errors....", diagnostics$errors)
    return(list(errors=data$errors, warnings=data$warnings))
  }

  
  ## returning estimation results
  if(data.back)
    data.back <- list(series=data$series, W=data$WW, X=data$XX, neighbours=data$neighbours)

  list(px=data$px, lon=data$lon, lat=data$lat, group=data$group, coeff.hat=fit$coeff.hat, fitted=res$fitted, resid=res$resid, mean.equation=mus,
       time_effects=data$time_effects, data=data.back, diagnostics=diagnostics$diagnostics$Mod.eigenA, model=data$model, warnings=data$warning)
}
