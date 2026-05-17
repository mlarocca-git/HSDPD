fit.sdpd.procedure <- function(series, model, check=FALSE, two.stage=FALSE, NAcovs="pairwise.complete.obs") 
{
  ## series can be an object of class sdpd-series (from function build.sdpd.series) or it can be the object dseries, where
  ## dseries is a matrix of dimension (pp, nn) with pp the number of univariate time dseries and nn the number of time observations;

  if(is.null(series)){
    cat("\n Errore: data is NULL")
    return(NULL)
  }

  ## checking validity of data...
  data <- check.sdpd.series(series=series$series, ww.index=series$ww.index, ww.values=series$ww.values, XX=series$X, model=model, 
                            px.neighbors=series$px.neighbors, px=series$px, lat=series$lat, lon=series$lon, group=series$group,
                            index.weights=series$index.weights, time.weights=series$time.weights)
    
  if(length(data$errors)>0){
    cat("\n There are some errors with the series....\n", data$errors)
    return(list(errors=data$errors, warnings=data$warnings))
  }
    
  COVs <- fit.sdpd.covs(series=data$series, X=data$XX, px.neighbors=data$px.neighbors, kk=data$kk, nn=data$nn, pp=data$pp, NAcovs=NAcovs)

  ## building spatial weights
  WW <- build.spatial.matrix(ww.index=data$ww.index, ww.values=data$ww.values)

  ## estimating model parameters...
  fit <- fit.sdpd.coefficients(W=WW, COVs=COVs, mu=data$mu, model=model)
  if(sum(is.na(fit$coeff.hat))>0){
    data$warnings <- c(data$warnings, "There are some NA in estimated coefficients.")
    cat("\n There are some NA in estimated coefficients.\n")
  }
    
  ## second stage estimation...
  if(two.stage)
    fit$coeff.hat <- fit.2nd.stage(dseries=data$series, W=WW, X.centr=data$X, model=model, coeff.hat=fit$coeff.hat)
    
  ## estimating fitted values...
  res <- fit.sdpd.series(dseries=data$series, W=WW, X.centr=data$X, model=model, coeff.hat=fit$coeff.hat, time_effects=data$time_effects)
    
  ## estimating the mean equation meta-model
  mus <- fit.sdpd.mean.equation.model(res=res, WW=WW, time.weights=data$time.weights, index.weights=data$index.weights)
    
  ## checking stationarity conditions for the estimated model
  res <- list(px=data$px, lon=data$lon, lat=data$lat, group=data$group, coeff.hat=fit$coeff.hat, fitted=res$fitted, resid=res$resid,
              mean.equation=mus, time_effects=data$time_effects, model=model, warnings=data$warning)
  if(check){
    check <- check.sdpd.model(res.fit=res, ww.index=data$ww.index, ww.values=data$ww.values)$diagnostics[1,"max.mod.eigenA"]
  }

  ## output
  res$diagnostics=check
  res
}

