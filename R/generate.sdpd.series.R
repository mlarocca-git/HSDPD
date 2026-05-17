generate.sdpd.series <- function(nn, rry, rrxx, model, markovian, num.steps, covariates.sim.model){
  
  ## This function simulates a multivariate (or spatio-temporal) time series from an SDPD model;
  if(is.null(rry) & (is.null(nn) | is.null(model$coeffs))){
    return(list(error="In assenza di rry, per simulare i dati, è necessario conoscere la lunghezza della serie e la matrice di coefficienti"))
  }
  if(is.null(nn))
    nn <- dim(rry)[2]
  if(is.null(rry)){
    pp <- dim(model$coeffs)[1]
    rry <- matrix(0, nrow=pp, ncol=nn)
    dimnames(rry)[[1]] <- dimnames(model$coeffs)[[1]]
    dimnames(rry)[[2]] <- seq(1, nn)
  }
  px <- dimnames(rry)[[1]]
  if(is.null(rrxx)){
    kk <- sum(model$beta_coeffs)
    beta.names <- names(model$beta_coeffs)[model$beta_coeffs]
  }
  data <- list(pp=pp, nn=nn, kk=kk)

  ## generation of errors
  if(is.null(model$eps)){
    if(is.null(model$sigma.eps))
      model$sigma.eps <- matrix(rep(1, data$pp), ncol=1)
    else
      model$sigma.eps <- as.matrix(model$sigma.eps, ncol=1)
    series.eps <- t(apply(model$sigma.eps, 1, FUN=function(xx, nn){rnorm(nn, sd=xx)}, nn=data$nn))
  }
  else if(dim(model$eps)[1]!=data$pp | dim(model$eps)[2]!=data$nn)
    return(list(error="La dimensione della matrice di errori in model$eps non è compatibile con la matrice dati"))
  else series.eps <- model$eps
  if(is.null(dimnames(series.eps)[[1]]))
    dimnames(series.eps)[[1]] <- dimnames(rry)[[1]]
  
  ## generation of the exogenous regressors as stationary ARMA process
  if(data$kk>0){
    if(is.null(rrxx)){
      if(data$kk==1){
        if(beta.names[1]=="trend")
          XX <- matrix(rep(seq(1, data$nn)/data$nn, data$pp), byrow=TRUE, nrow=data$pp)
        else
          XX <- matrix(arima.sim(n = data$nn*data$pp, model=covariates.sim.model), ncol=data$nn, nrow=data$pp, byrow=T)
        dimnames(XX)[[1]] <- dimnames(rry)[[1]]
        dimnames(XX)[[2]] <- dimnames(rry)[[2]]
      }
      else if(data$kk>1){
        XX <- array(0, dim=c(data$kk, data$pp, data$nn))
        for(jj in 1:data$kk){
          if(beta.names[jj]=="trend")
            XX[jj,,] <- matrix(rep(seq(1, data$nn)/data$nn, data$pp), byrow=TRUE, nrow=data$pp)
          else
            XX[jj,,] <- matrix(arima.sim(n = data$nn*data$pp, model=covariates.sim.model), ncol=data$nn, nrow=data$pp, byrow=T)
        }
        dimnames(XX)[[1]] <- beta.names
        dimnames(XX)[[2]] <- dimnames(rry)[[1]]
        dimnames(XX)[[3]] <- dimnames(rry)[[2]]
      }
    }
    else if((is.matrix(rrxx)|is.data.frame(rrxx)) & data$kk==1 & dim(rrxx)[1]==data$nn & dim(rrxx)[2]==data$pp & sum(is.na(rrxx))==0){
      XX 	<- as.matrix(rrxx)
    }
    else if(is.array(rrxx) & data$kk>1 & dim(rrxx)[1]==data$kk & dim(rrxx)[2]==nn & dim(rrxx)[3]==data$pp & sum(is.na(rrxx))==0){
      XX	<- rrxx
    }
    else return(error="Something wrong with rrxx")
  }
  else XX <- NULL
  
  ## checking validity of data...
  data <- check.sdpd.series(series=rry, XX=XX, model=model)

  ## building spatial weights
  WW <- build.spatial.matrix(ww.index=data$ww.index, ww.values=data$ww.values)
  
  rry <- fit.sdpd.series(dseries=data$series, W=WW, X.centr=data$XX, num.steps=num.steps, 
                  model=model, markovian=markovian, coeff.hat=model$coeffs, resid=series.eps)$series

  ## results....  
  list(rry=rry, rrxx=XX, eps=series.eps)
}


