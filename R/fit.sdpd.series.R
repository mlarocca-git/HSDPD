fit.sdpd.series <- function(dseries, W, X.centr=NULL, model, coeff.hat, time_effects=NULL, 
                            px.sim=NULL, resids=NULL, markovian=FALSE, num.steps=1) 
{
  beta.names <- names(model$beta_coeffs)[model$beta_coeffs]
  px <- dimnames(dseries)[[1]]
  tempo <- dimnames(dseries)[[2]]
  data <- list(pp=length(px), nn=length(tempo), kk=length(beta.names))
  
  ## building model components...
  ll0 <- ll1 <- ll2 <- matrix(0, nrow=data$pp, ncol=data$pp)
  Id <- diag(rep(1, data$pp))
  if(model$lambda_coeffs["lambda0"]) ll0 <- diag(coeff.hat[,"lambda0"])
  if(model$lambda_coeffs["lambda1"]) ll1 <- diag(coeff.hat[,"lambda1"])
  if(model$lambda_coeffs["lambda2"]) ll2 <- diag(coeff.hat[,"lambda2"])
  llk <- llc <- matrix(0, nrow=data$pp, ncol=data$nn)
  if(is.null(model$time_effects) | !model$time_effects | is.null(time_effects))
    time_effects <- numeric(data$nn)
  llv <- matrix(rep(time_effects, data$pp), byrow = TRUE, nrow=data$pp)
  if(model$fixed_effects){
    llc <- coeff.hat[,"fixed_effects"]%*%t(rep(1, data$nn))
  }
  if(data$kk==1){
    llk <- diag(coeff.hat[,beta.names[1]])%*%X.centr
  }
  else if(data$kk>1){
    for(jj in 1:data$kk){
      llk <- llk + diag(coeff.hat[,beta.names[jj]])%*%X.centr[jj,,]
    }
  }
  if(is.null(px.sim))
    px.sim <- as.character(px)
  else
    px.sim <- as.character(px.sim)

  ## estimation results...
  if(!is.null(resids) & markovian){
    ## questo rappresenta il caso autoregressivo di tipo markoviano
    fitteds <- newresid <- as.matrix(resids)
    tempI  <- solve(Id-ll0%*%W)
    AA		<- tempI%*%(ll1 + ll2%*%W)
    eps.star1 	<- tempI%*%(llk + llc + newresid)
    for(step in 1:num.steps){
      if(step>1)
        dseries[,1] <- fitteds[,2]+eps.star1[px.sim,1]
      for(tt in 2:data$nn){
        fitteds[,tt] <- AA%*%(dseries)[,tt-1]
        dseries[px.sim,tt] <- fitteds[px.sim,tt] + eps.star1[px.sim,tt]
      }
    }
    fitteds <- dseries - newresid
  }
  else if(!is.null(resids)){
    ## questo rappresenta il caso autoregressivo di tipo punto-fisso
    fitteds <- newresid <- as.matrix(resids)
    AA <- ll0%*%W
    BB <- (ll1+ll2%*%W)
    eps.star1 <- llk+llc+llv + newresid
    dimnames(eps.star1) <- dimnames(dseries)
    for(step in 1:num.steps){
      if(step>1)
        dseries[,1] <- fitteds[,2]
      for(tt in 2:data$nn){
        fitteds[,tt] <- AA%*%dseries[,tt]+BB%*%dseries[,tt-1]
        dseries[px.sim,tt] <- fitteds[px.sim,tt]+eps.star1[px.sim,tt]
      }
    }
    fitteds <- dseries - newresid
  }
  else{
    ## questo rappresenta il caso classico di calcolo dei fitted values
    fitteds <- as.matrix(dseries)
    fitteds[,2:data$nn] <- ll0%*%W%*%dseries[,2:data$nn]+(ll1+ll2%*%W)%*%dseries[,-data$nn]+llk[,-1]+llc[,-1]+llv[,-1]
    fitteds[,1] <- rep(NA, data$pp)
    newresid	<- dseries - fitteds
  }
  
  ## return...
  list(series=dseries, fitted=fitteds, resid=newresid)
}
