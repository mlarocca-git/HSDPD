
fit.sdpd.coefficients <- function (W, COVs, mu, model) 
{
  lambda.names <- names(model$lambda_coeffs)[model$lambda_coeffs]
  beta.names <- names(model$beta_coeffs)[model$beta_coeffs]
  fixed_effects.name <- c("fixed_effects")[model$fixed_effects]
  px <- dimnames(W)[[1]]
  data <- list(WW=W, pp=length(px), kk=length(beta.names))
  
  ## variable definition...
  coeff.hat	<- matrix(0, nrow=data$pp, ncol=sum(model$fixed_effects)+sum(model$lambda_coeffs)+sum(model$beta_coeffs))
  dimnames(coeff.hat)[[2]] <- c(lambda.names, beta.names, fixed_effects.name)
  dimnames(coeff.hat)[[1]] <- px
  ei <- numeric(data$pp); names(ei) <- px
  
  ## estimation of coefficients...
  invertible <- function(m) class(try(solve(m),silent=T))[1]=="matrix"
  if(data$kk==0){
    name.coeff <- lambda.names # sum(model$lambda_coeffs)
    for(ii in px){
      wi <- W[ii,px]
      ei[px] <- 0; ei[ii] <- 1
      indici <- as.character(na.exclude(COVs$index[ii,]))
      Yi <- t(COVs$cov12[px,indici])%*%ei
      Xi <- cbind(t(COVs$cov12[px,indici])%*%wi, t(COVs$cov11[px,indici])%*%ei, t(COVs$cov11[px,indici])%*%wi)[,model$lambda_coeffs]
      if(!invertible(t(Xi)%*%Xi)){
        coeff.hat[ii,name.coeff] <- NA
        next
      }
      coeff.hat[ii,name.coeff] <- solve(t(Xi)%*%Xi)%*%t(Xi)%*%Yi
    }
  }
  else if(data$kk==1){
    name.coeff <- c(lambda.names, beta.names)
    for(ii in px){
      wi <- W[ii,px]
      ei[px] <- 0; ei[ii] <- 1
      indici <- as.character(na.exclude(COVs$index[ii,]))
      Yi <- t(COVs$cov12[px,indici])%*%ei
      Xi <- cbind(t(COVs$cov12[px,indici])%*%wi, t(COVs$cov11[px,indici])%*%ei, t(COVs$cov11[px,indici])%*%wi, t(COVs$covX[px,indici])%*%ei)[,c(model$lambda_coeffs,T)]
      if(!invertible(t(Xi)%*%Xi)){
        coeff.hat[ii,name.coeff] <- NA
        next
      }
      coeff.hat[ii,name.coeff] <- (solve(t(Xi)%*%Xi)%*%t(Xi)%*%Yi)
    }
  }
  else if(data$kk>1){
    name.coeff <- c(lambda.names, beta.names)
    for(ii in px){
      wi <- W[ii,px]
      ei[px] <- 0; ei[ii] <- 1
      indici <- as.character(na.exclude(COVs$index[ii,]))
      Yi <- t(COVs$cov12[px,indici])%*%ei
      Xi <- cbind(t(COVs$cov12[px,indici])%*%wi, t(COVs$cov11[px,indici])%*%ei, t(COVs$cov11[px,indici])%*%wi)[,model$lambda_coeffs]
      Xi <- cbind(Xi, t(COVs$covX[beta.names,ii,indici]))
      if(!invertible(t(Xi)%*%Xi)){
        coeff.hat[ii,name.coeff] <- NA
        next
      }
      coeff.hat[ii,name.coeff] <- (solve(t(Xi)%*%Xi)%*%t(Xi)%*%Yi)
    }
  }
  if(model$fixed_effects){
    ll0 <- ll1 <- ll2 <- llv <- matrix(0, ncol=data$pp, nrow=data$pp)
    if(model$lambda_coeffs[1]) ll0 <- diag(coeff.hat[px,"lambda0"])%*%W
    if(model$lambda_coeffs[2]) ll1 <- diag(coeff.hat[px,"lambda1"])
    if(model$lambda_coeffs[3]) ll2 <- diag(coeff.hat[px,"lambda2"])%*%W
    B <- diag(rep(1, data$pp))-ll0-ll1-ll2
    coeff.hat[px,fixed_effects.name] <- B%*%mu[px]
  }
  
  ## estimation results...
  list(coeff.hat=coeff.hat)
}

