fit.2nd.stage <- function(dseries, X.centr, W, model, coeff.hat, data){
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
  ## second stage estimator for beta coefficients
  beta.names <- names(model$beta_coeffs)[model$beta_coeffs]
  YY2	<- dseries[,2:data$nn] - ll0%*%W%*%dseries[,2:data$nn] - ll1%*%dseries[,1:(data$nn-1)] - ll2%*%W%*%dseries[,1:(data$nn-1)]
  if(model$fixed_effects){
    if(data$kk>1){
      for(ii in 1:data$pp){
        Xi <- cbind(t(X.centr[beta.names,ii,-1]), rep(1, data$nn-1))
        coeff.hat[ii,c(beta.names, "fixed_effects")] <- (solve(t(Xi)%*%Xi)%*%t(Xi)%*%(YY2[ii,]))
      }
    }
    else{
      for(ii in 1:data$pp){
        Xi <- cbind(t(X.centr[ii,-1]), rep(1, data$nn-1))
        coeff.hat[ii,c(beta.names, "fixed_effects")] <- (solve(t(Xi)%*%Xi)%*%t(Xi)%*%(YY2[ii,]))
      }
    }
  }
  else if(data$kk>1){
    for(ii in 1:data$pp){
      Xi <- t(X.centr[beta.names,ii,-1])
      coeff.hat[ii,beta.names] <- (solve(t(Xi)%*%Xi)%*%t(Xi)%*%(YY2[ii,]))
    }
  }
  else if(data$kk==1){
    for(ii in 1:data$pp){
      Xi <- t(X.centr[ii,-1])
      coeff.hat[ii,beta.names] <- (solve(t(Xi)%*%Xi)%*%t(Xi)%*%(YY2[ii,]))
    }
  }
  coeff.hat
}

