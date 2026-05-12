check.sdpd.model <- function(res.fit=NULL, WW=NULL, coeffs=NULL, model=NULL, mus=NULL){

  vec.error <- vec.warning <- character(20)
  n.error <- n.warning <- 0
  diagnostics <- NULL
  nn <- pp <- kk <- mu <- NULL
  model.obj <- model

  if(!is.null(res.fit)){
    if(is.null(WW)) WW <- res.fit$data$W
    if(is.null(model)) model.obj <- res.fit$model
    if(is.null(coeffs)) coeffs <- res.fit$coeff.hat
    if(is.null(mus)) mu <- apply(res.fit$data$series, 1, mean)
    nn <- dim(res.fit$fitted)[2]
  }
  else if(!is.null(model)){
    coeffs <- model$coeffs
    model.obj <- model
    WW <- model$WW
    ## checking validity of the model
    if(is.null(model.obj$lambda_coeffs)|length(model.obj$lambda_coeffs)>3|length(model.obj$lambda_coeffs)==0|sum(model.obj$lambda_coeffs)==0|(!is.logical(model.obj$lambda_coeffs))){
      n.error <- n.error+1
      vec.error[n.error] <- "There are problems with the lambda-coefficients in the model to be estimated."
    }
    else if(length(model.obj$lambda_coeffs)<3){
      temp <- rep(FALSE, 3)
      names(temp) <- c("lambda0", "lambda1", "lambda2")
      if(is.null(names(model.obj$lambda_coeffs))){
        temp[1:length(model.obj$lambda_coeffs)] <- model.obj$lambda_coeffs
        model.obj$lambda_coeffs <- temp
        n.warning <- n.warning+1
        vec.warning[n.warning] <- "Some of the lambda-coefficients were missing, please check if now it is OK."
      }
      else{
        temp[names(model.obj$lambda_coeffs)] <- model.obj$lambda_coeffs
        model.obj$lambda_coeffs <- temp
        n.warning <- n.warning+1
        vec.warning[n.warning] <- "The names of some lambda-coefficients have been added, please check if now it is OK."
      }
    }
    if(is.null(names(model.obj$lambda_coeffs)))
      names(model.obj$lambda_coeffs) <- c("lambda0", "lambda1", "lambda2")
    
    if(is.null(model.obj$beta_coeffs)|length(model.obj$beta_coeffs)==0|sum(model.obj$beta_coeffs)==0|(!is.logical(model.obj$beta_coeffs))){
      n.warning <- n.warning+1
      vec.warning[n.warning] <- "The beta coefficients have been removed due to some problems, please check the model."
      model.obj$beta_coeffs <- FALSE
    }
    
    if(is.null(model.obj$fixed_effects) | (!is.logical(model.obj$fixed_effects)) | (!model.obj$fixed_effects)){
      model.obj$fixed_effects <- FALSE
      n.warning <- n.warning+1
      vec.warning[n.warning] <- "The model has not the fixed effects...please check if this is OK."
    }
    
    if(is.null(model.obj$time_effects)| !is.logical(model.obj$time_effects) | !(model.obj$time_effects))
      model.obj$time_effects <- FALSE
  }
  else{
    n.error <- n.error + 1
    vec.error[n.error] <- "There is no SDPD model to check"
  }
    
  kk <- sum(model.obj$beta_coeffs)
  
  if(!is.null(coeffs)){
    pp <- dim(coeffs)[1]
    ## calcoliamo autovalori delle componenti della matrice ridotta A
    ll0 <- ll1 <- ll2 <- matrix(0, ncol=pp, nrow=pp)
    vettore1 <- vettore2 <- eigenA <- rep(0, pp)
    if("lambda0"%in%dimnames(coeffs)[[2]]){
      ll0 <- diag(coeffs[,"lambda0"])%*%WW
      matrice1 <- solve(diag(rep(1, pp))-ll0)
      if(sum(is.na(matrice1))==0)
        vettore1 <- eigen(matrice1)$values
      else 
        vettore1 <- rep(NA, pp)
    }
    else matrice1 <- diag(vettore1)

    if("lambda1"%in%dimnames(coeffs)[[2]])
      ll1 <- diag(coeffs[,"lambda1"])

    if("lambda2"%in%dimnames(coeffs)[[2]])
      ll2 <- diag(coeffs[,"lambda2"])%*%WW

    matrice2 <- ll1+ll2
    if(sum(is.na(matrice2))==0)
      vettore2 <- eigen(matrice2)$values
    else
      vettore2 <- rep(NA, pp)
    if(sum(is.na(matrice1%*%matrice2))==0)
      eigenA <- eigen(matrice1%*%matrice2)$values
    else
      eigenA <- rep(NA, pp)
    diagnostics <- cbind(eigen1=Mod(vettore1), eigen2=Mod(vettore2), Mod.eigenA=Mod(eigenA))
    dimnames(diagnostics)[[1]] <- dimnames(coeffs)[[1]]
  }
  ## restituzione risultati
  list(model=model.obj, diagnostics=data.frame(diagnostics), coeffs=coeffs,
       errors=vec.error[vec.error!=""], warnings=vec.warning[vec.warning!=""], nn=nn, pp=pp, kk=kk, mu=mu)
}
