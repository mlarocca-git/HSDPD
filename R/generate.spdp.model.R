
generate.sdpd.model <- function(pp, model, ww.index, ww.values, px.neighbors, sim.options)
{
  ## generation of the (positive definite) weight matrix WW
  if(is.null(ww.index) | is.null(ww.values)){
    flag.WW 	<- TRUE
    while(flag.WW) {
      if(sim.options$ww.type==1 | sim.options$ww.type=="rook"){
        WW <- ww.index <- matrix(0, ncol=5, nrow=pp)
        ww.index[,1] <- seq(1, pp)
        ww.index[,2] <- ifelse(seq(1, pp)+1>pp, 0, seq(1, pp)+1)
        ww.index[,3] <- ifelse(seq(1, pp)+2>pp, 0, seq(1, pp)+2)
        ww.index[,4] <- ifelse(seq(1, pp)-1<1, 0, seq(1, pp)-1)
        ww.index[,5] <- ifelse(seq(1, pp)-2<1, 0, seq(1, pp)-2)
        WW[,1] <- 0
        WW[,2:5] <- 1
        flag.WW 	<- FALSE
      }
      else if(sim.options$ww.type==2 | sim.options$ww.type=="queen"){
        WW <- ww.index <- matrix(0, ncol=9, nrow=pp)
        ww.index[,1] <- seq(1, pp)
        ww.index[,2] <- ifelse(seq(1, pp)+1>pp, 0, seq(1, pp)+1)
        ww.index[,3] <- ifelse(seq(1, pp)+2>pp, 0, seq(1, pp)+2)
        ww.index[,4] <- ifelse(seq(1, pp)-1<1, 0, seq(1, pp)-1)
        ww.index[,5] <- ifelse(seq(1, pp)-2<1, 0, seq(1, pp)-2)
        ww.index[,6] <- ifelse(seq(1, pp)+3>pp, 0, seq(1, pp)+3)
        ww.index[,7] <- ifelse(seq(1, pp)+4>pp, 0, seq(1, pp)+4)
        ww.index[,8] <- ifelse(seq(1, pp)-3<1, 0, seq(1, pp)-3)
        ww.index[,9] <- ifelse(seq(1, pp)-4<1, 0, seq(1, pp)-4)
        WW[,1] <- 0
        WW[,2:9] <- 1
        flag.WW 	<- FALSE
      }
      else if(sim.options$ww.type==3 | sim.options$ww.type=="corr"){
        WW <- ww.index <- matrix(rnorm(pp * pp, sd=1), ncol = pp)
        WW <- WW %*% t(WW)
        diag(WW) <- 0 
        flag.WW <- (det(WW) == 0)
        ww.index <- matrix(rep(seq(1,pp), pp), nrow = pp, byrow = TRUE)
      }
    }
    dimnames(ww.index)[[1]] <- dimnames(WW)[[1]] <- seq(1, pp)
  }
  
  ## row-normalization of the spatial matrix by L1 norm
  ww.values 	<- t(apply(WW, 1, FUN = function(x) {x/sum(abs(x))}))

  ## generation of the matrix with neighbor points
  if(is.null(px.neighbors)){
    n.neigh.cols <- (2*sim.options$px.neighbors+1)^2-1
    neighbors <- matrix(NA, nrow=pp, ncol=n.neigh.cols+1)
    neighbors[,1:min(n.neigh.cols, dim(ww.index)[2])] <- ww.index[,1:min(n.neigh.cols, dim(ww.index)[2])]

    individua.neighbors <- function(x, ww.i, n.col){
      res <- rep(NA, n.col)
      vicini <- x[abs(x)>0][-1]
      vicini <- na.exclude(vicini)
      xx <- as.vector(ww.i[as.character(vicini),])
      vicini.nuovi <- xx[abs(xx)>0]
      vicini.nuovi <- unique(vicini.nuovi)
      vicini.nuovi <- vicini.nuovi[!(vicini.nuovi%in%c(x[1],vicini))]
      vicini.nuovi <- c(vicini, vicini.nuovi)
      if(length(vicini.nuovi)>=n.col)
        res[1:n.col] <- vicini.nuovi[1:n.col]
      else
        res[1:length(vicini.nuovi)] <- vicini.nuovi[1:length(vicini.nuovi)]
      ## output
      c(x[1], res)
    }
    # flag <- FALSE
    # while(!flag){
    for(zz in 1:100){
      old.neighbors <- neighbors
      neighbors <- t(apply(old.neighbors, 1, FUN=individua.neighbors, ww.i=ww.index, n.col=n.neigh.cols))
      # flag <- sum(abs(old.neighbors-neighbors), na.rm=TRUE)<1
    }
    neighbors <- neighbors[,-1]
    dimnames(neighbors)[[1]] <- dimnames(ww.index)[[1]]
  }
    
  ## generation of model parameters lambda and beta
  lambda.names <- names(model$lambda_coeffs)[model$lambda_coeffs]
  beta.names <- names(model$beta_coeffs)[model$beta_coeffs]
  fixed.name <- c("fixed_effects")[model$fixed_effects]
  param 	<- matrix(0, nrow=pp, ncol=sum(model$lambda_coeffs)+sum(model$beta_coeffs)+sum(model$fixed_effects))
  dimnames(param)[[2]] <- c(lambda.names, beta.names, fixed.name)
  
  ## fixed effects values
  if(model$fixed_effects){
    if(sim.options$fixed.effects$constant[1])
      param[,"fixed_effects"]	<- rep(runif(1, sim.options$fixed.effects$min[1], sim.options$fixed.effects$max[1]), pp)
    else 
      param[,"fixed_effects"]	<- runif(pp, sim.options$fixed.effects$min[1], sim.options$fixed.effects$max[1])
  }
  
  
  ## setting some default parameters
  if(is.null(sim.options$index.weights))
    sim.options$index.weights <- rep(1, pp)
  
  ## beta coefficients for the exogeneous regressors
  kk <- length(beta.names)
  if(kk>0){
    if(is.null(sim.options$betas$constant))
      sim.options$betas$constant <- rep(F, kk)
    if(length(sim.options$betas$constant)==1 | length(sim.options$betas$constant)!=kk)
      sim.options$betas$constant <- rep(sim.options$betas$constant[1], kk)
    if(length(sim.options$betas$min)==1 | length(sim.options$betas$min)!=kk)
      sim.options$betas$min <- rep(sim.options$betas$min[1], kk)
    if(length(sim.options$betas$max)==1 | length(sim.options$betas$max)!=kk)
      sim.options$betas$max <- rep(sim.options$betas$max[1], kk)
    for(ss in 1:kk){
      if(sim.options$betas$constant[ss]){
        for(jj in 1:max(sim.options$index.weights))
          param[sim.options$index.weights==jj,beta.names[ss]] <- rep(round((runif(1, min=sim.options$betas$min, max=sim.options$betas$max)), digits = 4), sum(sim.options$index.weights==jj))
      }
      else
        param[,beta.names[ss]] <- (round(runif(pp, min=sim.options$betas$min[ss], max=sim.options$betas$max[ss]), digits = 4))
    }
  }
  
  ## generation of lambda coefficients
  if(model$lambda_coeffs["lambda0"] & is.null(sim.options$lambda0$constant)) sim.options$lambda0$constant <- FALSE
  if(model$lambda_coeffs["lambda1"] & is.null(sim.options$lambda1$constant)) sim.options$lambda1$constant <- FALSE
  if(model$lambda_coeffs["lambda2"] & is.null(sim.options$lambda2$constant)) sim.options$lambda2$constant <- FALSE
  if(sim.options$lambda0$constant & model$lambda_coeffs["lambda0"]){
    for(jj in 1:max(sim.options$index.weights))
      param[sim.options$index.weights==jj,"lambda0"] <- rep(round((runif(1, min=sim.options$lambda0$min, max=sim.options$lambda0$max)), digits = 4), sum(sim.options$index.weights==jj))
  }
  else
    param[,"lambda0"] <- (round(runif(pp, min=sim.options$lambda0$min, max=sim.options$lambda0$max), digits = 4))
  if(sim.options$lambda1$constant){
    for(jj in 1:max(sim.options$index.weights))
      param[sim.options$index.weights==jj,"lambda1"] <- rep(round(runif(1, min=sim.options$lambda1$min, max=sim.options$lambda1$max), digits = 4), sum(sim.options$index.weights==jj))
  }
  else
    param[,"lambda1"] <- (round(runif(pp, sim.options$lambda1$min, sim.options$lambda1$max), digits = 4))
  if(sim.options$lambda2$constant){
    for(jj in 1:max(sim.options$index.weights))
      param[sim.options$index.weights==jj,"lambda2"] <- rep(round((runif(1, min=sim.options$lambda2$min, max=sim.options$lambda2$max)), digits = 4), sum(sim.options$index.weights==jj))
  }
  else 
    param[,"lambda2"] <- (round(runif(pp, sim.options$lambda2$min, sim.options$lambda2$max), digits = 4))
  
  ## generation of sigma.eps values
  sigma.eps <- numeric(pp)
  if(is.null(sim.options$sigma.eps$constant))
    sim.options$sigma.eps$constant <- rep(T, kk)
  if(length(sim.options$sigma.eps$constant)==1 | length(sim.options$sigma.eps$constant)!=kk)
    sim.options$sigma.eps$constant <- rep(sim.options$sigma.eps$constant[1], kk)
  if(length(sim.options$sigma.eps$min)==1 | length(sim.options$sigma.eps$min)!=kk)
    sim.options$sigma.eps$min <- rep(sim.options$sigma.eps$min[1], kk)
  if(length(sim.options$sigma.eps$max)==1 | length(sim.options$sigma.eps$max)!=kk)
    sim.options$sigma.eps$max <- rep(sim.options$sigma.eps$max[1], kk)
  for(ss in 1:kk){
    if(sim.options$sigma.eps$constant[ss]){
      for(jj in 1:max(sim.options$index.weights))
        sigma.eps[sim.options$index.weights==jj] <- rep(round((runif(1, min=sim.options$sigma.eps$min[ss], max=sim.options$sigma.eps$max[ss])), digits = 4), sum(sim.options$index.weights==jj))
    }
    else
      sigma.eps[sim.options$index.weights==jj] <- (round(runif(pp, min=sim.options$sigma.eps$min[ss], max=sim.options$sigma.eps$max[ss]), digits = 4))
  }

  ### returning the parameters of the model
  dimnames(param)[[1]] <- names(sigma.eps) <- names(sim.options$index.weights) <- seq(1, pp)
  list(coeffs = param, ww.values = ww.values, ww.index=ww.index, px.neighbors=list(index=neighbors, seriesBoundary=NULL),
       sigma.eps=sigma.eps, index.weights=sim.options$index.weights)
}

