fit.sdpd.mean.equation.model <- function(res, WW, time.weights=NULL, index.weights=NULL){
  
  ## computing the mean over the whole period
  mu.i <- apply(res$fitted, 1, mean, na.rm=T)
  mu.around.i <- as.vector(WW%*%mu.i)
  mus <- cbind(mu.i=mu.i, mu.around.i=mu.around.i)
  if(!is.null(index.weights)){
    pesi.R <- as.vector(WW%*%index.weights)
    muR.around.i <- as.vector(WW%*%(mu.i*index.weights))/pesi.R
    mus <- cbind(mus, mu_j_weighted.around.i=muR.around.i)
  }
  
  ## computing the means over time groups (for example, for seasons) or by time-weighting
  mu.group.i <- mu.group.around.i <- muR.group.around.i <- mus_t_weighted <- NULL
  if(!is.null(time.weights)){
    if(is.factor(time.weights)){
      groups <- levels(time.weights)
      for(group in groups){
        mu.temp <- apply(res$fitted[,time.weights==group], 1, mean, na.rm=T)
        mu.group.i <- cbind(mu.group.i, mu.temp)
        mu.group.around.i <- cbind(mu.group.around.i, as.vector(WW%*%mu.temp))
        if(!is.null(index.weights))
          muR.group.around.i <- cbind(muR.group.around.i, as.vector(WW%*%(mu.temp*index.weights))/pesi.R)
      }
      dimnames(mu.group.i)[[2]] <- dimnames(mu.group.around.i)[[2]] <- groups
      #mus_t_weighted$mu.i <- data.frame(mu.group.i)
      #mus_t_weighted$mu.around.i <- data.frame(mu.group.around.i)
      if(!is.null(index.weights)){
        dimnames(muR.group.around.i)[[2]] <- groups
        #mus_t_weighted$mu_j_weighted.around.i <- data.frame(muR.group.around.i)
      }
      mus_t_weighted <- data.frame(mu.i=data.frame(mu.group.i), mu.around.i=data.frame(mu.group.around.i),
                                  mu_j_weighted.around.i=data.frame(muR.group.around.i))
    }
    else if(is.numeric(time.weights)){
      pesi.T <- sum(time.weights)
      mus_t_weighted$mu.i <- apply(res$fitted, 1, FUN=function(x,w){mean(x*w, na.rm=T)}, w=time.weights)
      mus_t_weighted$mu.around.i <- as.vector(WW%*%mus_t_weighted$mu.i)
      if(!is.null(index.weights))
        mus_t_weighted$mu_j_weighted.around.i <- as.vector(WW%*%(mus_t_weighted$mu.i*index.weights))/pesi.R
    }
  }

  ## returning results
  list(mus=data.frame(mus), mus_t_weighted=mus_t_weighted)
}

