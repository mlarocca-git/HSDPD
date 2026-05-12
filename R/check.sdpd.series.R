check.sdpd.series <- function(series, WW=NULL, XX=NULL, model, index.weights=NULL, time.weights=NULL, vec.options=NULL) 
{
	## This function receives a multivariate (spatio-temporal) time series and a spatial weight matrix and some regressors,
	## then it verifies if the structure of the dataset is correct
	## series is a matrix of dimension (pp, nn) where pp is the number of univariate time series and nn the number of time observations
	## WW is a spatial weight matrix of dimension (pp,pp)
	## XX is an array of dim=c(kk, pp, nn) which includes the values for kk exogenous regressors. If kk=1 then XX is a matrix of dim=(nn,pp)

	vec.error <- vec.warning <- character(20)
	n.error <- n.warning <- 0
	
	## checking validity of series
	if(is.list(series)){
	  dseries <- series$series
	  px <- series$px
	  lon <- series$lon
	  lat <- series$lat
	  group <- series$group
    if(is.null(series$px.neighbors)){
      n.error <- n.error + 1
      vec.error[n.error] <- "\n The object px.neighbors is NULL!"
    }
	  px.neighbors <- series$px.neighbors
	  if(is.null(XX))
	    XX <- series$X
	  if(is.null(WW))
	    WW <- series$W
    if(is.null(WW))
      WW <- list(ww.index=series$ww.index, ww.values=series$ww.values)
  }
	else
	  dseries <- series

	if(is.matrix(dseries)|is.data.frame(dseries)){
	  dseries <- as.matrix(dseries)
	  nn <- dim(dseries)[2]
		pp <- dim(dseries)[1]
		tempo <- dimnames(dseries)[[2]]
		if(is.null(dimnames(dseries)[[1]])){
		  ## note that names of pixels should be numbers, so that they can be also managed by raster objects
		  dimnames(dseries)[[1]] <- seq(1, pp)
		}
		if(is.null(dimnames(dseries)[[2]])){
		  ## note that names of time units should be numbers (or dates), so that they can converted in Date-format
		  dimnames(dseries)[[2]] <- tempo <- seq(1, nn)
		}
		if(!is.null(time.weights)){
      if(length(time.weights)!=nn){
        time.weights <- NULL
        n.warning <- n.warning + 1
        vec.warning[n.warning] <- "\n There is a problem with time.weights...they been set to NULL"
      }	  
		}
		if(!is.null(index.weights)){
		  if((is.matrix(index.weights) | is.data.frame(index.weights)) & dim(index.weights)[2]>1){
		    ix.weights <-  index.weights[dimnames(dseries)[[1]],2]
		    names(ix.weights) <- index.weights[dimnames(dseries)[[1]],1]
		    index.weights <- ix.weights
		  }
		  else{
		    index.weights <- NULL
		    n.warning <- n.warning + 1
		    vec.warning[n.warning] <- "\n There is a problem with argument index.weights...it has been set to NULL"
		  }	  
		}
	}
	else{
	  n.error <- n.error + 1
	  vec.error[n.error] <- "The series is not a matrix or dataframe"
	}

	## checking validity of regressors
	if(is.list(model)){
	  model.obj <- model
	  nomi.covariate <- names(model$beta_coeffs)[model$beta_coeffs]
	}
	else{
	  n.error <- n.error + 1
	  vec.error[n.error] <- "Something wrong with the model, which is not a list. Please check."
	}
	if(is.null(model.obj$kk))
	  model.obj$kk <- sum(model.obj$beta_coeffs)
	else if(model.obj$kk!=sum(model.obj$beta_coeffs)){
	  n.warning <- n.warning + 1
	  vec.warning[n.warning] <- "The passed model$kk value is not coerent with model$beta_coeffs. It has been corrected, please check."
	}
	if(is.null(XX)){
	  if(model.obj$kk>0){
	    n.error <- n.error + 1
	    vec.error[n.error] <- paste("The model has exogenous coavariates, but there is no data passed in XX." )
	  }
	  kk <- 0
	}
	else if(model.obj$kk==0){
	  n.warning <- n.warning + 1
	  vec.warning[n.warning] <- "The model has no exogenous covariates, so the object passed in X has been ignored."
	  
	}
	else if(is.matrix(XX)|is.data.frame(XX)){
	  if(dim(XX)[2]!=nn){
	    n.error <- n.error + 1
	    vec.error[n.error] <- "The regressor must have the same number of observations (=columns) as the series."
	  }
	  if(dim(XX)[1]!=pp){
	    n.error <- n.error + 1
	    vec.error[n.error] <- "The regressor must have the same number of locations (=rows) as the series."
	  }
	  XX <- as.matrix(XX)
	  XX <- apply(XX, 1, FUN=function(x){x-mean(x)})
	  XX <- t(XX)
	  n.warning <- n.warning + 1
	  vec.warning[n.warning] <- "The regressor has been mean-centered."
	  if(is.null(dimnames(XX)[[1]]))
	    dimnames(XX)[[1]] <- dimnames(dseries)[[1]]
	  if(is.null(dimnames(XX)[[2]]))
	    dimnames(XX)[[2]] <- dimnames(dseries)[[2]]
	  kk <- 1
	}
	else if(is.array(XX) & length(dim(XX)==3)){
	  if(dim(XX)[3]!=nn){
	    n.error <- n.error + 1
	    vec.error[n.error] <- paste("The regressors must have the same number of observations as the series (=", nn, " instead of ", dim(XX)[3], ").", sep="")
	  }
	  if(dim(XX)[2]!=pp){
	    n.error <- n.error + 1
	    vec.error[n.error] <- paste("The regressors must have the same number of locations as the series (=", pp, " instead of ", dim(XX)[2], ").", sep="")
	  }
	  XX <- apply(XX, c(1,2), FUN=function(x){x-mean(x)})
	  XX <- aperm(XX, c(2,3,1))
	  n.warning <- n.warning + 1
	  vec.warning[n.warning] <- "The regressors have been mean-centered."
	  if(is.null(dimnames(XX)[[3]]))
	    dimnames(XX)[[3]] <- dimnames(dseries)[[2]]
	  if(is.null(dimnames(XX)[[2]]))
	    dimnames(XX)[[2]] <- dimnames(dseries)[[1]]
	  if(is.null(dimnames(XX)[[1]])){
	    if(dim(XX)[1]<length(nomi.covariate)){
	      n.error <- n.error + 1
	      vec.error[n.error] <- "Some covariates of the model are missing in the array X."
	    }
	    else if(dim(XX)[1]>length(nomi.covariate)){
	      n.warning <- n.warning + 1
	      vec.warning[n.warning] <- paste("There are more covariates in the array X than in the model. The first",
	                                  length(nomi.covariate), "have been used.")
	      XX <- XX[1:length(nomi.covariate),,]
	      dimnames(XX)[[1]] <- nomi.covariate
	    }
	    else{
	      n.warning <- n.warning + 1
	      vec.warning[n.warning] <- paste("The array X has no named covariates. They have been used in the same order as in the passed model.")
	      dimnames(XX)[[1]] <- nomi.covariate
	    }
	  }
	  kk <- dim(XX)[1]
	}
	else{
	  n.error <- n.error + 1
	  vec.error[n.error] <- "Something wrong with the regressor X. It must be a matrix (or list of matrices), a dataframe (or list of dataframes), or an array of order 3."
	}
	
  if(model.obj$kk>0 & length(dim(XX))==3){
	  covariate.temp <- dimnames(XX)[[1]] %in% nomi.covariate
	  XX <- XX[covariate.temp,,]
	  kk <- dim(XX)[1]
	  if(kk!=model.obj$kk){
	    n.error <- n.error+1
	    vec.error[n.error] <- "Some covariates in the model are missing from X."
	  }
	}
	if(kk==1 & length(dim(XX))==2){
	  if(is.null(names(model.obj$beta_coeffs)))
	    names(model.obj$beta_coeffs) <- "varX1"
	}
	
	## checking validity of the spatial matrix
	if(is.null(WW)){
	  n.error <- n.error + 1
	  vec.error[n.error] <- "The spatial matrix W is missing."
	}
	else if(is.list(WW)){
	  if(is.null(WW$ww.index)|is.null(WW$ww.values)){
	    n.error <- n.error + 1
	    vec.error[n.error] <- "The spatial matrix W is missing or not complete."
	  }
	  else{
	    WW.temp <- matrix(0, nrow=pp, ncol=pp)
	    dimnames(WW.temp)[[1]] <- dimnames(WW.temp)[[2]] <- dimnames(dseries)[[1]]
	    for(ii in 1:pp){
	      indici <- as.character(WW$ww.index[ii,WW$ww.index[ii,]>0])
	      WW.temp[ii, indici] <- WW$ww.values[ii,WW$ww.index[ii,]>0]
	    }
	    WW <- WW.temp
	  }
	}
	else if(is.matrix(WW)){
	  if(dim(WW)[1]!=pp){
	    n.error <- n.error + 1
	    vec.error[n.error] <- "The spatial matrix has not the correct number of rows."
	  }
	  if(dim(WW)[2]!=pp){
	    n.error <- n.error + 1
	    vec.error[n.error] <- "The spatial matrix has not the correct number of columns."
	  }
	  dimnames(WW)[[1]] <- dimnames(WW)[[2]] <- dimnames(dseries)[[1]]
	}
	else{
	  n.error = n.error+1
	  vec.error[n.error] <- "Something wrong with the spatial matrix."
	} 
	
	if(sum(is.na(WW))>0){
	  n.error <- n.error + 1
	  vec.error[n.error] <- "The spatial matrix cannot have NA values."
	}
	if(sum(diag(WW))!=0 & var(diag(WW))!=0){
	  n.error <- n.error + 1
	  vec.error[n.error] <- "The spatial matrix has not zero diagonal."
	}
	indici.w <- apply(abs(WW), 1, sum)==0
	if(sum(indici.w)>0){
	  n.error <- n.error + 1
	  vec.error[n.error] <- "The spatial matrix has one (or more than one) zero-row(s)."
	}

	## checking the missing values	
	na <- sum(is.na(dseries))
	if(na>0){
	  n.error = n.error+1
	  vec.error[n.error] <- "The series has NA values. They are replaced with the mean values"
	}
	if(kk==0)
	  naX <- 0
	else if(kk==1)
	  naX <- sum(is.na(XX))
	else if(kk>1)
	  naX <- apply(XX, 1, FUN=function(x){sum(is.na(x))})
	if(sum(naX)!=0){
	  n.warning = n.warning+1
	  vec.warning[n.warning] <- "There are some missing values in the covariates X."
	  na <- c(na, naX)
	  names(na) <- c("naY", paste("naX", seq(1,kk), sep=""))	  
	}

	if(model.obj$time_effects){
	  time_effects <- apply(dseries, 2, mean)
	  dseries <- dseries - matrix(rep(time_effects, pp), nrow=pp, byrow = TRUE)
	}
	else time_effects <- numeric(nn)
	
	if(model.obj$fixed_effects){
	  mu <- apply(dseries, 1, mean)
	}
	else{
	  mu <- numeric(pp)
	  names(mu) <- dimnames(dseries)[[1]]
	} 
	
	## returning the structure of data
	list(series=dseries, WW=WW, XX=XX, neighbours=px.neighbors, mu=mu, time_effects=time_effects, 
	     nn=nn, pp=pp, kk=kk, na=na, model=model.obj, px=px, lat=lat, lon=lon, group=group,
	     index.weights=index.weights, time.weights=time.weights, 
	     errors=vec.error[vec.error!=""], warnings=vec.warning[vec.warning!=""])
}
