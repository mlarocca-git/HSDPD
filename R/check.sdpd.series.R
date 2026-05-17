check.sdpd.series <- function(series, XX=NULL, model=NULL, ww.index=NULL, ww.values=NULL,
                              px.neighbors=NULL, px=NULL, lat=NULL, lon=NULL,  group=NULL,
                              index.weights=NULL, time.weights=NULL) 
{
  ## Questa funzione riceve un oggetto sdpd-series e un oaggetto sdpd-model (costruiti mediante le funzioni
  ## build.sdpd.series e build.sdpd.model), oppure riceve le singole componenti di sdpd-series.
  ## La funzione, quindi, controlla la consistenza dei componenti della serie con il modello e
  ## restituisce gli oggetti "controllati" insieme alla lista di errori e warnings
  
	vec.error <- vec.warning <- character(20)
	n.error <- n.warning <- 0
  model.obj <- model
  ww.index.temp <- ww.values.temp <- NULL

	## checking validity of series
	if(is.list(series)){
	  ## in questo caso tutti gli oggetti vengono passati attraverso una lista (procedura parallelizzata)
	  dseries <- series$series
	  px <- series$px
	  lon <- series$lon
	  lat <- series$lat
	  group <- series$group
	  px.neighbors <- series$px.neighbors
	  ww.index.temp <- series$ww.index
	  ww.values.temp <- series$ww.values
	  if(is.null(XX))
	    XX <- series$X
	}
	else{
	  dseries <- series
	}
	
	if(is.matrix(dseries)|is.data.frame(dseries)){
	  dseries <- as.matrix(dseries)
	  nn <- dim(dseries)[2]
		pp <- dim(dseries)[1]
		if(is.null(dimnames(dseries)[[1]]) | is.null(dimnames(dseries)[[2]])){
		  n.error <- n.error + 1
		  vec.error[n.error] <- "\n Missing names for dseries (they should be numbers, for compatibility with plot functions)"
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
		    if(is.null(dimnames(index.weights)[[1]]))
		      dimnames(index.weights)[[1]] <- dimnames(dseries)[[1]]
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
	if(is.null(model.obj)){
	  n.error <- n.error + 1
	  vec.error[n.error] <- "The model is missing. Please check."
	  nomi.covariate <- NULL
	}
	else{
	  model.obj$kk <- sum(model.obj$beta_coeffs)
	  nomi.covariate <- names(model$beta_coeffs)[model$beta_coeffs]
	}
	if(is.null(XX)){
	  if(!is.null(model.obj$kk) & model.obj$kk>0){
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
	  if(is.null(dimnames(XX)[[1]]) | is.null(dimnames(XX)[[2]])){
	    n.error <- n.error + 1
	    vec.error[n.error] <- "Dimnames for XX are missing."
	  }
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
	  if(is.null(dimnames(XX)[[1]]) | is.null(dimnames(XX)[[2]]) | is.null(dimnames(XX)[[3]])){
	    n.error <- n.error + 1
	    vec.error[n.error] <- paste("Dimnames for X are missing.")
	  }
	  kk <- dim(XX)[1]
	}
	else{
	  n.error <- n.error + 1
	  vec.error[n.error] <- "Something wrong with regressor X. It must be a matrix (or list of matrices), a dataframe (or list of dataframes), or an array of order 3."
	  kk <- 0
	}
	
	if(kk>0 & length(dim(XX))==3){
	  covariate.temp <- dimnames(XX)[[1]] %in% nomi.covariate
	  XX <- XX[covariate.temp,,]
	  kk <- dim(XX)[1]
	  if(kk!=model.obj$kk){
	    n.error <- n.error+1
	    vec.error[n.error] <- "Some covariates in the model are missing from X."
	  }
	}
	if(kk==1 & length(dim(XX))==2){
	  if(is.null(names(model.obj$beta_coeffs))){
	    n.error <- n.error+1
	    vec.error[n.error] <- "The name of the covariate in the model is missing."
	  }
	}
	
	## checking validity of the spatial matrix
	if(is.null(ww.index))
	  ww.index <- ww.index.temp
	if(is.null(ww.values))
	  ww.values <- ww.values.temp
	if(!is.null(model) & (is.null(ww.index) | is.null(ww.values))){
	  ww.index <- model$ww.index
	  ww.values <- model$ww.values
	}
	if(is.null(ww.index)|is.null(ww.values)){
	  n.error <- n.error + 1
	  vec.error[n.error] <- "The spatial matrix components are missing."
	}
	else if(is.null(dimnames(ww.index)[[1]]) | is.null(dimnames(ww.values)[[1]])){
	  n.error = n.error+1
	  vec.error[n.error] <- "Names for spatial matrix are missing"
	}

	## checking the missing values	
	na <- sum(is.na(dseries))
	if(na>0){
	  n.error = n.error+1
	  vec.error[n.error] <- "The series has NA values."
	}
	if(kk==0)
	  naX <- 0
	else if(kk==1)
	  naX <- sum(is.na(XX))
	else if(kk>1)
	  naX <- apply(XX, 1, FUN=function(x){sum(is.na(x))})
	if(sum(naX)!=0){
	  n.error = n.error+1
	  vec.error[n.error] <- "There are some missing values in the covariates X."
	  na <- c(na, naX)
	  names(na) <- c("naY", paste("naX", seq(1,kk), sep=""))	  
	}

	if(model.obj$time_effects){
	  time_effects <- apply(dseries, 2, mean)
	  dseries <- dseries - matrix(rep(time_effects, pp), nrow=pp, byrow = TRUE)
	}
	else time_effects <- numeric(nn)
	
	## checking lon, lat and group
	if(!is.null(lon)){
	  if(!is.vector(lon) | !is.numeric(lon) | is.null(names(lon))){
	    n.error = n.error+1
	    vec.error[n.error] <- "Something wrong with the vector of longitudes (missing names?)"
	  }
	}
	if(!is.null(lat)){
	  if(!is.vector(lat) | !is.numeric(lat) | is.null(names(lat))){
	    n.error = n.error+1
	    vec.error[n.error] <- "Something wrong with the vector of latitudes (missing names?)."
	  }
	}
	if(!is.null(group)){
	  if(length(dim(group))!=2 | dim(group)[2]!=2 | is.null(dimnames(group)[[1]])){
	    n.error = n.error+1
	    vec.error[n.error] <- "Something wrong with the object group (missing names?)."
	  }
	  else if(sum(c("COD", "LABEL") %in% dimnames(group)[[2]])<2){
	    n.error = n.error+1
	    vec.error[n.error] <- "L'oggetto groups deve contenere le colonne COD e LABEL"
	  }
	}
	
	if(model.obj$fixed_effects){
	  mu <- apply(dseries, 1, mean)
	}
	else{
	  mu <- numeric(length(px))
	  names(mu) <- dimnames(dseries)[[1]]
	} 
	
	if(is.null(px.neighbors)){
	  n.error <- n.error + 1
	  vec.error[n.error] <- "\n The object px.neighbors is missing."
	}
	

	## returning the structure of data
	list(series=dseries, ww.index=ww.index, ww.values=ww.values, XX=XX, px.neighbors=px.neighbors, time_effects=time_effects, 
	     nn=nn, pp=pp, kk=kk, na=na, px=px, lat=lat, lon=lon, group=group,
	     index.weights=index.weights, time.weights=time.weights, mu=mu, 
	     errors=vec.error[vec.error!=""], warnings=vec.warning[vec.warning!=""])
}
