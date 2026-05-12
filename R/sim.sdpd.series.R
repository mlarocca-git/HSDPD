sim.sdpd.series <- function(series=NULL, model, markovian=FALSE, only_inside=FALSE, num.steps=5, sigma.eps=1) 
{
  ## dseries is a matrix of dimension (nn,pp) where pp is the number of univariate time dseries and data$nn the number of time observations
  ## W is a spatial weight matrix of dimension (pp,pp)
  ## X is an array of dim=c(kk,nn,pp) which includes the data for kk exogenous regressors. If kk=1 then X is a matrix of dim=(nn,pp)
  
  if(is.null(series)){
    cat("\n Errore: df.data is NULL")
    return(NULL)
  }
  
  
  ## checking validity of data...
  data <- check.sdpd.series(series=series, model=model)
  
  if(length(data$errors)>0){
    cat("\n There are some errors with the series....\n", data$errors)
    return(list(errors=data$errors, warnings=data$warnings))
  }
  if(is.null(series$coeffs)){
    cat("\n I dati non possono essere simulati senza la matrice dei coefficienti, da passare in model$coeffs")
    return(list(error="I dati non possono essere simulati senza la matrice dei coefficienti, da passare in model$coeffs"))
  }
  if(is.null(series$eps)){
    series$eps <- matrix(rnorm(data$nn*data$pp, sd=sigma.eps), nrow=data$pp, ncol=data$nn)
    dimnames(series$eps) <- dimnames(series$series)
  }
  else if(dim(series$eps)[1]!=data$pp | dim(series$eps)[2]!=data$nn)
    return(list(error="La dimensione della matrice di errori passata in model$eps non è compatibile con la matrice dati passata"))
  
  
  ## simulating values...
  px.sim <- NULL
  if(only_inside)
    px.sim <- as.character(series$px)
  
  
  series$model <- model
  series$series <- fit.sdpd.series(dseries=data$series, W=data$WW, X.centr=data$XX, px.sim=px.sim, num.steps=num.steps, 
                        model=data$model, markovian=markovian, coeff.hat=series$coeffs, resid=series$eps)$series

  series
}
