build.sdpd.series <- function(df.obj=NULL, px=NULL, rry=NULL, rrxx=NULL, lon=NULL, lat=NULL,
                              rrgroups=NULL, label_groups=NULL,
                              model, check=TRUE,
                              ww.index=NULL, ww.values=NULL,
                              px.neighbors=NULL,
                              SIM=FALSE, nn=NULL,
                              vec.options=list(px.core=1, px.neighbors=6, na.rm=T, NAcovs="pairwise.complete.obs",
                                               covariates.sim.model=list(ar = c(0.8), sd = 1),
                                               markovian=TRUE, num.steps=10)){
## Questa funzione costruisce l'oggetto sdpd-series, così come atteso dalle successive funzioni di stima.
## I dati vengono simulati quando non sono forniti in input, a condizione che SIM=TRUE.
## Si noti che l'unico argomento non predefinito è l'oggetto sdpd-model, costruito con la funzione build.sdpd.model()
## Importante: il primo oggetto df.obj serve per permettere l'utilizzo di questa funzione in modo parallelizzato
## Si noti che, se SIM=TRUE, la serie rry eventualmente passata viene utilizzata soltanto per inizializzare la serie simulata
  
  if(is.null(px))
    px <- df.obj$px

  if(!SIM & is.null(rry))
    return(list(error="L'oggetto rry non può mancare, a meno che non si imposti il parametro SIM=TRUE"))
  
  if(SIM){
    new.data <- generate.sdpd.series(nn=nn, rry=rry, rrxx=rrxx, model=model, markovian=vec.options$markovian, 
                                num.steps=vec.options$num.steps, covariates.sim.model=vec.options$covariates.sim.model) 
    rry <- new.data$rry
    rrxx <- new.data$rrxx
    ww.index <- model$ww.index
    ww.values <- model$ww.values
    px.neighbors <- model$px.neighbors
    rrgroups <- model$groups
  }
  if(is.data.frame(rry) | is.matrix(rry))
    res <- read.data.from.dataframe(px=px, rry=rry, rrXX=rrxx, model=model, longit=lon, latit=lat,
                rrgroups=rrgroups, ww.index=ww.index, ww.values=ww.values, px.neighbors=px.neighbors)
  else if(inherits(rry, "SpatRaster")){
    res <- read.data.from.raster(px=df.obj$px, rry=rry, rrXX=rrxx, lon=lon, lat=lat,
                rrgroups=rrgroups, model=model, label_groups=label_groups, vec.options=vec.options)
  }
  else return(list(error="Gli oggetti rry e rrxx non sono del formato atteso (dataframe, matrix oppure raster)"))
  
  
  ## diagnostics
  # res$model <- model
  if(check){
    check <- check.sdpd.series(series=res, model=model)
    if(length(check$errors)>0){
      return(list(errors=check$errors, warnings=check$warnings))
    }
    res$warnings <- check$warnings
  }
  
  ## output  
  res
}

