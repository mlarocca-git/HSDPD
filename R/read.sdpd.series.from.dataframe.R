read.sdpd.series.from.dataframe <- function(px=NULL, latit=NULL, longit=NULL,
                                          rry, rrXX=NULL, rrgroups=NULL, model,
                                          ww.index, ww.values, px.neighbors){
  ### estrazione e organizzazione dei dati input di formato matrix o dataframe, con identificazione dei valori mancanti
  ### restituisce la serie spazio-temporale, aggiungendo i vicini dei pixel selezionati, il vettore w_i dei pesi spaziali
  ### ed eventualmente la matrice X dei regressori
  
  if((!is.data.frame(rry)) & (!is.matrix(rry)))
    return(list(error="L'argomento rry deve essere un dataframe oppure una matrice"))

  rry <- as.matrix(rry)
  if(is.null(dimnames(rry)[[1]]))
    dimnames(rry)[[1]] <- seq(1, dim(rry)[1])
  if(is.null(px)){
    px <- dimnames(rry)[[1]]
  }
  indici <- dimnames(rry)[[1]]

  lat <- which(dimnames(rry)[[2]]=="latitude")
  lon <- which(dimnames(rry)[[2]]=="longitude")
  temp.latit <- temp.longit <- NULL
  if(length(lat)!=0 & length(lon)!=0){
    temp.latit <- rry[as.character(px),lat]
    temp.longit <- rry[as.character(px),lon]
    names(temp.latit) <- names(temp.longit) <- px
    rry <- rry[,-c(lat, lon)]
  }
  else if(length(lat)!=0 | length(lon)!=0){
    rry <- rry[,-c(lat, lon)]
  }
  if(is.null(latit) | is.null(longit)){
    latit <- temp.latit
    longit <- temp.longit
  }
  else if(is.numeric(latit) & is.numeric(longit)){
    if(length(latit)==1 & length(longit)==1){
      temp.latit <- rry[as.character(px),latit]
      temp.longit <- rry[as.character(px),longit]
      rry <- rry[,-c(latit, longit)]
      latit <- temp.latit
      longit <- temp.longit
    }
    else{
      latit <- latit[as.character(px)]
      longit <- longit[as.character(px)]
    } 
  }
  else return(list(error="I valori di latitudine e longitudine non hanno il formato o valore atteso"))

  if(sum(as.character(px) %in% indici) < length(px))  
    return(list(error="Alcuni valori di px non sono contenuti nel dataframe rry"))
  
  ### derivazione serie dei vicini-stretti
  if(is.null(ww.index) | is.null(ww.values))
    return(list(error="La matrice spaziale non può essere costruita senza ww.index e ww.values..."))
  if(is.data.frame(ww.index) | is.matrix(ww.index)){
    vicini.stretti <- as.numeric(ww.index[as.character(px),])
    vicini.stretti <- ifelse(vicini.stretti==0, NA, vicini.stretti)
    vicini.stretti <-  na.exclude(unique(vicini.stretti))
    indici <- as.character(unique(c(px, vicini.stretti)))
    if(length(indici)>1)
      serie <- as.matrix(rry[indici,])
    else if(length(indici)==1){
      serie <- matrix(rry[indici,], nrow=1)
      dimnames(serie)[[1]] <- indici      
    }
    else return(list(error="Non vi sono dati nella matrice series"))
    tempo <- dimnames(serie)[[2]]
    tt <- length(tempo)
    pp <- length(indici)
    n.px <- length(px)
  }
  else  return(list(error="L'argomento ww.index deve essere un dataframe o matrice"))

  ### derivazione dei gruppi
  if(is.null(rrgroups)){
    gruppi <- rep(1, nrow=pp)
    labels <- rep("Common group", nrow=pp)
    gruppi <- cbind(COD=gruppi, LABEL=labels)
  }
  else{
    if((!is.data.frame(rrgroups)) & (!is.matrix(rrgroups)))
      return(list(error="L'argomento rrgroups deve essere un dataframe o matrice"))
    if(dim(rrgroups)[1] != dim(rry)[1])
      return(list(error="L'oggetto rrgroups deve avere la stessa dimensione (=numero righe) di rry"))
    if(sum(c("COD", "LABEL") %in% dimnames(rrgroups)[[2]])<2)
      return(list(error="L'oggetto rrgroups deve contenere le colonne COD e LABEL"))
    gruppi <- as.matrix(rrgroups[as.character(px),c("COD", "LABEL")])
  }
  
  ### derivazione matrici dei pesi spaziali e dei punti vicini-lontani
  if(is.data.frame(ww.values) | is.matrix(ww.values)){
    ww.values <- as.matrix(ww.values[indici,])
    ww.index <- as.matrix(ww.index[indici,])
  }  
  else  return(list(error="L'argomento ww.index deve essere un dataframe o matrice"))
    
  ### correzione indici e pesi della matrice spaziale e della matrice dei vicini, introducendo effetti boundary
  if(is.null(px.neighbors$index))
    return(list(error="La matrice spaziale non può essere costruita senza px.neighbors$index..."))
  indici.intorno <- unique(na.omit(as.vector(px.neighbors$index[as.character(px),])))
  neighbors <- as.matrix(px.neighbors$index[indici,])
  rimanenti <- indici[!(as.numeric(indici) %in% px)]
  for(riga in rimanenti){
    presenti <- ww.index[riga,]%in%px
    ww.index[riga, !presenti] <- 0
    ww.values[riga, !presenti] <- 0
    ww.values[riga,] <- ww.values[riga,]/sum(abs(ww.values[riga,]))
    presenti <- neighbors[riga,]%in%as.numeric(indici.intorno)
    neighbors[riga, !presenti] <- NA
  }
  
  ### derivazione serie di intorno dei vicini-lontani
  indici.intorno <- as.character(unique(na.omit(as.vector(neighbors))))
  indici.esterni <- indici.intorno[!(indici.intorno %in% indici)]
  indici.residui <- indici.esterni[!(indici.esterni %in% dimnames(rry)[[1]])]
  indici.esterni <- indici.esterni[indici.esterni %in% dimnames(rry)[[1]]]
  if(length(indici.esterni)>1)
    seriesBoundary_first <- as.matrix(rry[as.character(indici.esterni),])
  else if(length(indici.esterni)==1){
    seriesBoundary_first <- matrix(rry[as.character(indici.esterni),], nrow=1)
    dimnames(seriesBoundary_first)[[1]] <- indici.esterni
  }
  else
    seriesBoundary_first <- NULL
  
  
  seriesBoundary_second <- NULL
  if(length(indici.residui)>0 & is.null(px.neighbors$seriesBoundary)){
    temp <- dimnames(neighbors)
    neighbors <- t(apply(neighbors, 1, FUN=function(x, ind){ifelse(x%in%ind, NA, x)}, ind=as.numeric(indici.residui)))
    dimnames(neighbors) <- temp
  }
  else if(length(indici.residui)>1)
    seriesBoundary_second <- as.matrix(px.neighbors$seriesBoundary[as.character(indici.residui),])
  else if(length(indici.residui)==1){
    seriesBoundary_second <- matrix(px.neighbors$seriesBoundary[as.character(indici.residui),], nrow=1)
    dimnames(seriesBoundary_second)[[1]] <- indici.residui
  }
  seriesBoundary <- rbind(seriesBoundary_first, seriesBoundary_second)

  ### creazione regressori
  if(is.null(rrXX) | sum(model$beta_coeffs)==0){
    kk <- 0
    XX <- varX <- NULL
    if(sum(model$beta_coeffs)>0)
      return(list(error="Il modello prevede dei regressori esogeni, tuttavia l'argomento rrXX non è stato valorizzato"))
  }
  else kk <- sum(model$beta_coeffs)

  if(kk>0 & is.list(rrXX)){
    varX <- names(model$beta_coeffs)[model$beta_coeffs]
    XX <- array(0, dim=c(kk, pp, tt))
    dimnames(XX) <- list(varX, indici, tempo) 
    for(rr in varX){
      if(rr=="trend"){
        XX[rr,,] <- matrix(rep(seq(1, tt)/tt, pp), byrow=TRUE, nrow=pp)
      }
      else if(rr %in% names(rrXX)){
        if(is.null(rrXX[[rr]]))
          return(list(error="C'è un valore NULL nella lista di covariate rrXX"))
        XX[rr,,] <- as.matrix(rrXX[[rr]][indici,as.character(tempo)])
      }
      else
        return(list(error="L'argomento rrXX non contiene tutte le covariate previste dal modello SDPD"))
    }
  }
  else if(kk>0 & is.array(rrXX)){
    varX <- names(model$beta_coeffs)[model$beta_coeffs]
    XX <- rrXX[varX, indici, tempo] 
  }
  
  ### output
  res <- list(series=serie, X=XX, ww.index=ww.index, ww.values=ww.values, px.neighbors=list(index=neighbors, seriesBoundary=seriesBoundary), 
              px=px, lon=longit, lat=latit, group=gruppi)
  
  ## in case of simulations, we also return the true coefficients and the error series
  coeffs <- model$coeffs
  if(!is.null(coeffs)){
    res$coeffs <- coeffs[indici,]
  }
  series.eps <- model$eps
  if(!is.null(series.eps)){
    res$eps <- model$eps[indici,]
  }
  
  res
}
