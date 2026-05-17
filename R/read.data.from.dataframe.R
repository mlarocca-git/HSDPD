read.data.from.dataframe <- function(px, latit=NULL, longit=NULL,
                                     rry, rrXX=NULL, rrgroups=NULL, model,
                                     ww.index, ww.values, px.neighbors){
  ### estrazione e organizzazione dei dati input di formato matrix o dataframe, con identificazione dei valori mancanti
  ### restituisce la serie spazio-temporale, aggiungendo i vicini dei pixel selezionati, il vettore w_i dei pesi spaziali
  ### ed eventualmente la matrice X dei regressori
  
  if((!is.data.frame(rry)) & (!is.matrix(rry)))
    return(list(error="L'argomento rry deve essere un dataframe oppure una matrice."))
  
  rry <- as.matrix(rry)
  if(is.null(dimnames(rry)[[1]]))
    return(error="I nomi di riga (=locations) di rry devono essere definiti.")
  if(is.null(px)){
    px <- dimnames(rry)[[1]]
  }
  indici <- dimnames(rry)[[1]]
  
  if(is.null(dimnames(rry)[[2]]))
    lat <- lon <-  numeric(0)
  else{
    lat <- which(dimnames(rry)[[2]]=="latitude")
    lon <- which(dimnames(rry)[[2]]=="longitude")
  }
  if(length(lat)!=0 & length(lon)!=0){
    latit <- rry[as.character(px),lat]
    longit <- rry[as.character(px),lon]
    names(latit) <- names(longit) <- px
    rry <- rry[,-c(lat, lon)]
  }
  else if(length(lat)!=0 | length(lon)!=0){
    rry <- rry[,-c(lat, lon)]
  }
  if(is.null(latit) | is.null(longit)){
    latit <- longit <- NULL
  }
  else if(is.numeric(latit) & is.numeric(longit)){
    if(length(latit)==1 & length(longit)==1){
      temp.latit <- rry[as.character(px),latit]
      temp.longit <- rry[as.character(px),longit]
      names(temp.latit) <- names(temp.longit) <- px
      rry <- rry[,-c(latit, longit)]
      latit <- temp.latit
      longit <- temp.longit
    }
    else if(length(latit)==dim(rry)[1] & length(longit)==dim(rry)[1]){
      if(is.null(names(latit)))
        names(latit) <- dimnames(rry)[[1]]
      if(is.null(names(longit)))
        names(longit) <- dimnames(rry)[[1]]
      latit <- latit[as.character(px)]
      longit <- longit[as.character(px)]
    } 
    else return(list(error="I valori di latitudine e longitudine non hanno il formato o valore atteso"))
  }
  else return(list(error="I valori di latitudine e longitudine non hanno il formato o valore atteso"))
  
  if(sum(as.character(px) %in% indici) < length(px))  
    return(list(error="Alcuni valori di px non sono contenuti nel dataframe rry"))
  
  ### derivazione serie dei vicini-stretti
  if(is.null(ww.index) | is.null(ww.values))
    return(list(error="La matrice spaziale non può essere costruita senza ww.index e ww.values..."))
  if(is.matrix(ww.index)){
    vicini.stretti <- as.numeric(ww.index[as.character(px),])
    vicini.stretti <- vicini.stretti[vicini.stretti>0]
    vicini.stretti <-  unique(vicini.stretti)
    indici <- as.character(unique(c(px, vicini.stretti)))
    if(length(indici)>1)
      serie <- as.matrix(rry[indici,])
    else if(length(indici)==1){
      serie <- matrix(rry[indici,], nrow=1)
      dimnames(serie)[[1]] <- indici      
    }
    else return(list(error="Non vi sono dati nella matrice series"))
    tempo <- dimnames(serie)[[2]]
    tt <- dim(serie)[2]
    pp <- length(indici)
    n.px <- length(px)
  }
  else  return(list(error="L'argomento ww.index deve essere una matrice"))
  
  ### derivazione dei gruppi
  if(is.null(rrgroups) & is.null(model$groups)){
    gruppi <- rep(1, length(px))
    labels <- rep("group_1", length(px))
    gruppi <- cbind(COD=gruppi, LABEL=labels)
    dimnames(gruppi)[[1]] <- px
  }
  else if(!is.null(rrgroups)){
    if((!is.data.frame(rrgroups)) & (!is.matrix(rrgroups)))
      return(list(error="L'argomento rrgroups deve essere un dataframe o matrice"))
    if(dim(rrgroups)[1] != dim(rry)[1])
      return(list(error="L'oggetto rrgroups deve avere la stessa dimensione (=numero righe) di rry"))
    if(sum(c("COD", "LABEL") %in% dimnames(rrgroups)[[2]])<2)
      return(list(error="L'oggetto rrgroups deve contenere le colonne COD e LABEL"))
    gruppi <- as.matrix(rrgroups[as.character(px),c("COD", "LABEL")])
  }
  else if(!is.null(model$groups))
    gruppi <- as.matrix(model$groups[indici,])
  
  ### derivazione matrici dei pesi spaziali e dei punti vicini-lontani
  if(is.matrix(ww.values)){
    ww.values <- ww.values[indici,]
    ww.index <- ww.index[indici,]
  }  
  else  return(list(error="L'argomento ww.values deve essere una matrice."))
  
  ### correzione indici e pesi della matrice spaziale e della matrice dei vicini, introducendo effetti boundary
  if(is.null(px.neighbors))
    return(error="La matrice px.neighbours non può mancare.")
  indici.intorno <- unique(na.omit(as.vector(px.neighbors$index[as.character(px),])))
  neighbors <- as.matrix(px.neighbors$index[indici,])
  rimanenti <- indici[!(indici %in% px)]
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
  varX <- names(model$beta_coeffs)[model$beta_coeffs]
  if(kk>0 & is.list(rrXX)){
    XX <- array(0, dim=c(kk, pp, tt))
    dimnames(XX) <- list(varX, indici, tempo) 
    for(rr in varX){
      if(rr=="trend"){
        XX[rr,,] <- matrix(rep(seq(1, tt)/tt, pp), byrow=TRUE, nrow=pp)
      }
      else if(rr %in% names(rrXX)){
        if(is.null(rrXX[[rr]]))
          return(list(error="C'è un oggetto NULL nella lista di covariate rrXX"))
        XX[rr,,] <- as.matrix(rrXX[[rr]][indici,as.character(tempo)])
      }
      else
        return(list(error="L'argomento rrXX non contiene tutte le covariate previste dal modello SDPD"))
    }
  }
  else if(kk>0 & is.array(rrXX)){
    XX <- rrXX[varX, indici, ]
  }
  
  ### output
  res <- list(series=serie, X=XX, ww.index=ww.index[indici,], ww.values=ww.values[indici,], px.neighbors=list(index=neighbors, seriesBoundary=seriesBoundary), 
              px=px, lon=longit, lat=latit, group=gruppi)
  res
}

