read.data.from.raster <- function(px=NULL, latit=NULL, longit=NULL, 
                                  rry, rrXX=NULL, 
                                  rrgroups=NULL, label_groups=NULL,
                                  model,
                                  vec.options,
                                  type.w="distance"){
  # Estrazione e costruzione serie spazio-temporale a partire da variabili in formato raster (mediante package terra)
  # Attenzione: l'argomento px deve essere in prima posizione, per poter parallelizzare la procedura
  # Se px="all"
  #   - costruisce la serie spazio-temporale completa, controllando tutti i pixel presenti ed eliminando quelli con valori missing
  # Se px=vettore (oppure se sono passati i vettori latit e longit)
  #   - costruisce la serie spazio-temporale per i soli pixel selezionati, aggiungendo il "bordo-vicini" dei pixel selezionati; 
  # Restituisce: i dati della endogena, il vettore w_i dei pesi spaziali e la matrice X dei regressori
  
  if(is.null(px)){
    if(is.null(latit)|is.null(longit)){
      px <- "all"
    }
    else{
      if(length(latit)==length(longit))
        n.px <- length(latit)
      else
        return(list(error="Latitudine e longitudine sono di lunghezza diversa."))
      px <- cellFromXY(rry, cbind(longit, latit))
    }
  }
  if(is.numeric(px)){
    n.px <- length(px)
    temp <- xyFromCell(rry, px)
    latit <- temp[,"y"]
    longit <- temp[,"x"]
    names(latit) <- names(longit) <- px
  }
  else if(is.character(px)){
    n.px <- length(px)
    if(px[1]=="all"){
      px <- seq(1, dim(values(rry))[1])
    }
    else
      px <- as.numeric(px)
    temp <- xyFromCell(rry, px)
    latit <- temp[,"y"]
    longit <- temp[,"x"]
    names(latit) <- names(longit) <- px
  }
  if(n.px==0)
    return(list(error="Non ci sono serie da analizzare."))
  ### inizializzazione serie includenti i punti dei vicini stretti (=mat.core, che considera px.core cerchi intorno ad ogni pixel)
  if(vec.options$px.core==0)
    vicini.stretti <- ww.index <- ww.values <- px.neighbors <- NULL
  else{
    mat.core <- matrix(1, ncol=2*vec.options$px.core+1, nrow=2*vec.options$px.core+1)
    mat.core[vec.options$px.core+1,vec.options$px.core+1]<- 0
    vicini.stretti <- adjacent(rry, px, pairs=TRUE, directions=mat.core, include=TRUE, symmetrical=FALSE)    
  }
  indici <- na.exclude(unique(c(px, vicini.stretti)))
  tempo <- time(rry)
  serie <- values(rry)[indici,]
  dimnames(serie)[[1]] <- indici
  dimnames(serie)[[2]] <- as.character(tempo)
  tt <- length(tempo)
  pp <- dim(serie)[1]
  
  ### creazione regressori
  coordinate <- xyFromCell(rry, indici)
  if(is.null(rrXX)){
    kk <- n.reg <- 0
    XX <- varX <- NULL
    if(sum(model$beta_coeffs)>0)
      return(list(error="Il modello prevede dei regressori esogeni, tuttavia l'argomento rrXX non è stato valorizzato"))
  }
  else if(!is.list(rrXX)){
    kk <- n.reg <- 1
    if(sum(model$beta_coeffs)>1)
      return(list(error="Il modello prevede più di un regressore esogeno, tuttavia l'argomento rrXX non contiene variabili sufficienti"))
  }
  else if(is.list(rrXX)){
    kk <- length(rrXX)
    if(sum(model$beta_coeffs)>kk)
      return(list(error=paste("Il modello prevede più di", kk, "regressori esogeni, tuttavia l'argomento rrXX non contiene variabili sufficienti")))
    if(!is.null(names(model$beta_coeffs))){
      if(sum(names(model$beta_coeffs)[model$beta_coeffs]%in%names(rrXX))<sum(model$beta_coeffs))
        return(list(error="Il modello prevede alcuni regressori esogeni che non sono contenuti in rrXX"))
    }
  }
  if(kk>0){
    n.reg <- 0
    varX <- character(kk)
    XX <- array(0, dim=c(kk, pp, tt))
    for(rr in 1:dim(XX)[1]){
      if(!is.null(rrXX[[rr]])){
        n.reg <- n.reg+1
        if(is.character(rrXX[[rr]])){
          if(rrXX[[rr]]=="trend"){
            XX[n.reg,,] <- matrix(rep(seq(1, tt)/tt, pp), byrow=TRUE, nrow=pp)
            varX[n.reg] <- "trend"
          }
          else n.reg <- n.reg-1
        }
        else if(!is.null(names(model$beta_coeffs))){
          if(names(rrXX)[rr] %in% names(model$beta_coeffs)){
            varX[n.reg] <- names(rrXX)[rr]
            pXX <- cellFromXY(rrXX[[rr]], coordinate)
            XX[n.reg,,] <- values(rrXX[[rr]])[pXX,]
          }
          else n.reg <- n.reg-1
        }else{
          varX[n.reg] <- names(rrXX)[rr]
          pXX <- cellFromXY(rrXX[[rr]], coordinate)
          XX[n.reg,,] <- values(rrXX[[rr]])[pXX,]
        }
      }
    }
    if(n.reg>0){
      kk <- n.reg
      dimnames(XX) <- list(varX, indici, as.character(tempo)) 
    }
  }
  ### definizione dei gruppi (in caso di alta dimensione, per frazionare la stima su sottogruppi)
  if(is.null(rrgroups)){
    gruppi <- rep(1, nrow=pp)
    labels <- rep("Common group", nrow=pp)
    gruppi <- cbind(COD=gruppi, LABEL=labels)
  }
  else{
    gruppi <- values(rrgroups)[indici,1]
    labels <- label_groups[as.character(gruppi)]
    gruppi <- data.frame(COD=gruppi, LABEL=labels)
  }
  dimnames(gruppi)[[1]] <- indici
  
  ## eliminazione dei pixel con valori NA
  n.NAY <- apply(serie, 1, FUN=function(x){sum(is.na(x))})
  n.NAg <- is.na(gruppi[,1])
  n.NAgor <- xor(n.NAY>0,n.NAg>0)
  if(kk>1){
    n.NAX <- t(apply(XX, c(1,2), FUN=function(x){sum(is.na(x))}))
    n.NAxor <- apply(n.NAX, 2, FUN=function(x,y){xor(x>0,y>0)}, y=n.NAY)
    tot.NAX <- apply(n.NAX, 1, sum)
  }
  else if(kk==1){
    n.NAX <- tot.NAX <- apply(XX, 1, FUN=function(x){sum(is.na(x))})
    n.NAxor <- xor(n.NAX>0,n.NAY>0)
  }
  else n.NAX <- n.NAxor <- tot.NAX <- rep(0, pp)
  if(vec.options$na.rm)
    da.mantenere <- n.NAY==0 & n.NAg==0 & tot.NAX==0
  else
    da.mantenere <- n.NAY<tt & n.NAg==0  & tot.NAX<tt
  indici.na <- indici[!da.mantenere]
  if(sum(da.mantenere)==0)
    return(list(error="Non ci sono serie da analizzare."))
  varY <- varnames(rry)[1]
  if(kk>0){
    n.NA <- cbind(coordinate, n.NAY, n.NAX, n.NAxor, n.NAgor)
    dimnames(n.NA) <- list(indici, c("lat", "lon", varY, varX, paste(varY, "XOR", varX, sep=""), paste(varY, "XORGroups")))
  }
  else if(kk==0){
    n.NA <- cbind(coordinate, n.NAY, n.NAgor)
    dimnames(n.NA) <- list(indici, c("lat", "lon", varY, paste(varY, "XORGroups")))
  }
  
  ### creazione dei vettori dei pesi spaziali w_i e verifica dei punti isolati
  if(vec.options$px.core>0){
    coordinate1 <- xyFromCell(rry, vicini.stretti[,1])
    coordinate2 <- xyFromCell(rry, vicini.stretti[,2])
    ww.index <- ww.values <- matrix(0, nrow=sum(da.mantenere), ncol=(2*vec.options$px.core+1)^2)
    dimnames(ww.index)[[1]] <- dimnames(ww.values)[[1]] <- indici[da.mantenere]
    if(type.w=="distance"){
      WW <- distance(x=coordinate1, y=coordinate2, lonlat=TRUE, pairwise=TRUE)
      WW <- ifelse(WW<0.01, 0, 1/WW)
    }
    else return(list(error="The values set for type.w is not allowed."))
    punti.isolati <- logical(pp)
    names(punti.isolati) <- indici
    for(ii in 1:sum(da.mantenere)){
      pixel <- indici[da.mantenere][ii]
      # controlliamo prima nella prima colonna di vicini.stretti
      ww1 <- vicini.stretti[,1]==pixel
      ww2 <- vicini.stretti[ww1,2]
      to.exclude <- ww2 %in% indici.na
      ww2 <- ww2[!to.exclude]
      if(sum(ww1)-sum(to.exclude)>0){
        ww.index[ii,1:length(ww2)] <- ww2
        ww.values[ii,1:length(ww2)] <- WW[ww1][!to.exclude]
      }
      else{
        # per completezza, controlliamo anche eventuali punti aggiuntivi nella seconda colonna di vicini.stretti
        ww1 <- vicini.stretti[,2]==pixel
        ww2 <- vicini.stretti[ww1,1]
        to.exclude <- ww2 %in% indici.na
        ww2 <- ww2[!to.exclude]
        if(sum(ww1)-sum(to.exclude)>0){
          ww.index[ii,1:length(ww2)] <- ww2
          ww.values[ii,1:length(ww2)] <- WW[ww1][!to.exclude]
        }
      }
      if(sum(abs(ww.values[ii,]))>0)
        ww.values[ii,] <- ww.values[ii,]/sum(abs(ww.values[ii,]))
      else
        punti.isolati[da.mantenere][ii] <- TRUE
    }
    ww.values <- ww.values[!(punti.isolati[da.mantenere]),]
    ww.index <- ww.index[!(punti.isolati[da.mantenere]),]
    da.mantenere <- da.mantenere & !punti.isolati
    n.NA <- cbind(n.NA, punti.isolati=punti.isolati)
  }
  
  ## ridefinizione delle quantità al netto dei punti isolati e missing values
  indici <- as.character(indici[da.mantenere])
  serie <- serie[indici,]
  XX <- XX[,indici,]
  i.px <- as.character(px) %in% indici
  latit <- latit[i.px]
  longit <- longit[i.px]
  px <- px[i.px]
  n.px <- length(px)
  if(n.px==0)
    return(list(error="Tutte le serie scelte hanno valori NA nella Y o nelle covariate X o nei dintorni."))
  if(!is.null(gruppi))
    gruppi <- gruppi[as.character(px),]
  
  ### creazione matrice index con gli indici delle serie ricadenti nell'intorno dei vicini-lontani
  if(vec.options$px.neighbors>0){
    mat.intorno <- matrix(1, ncol=2*vec.options$px.neighbors+1, nrow=2*vec.options$px.neighbors+1)
    mat.intorno[vec.options$px.neighbors+1,vec.options$px.neighbors+1]<- 0
    px.neighbors <- matrix(NA, nrow=length(indici), ncol=(2*vec.options$px.neighbors+1)^2-1)
    dimnames(px.neighbors)[[1]] <- indici
    for(ii in 1:n.px){
      temp <- adjacent(rry, cells=px[ii], directions=mat.intorno)
      temp <- temp[temp>0]
      px.neighbors[as.character(px[ii]), 1:length(temp)] <- temp
    }
    indici.intorno <- na.exclude(unique(as.vector(px.neighbors)))
    rimanenti <- indici[!(as.numeric(indici) %in% px)]
    for(ii in rimanenti){
      temp <- adjacent(rry, cells=as.numeric(ii), directions=mat.intorno)
      temp <- temp[temp>0]
      temp <- temp[temp %in% indici.intorno]
      px.neighbors[as.character(ii), 1:length(temp)] <- temp
    }
    ## controllo se alcune serie dell'intorno dei vicini-lontani includono NA, in tal caso le elimino dall'intorno
    serie.intorno <- values(rry)[indici.intorno,]
    dimnames(serie.intorno)[[1]] <- indici.intorno
    dimnames(serie.intorno)[[2]] <-  as.character(tempo)
    na.intorno <- apply(serie.intorno, 1, FUN=function(x){sum(is.na(x))})
    px.neighbors <- t(apply(px.neighbors, 1, FUN=function(x, ind){ifelse(x%in%ind, NA, x)}, ind=indici.intorno[na.intorno>0]))
    ## infine, estraggo la serie dei vicini-lontani, che sarà utilizzata per il calcolo delle covarianze
    indici.intorno <- na.exclude(unique(as.vector(px.neighbors)))
    esterni <- indici.intorno[!(indici.intorno %in% indici)]
    serie.intorno <- serie.intorno[as.character(esterni),]
    px.neighbors <- list(index=px.neighbors, seriesBoundary=serie.intorno)
  }
  
  ### output
  res <- list(series=serie, X=XX, ww.index=ww.index, ww.values=ww.values, px.neighbors=px.neighbors, n.NA=data.frame(n.NA), 
              px=px, lon=longit, lat=latit, group=gruppi)
  res
}
