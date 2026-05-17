build.sdpd.model <- function(endogenous="var_y", lambda0=T, lambda1=T, lambda2=T, 
                            covariates=0, fixed.effects=T, time.effects=F, 
                            coeffs=NULL, check=TRUE,
                            ww.index=NULL, ww.values=NULL, px.neighbors=NULL,
                            rrgroups=NULL, 
                            SIM=FALSE, pp=NULL, d.i=100,
                            sim.options=list(ww.type="queen", px.neighbors=8, index.group=rep(1, pp), 
                                             fixed.effects=list(constant=F, min=10, max=20),
                                             lambda0=list(constant=F, min=-0.5, max=0.5),
                                             lambda1=list(constant=F, min=-0.5, max=0.5),
                                             lambda2=list(constant=F, min=-0.5, max=0.5),
                                             betas=list(constant=F, min=c(0,0,0,1), max=c(1,1,1,4)),
                                             sigma.eps=list(constant=T, min=1, max=1))){
  
  ## Questa funzione costruisce l'oggetto sdpd-model, che sarà il principale input alle successive funzioni di stima;
  ## se SIM=TRUE, allora tutte le componenti del modello vengono generate in modo casuale, utilizzando i parametri passati in sim.options;
  ## se check=TRUE, allora effettua a posteriori anche una verifica della struttura e della stazionarietà del modello ridotto

  model <- list()

  # le seguenti tre condizioni, apparentemente superflue, permettono di verificare che vengano passati valori booleani
  if(lambda0) 
    lambda0 <- TRUE
  if(lambda1) 
    lambda1 <- TRUE
  if(lambda2) 
    lambda2 <- TRUE
  # definizione componenti lambda
  if(lambda0+lambda1+lambda2==0)
    return(list(error="Deve esserci almeno un componente lambda attivo"))
  model$lambda_coeffs <- c(lambda0, lambda1, lambda2)
  names(model$lambda_coeffs) <- lambda.names <- c("lambda0", "lambda1", "lambda2")
  lambda.names <- lambda.names[model$lambda_coeffs]
  
  # inizializzazione nomi di variabili eventualmente passati mediante la matrice coeffs
  if(is.null(coeffs))
    nomi <- NULL
  else if(is.data.frame(coeffs) | is.matrix(coeffs)){
    coeffs <- as.matrix(coeffs)
    nomi <- dimnames(coeffs)[[2]]
    if(is.null(nomi))
      return(list(error="La matrice dei coefficienti deve avere i nomi di colonna definiti"))
  }
  else return(list(error="Il parametro coeffs non è del formato atteso (matrice o dataframe di numeri)"))

  # controllo nome/indice variabile endogena 
  if(is.vector(endogenous) & is.character(endogenous)){
    if(length(endogenous)>1)
      return(list(error="La variabile endogena deve essere unica, non è ammesso più di un valore"))
    if(sum(endogenous%in%covariates)>0)
      return(list(error="La variabile endogena non può risultare tra le covariate"))
    model$name_endogenous <- endogenous
  }
  else return(list(error="La variabile endogena non è del formato atteso (carattere)"))
  
  # controllo nomi/indici variabili esogene
  if(is.vector(covariates) & (is.numeric(covariates)|is.character(covariates))){
    if(covariates[1]==0){
      model$beta_coeffs <- beta.names <- NULL
    }
    else if(sum(covariates%in%nomi)>0){
      beta.names <- covariates[covariates%in%nomi]
      if(length(beta.names)<length(covariates))
        cat("\nWarning: alcune covariate sono state eliminate dal modello perché non presenti nella matrice coeffs\n")
      model$beta_coeffs <- rep(T, length(beta.names))
      names(model$beta_coeffs) <- beta.names
    }
    else if(is.null(nomi)){
      model$beta_coeffs <- rep(T, length(covariates))
      beta.names <- names(model$beta_coeffs) <- covariates        
    } 
  }
  else return(list(error="Il parametro delle variabili esogene non è del formato atteso (vettore numerico o carattere)"))
  
  # definizione della componente effetti fissi (con contestuale controllo del formato booleano)
  if(!is.null(nomi) & fixed.effects){
    if("fixed_effects"%in%nomi){
        model$fixed_effects <- TRUE
        fixed.name <- "fixed_effects"
    }
    else{
      model$fixed_effects <- FALSE
      cat("\nWarning: gli effetti fissi sono stati eliminati dal modello poiché non presenti nella matrice coeffs.\n")
      fixed.name <- NULL
    }
  }
  else{
    model$fixed_effects <- fixed.effects
    fixed.name <- c("fixed_effects")[fixed.effects]
  }

  # definizione della componente effetti temporali (con contestuale controllo del formato)
  if(is.logical(time.effects) & length(time.effects)==1)
    model$time_effects <- time.effects
  else return(list(error="Il parametro time.effects non è del formato previsto (valore booleano)"))

  ## in caso di modello simulato
  if(SIM){
    if(is.null(pp))
      return(list(error="In caso di coefficienti da simulare, è necessario conoscere il valore della dimensione pp"))
    if(is.null(rrgroups)){
      n.groups=pp%/%d.i
      gruppi <- rep(n.groups, pp)
      gruppi[1:length(rep(seq(1, n.groups), each=pp%/%n.groups))] <- rep(seq(1, n.groups), each=pp%/%n.groups)
      labels <- paste("group", gruppi, sep="_")
      gruppi <- cbind(COD=gruppi, LABEL=labels)
      dimnames(gruppi)[[1]] <- seq(1, pp)
    }
    new.model <- generate.sdpd.model(pp=pp, model=model, ww.index=ww.index, ww.values=ww.values, px.neighbors=px.neighbors,
                                     sim.options=sim.options) 
    model$coeffs <- new.model$coeffs
    model$ww.index <- new.model$ww.index
    model$ww.values <- new.model$ww.values
    model$px.neighbors <- new.model$px.neighbors
    model$sigma.eps <- new.model$sigma.eps
    model$groups <- gruppi
  }
  
  ## diagnostics
  if(!is.null(model$coeffs) & !is.null(model$ww.index) & !is.null(model$ww.values) & check){
    model$diagnostics <- check.sdpd.model(model=model)$diagnostics[1,]
  }
  
  ## output
  model
}
