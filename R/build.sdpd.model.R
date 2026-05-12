build.sdpd.model <- function(endogenous, lambda0=T, lambda1=T, lambda2=T, 
                             covariates=NULL, fixed.effects=T, time.effects=F, 
                             var.names=NULL, coeffs=NULL, pp=NULL){
  # questa funzione costruisce la struttura del modello SDPD così come atteso dalle funzioni di stima
  # l'unico parametro obbligatorio è la variabile endogena, che indica:
  #     - se numerico, l'indice del vettore nomi (che deve in tal caso essere presente)
  #     - se stringa, il nome di variabile (anche in assenza del parametro nomi)

  model <- list()
  
  # verifica del formato di nomi.variabili
  if((!is.null(var.names)) & (!is.character(var.names) | !is.vector(var.names)))
    return(list(error="L'attributo var.names non è del formato atteso (character)"))
  
  # verifica della consistenza del parametro coeffs e definizione nomi variabili
  if(!is.null(coeffs) & (is.data.frame(coeffs)|is.matrix(coeffs))){
    nomi <- dimnames(coeffs)[[2]]
    if(is.null(nomi))
      return(list(error="L'attributo names della matrice dei coefficienti deve essere definito"))
    if("lambda0"%in%dimnames(coeffs)[[2]]){
      if(!lambda0)
        coeffs <- coeffs[-c("lambda0")]
    }
    else lambda0 <- FALSE
    if("lambda1"%in%dimnames(coeffs)[[2]]){
      if(!lambda1)
        coeffs <- coeffs[-c("lambda1")]
    }
    else lambda1 <- FALSE
    if("lambda2"%in%dimnames(coeffs)[[2]]){
      if(!lambda2)
        coeffs <- coeffs[-c("lambda2")]
    }
    else lambda2 <- FALSE
    model$coeffs <- as.matrix(coeffs)
    model$pp <- dim(coeffs)[1]
  }
  else if(!is.null(coeffs) & is.vector(coeffs)  & is.numeric(coeffs)){
    nomi <- names(coeffs)    
    if(is.null(nomi))
      return(list(error="L'attributo names del vettore dei coefficienti non può essere NULL"))
    if(is.null(pp))
      return(list(error="In presenza di un vettore di coefficienti il parametro pp non può essere NULL"))
    else if(!numeric(pp) | length(pp)>1 | pp<0)
      return(list(error="Il parametro pp non è del formato atteso (valore positivo)"))
    coeffs <- matrix(rep(coeffs, pp), nrow=pp, byrow = T)
    if("lambda0"%in%dimnames(coeffs)[[2]]){
      if(!lambda0)
        coeffs <- coeffs[-c("lambda0")]
    }
    else lambda0 <- FALSE
    if("lambda1"%in%dimnames(coeffs)[[2]]){
      if(!lambda1)
        coeffs <- coeffs[-c("lambda1")]
    }
    else lambda1 <- FALSE
    if("lambda2"%in%dimnames(coeffs)[[2]]){
      if(!lambda2)
        coeffs <- coeffs[-c("lambda2")]
    }
    else lambda2 <- FALSE
    model$coeffs <- coeffs
    model$pp <- dim(coeffs)[1]
  }
  else if(is.null(coeffs)){
    nomi <- NULL
    # le seguenti condizioni sui coefficienti lambda, apparentemente superflue, permettono di verificare che vengano passati valori booleani
    if(lambda0) 
      lambda0 <- TRUE
    if(lambda1) 
      lambda1 <- TRUE
    if(lambda2) 
      lambda2 <- TRUE
  }
  else return(list(error="I coefficienti passati non sono del formato atteso (dataframe, matrice o vettore numerico)"))

  # definizione componenti lambda
  if(lambda0+lambda1+lambda2==0)
    return(list(error="Deve esserci almeno un componente lambda attivo"))
  model$lambda_coeffs <- c(lambda0, lambda1, lambda2)
  names(model$lambda_coeffs) <- c("lambda0", "lambda1", "lambda2")
  
  # controllo nome/indice variabile endogena 
  if(is.numeric(endogenous)|is.character(endogenous)){
    if(length(endogenous)>1)
      return(list(error="La variabile endogena deve essere unica, non è ammesso più di un valore"))
    if(!is.null(var.names)){
      if(is.character(endogenous) & sum(endogenous%in%var.names)==0)
        return(list(error="La variabile endogena non è presente nel vettore di variabili"))
      if(is.numeric(endogenous) & endogenous>length(var.names))
        return(list(error="La variabile endogena non è presente nel vettore di variabili"))
      else if(is.numeric(endogenous))
        endogenous <- var.names(endogenous)
    }
    model$name_endogenous <- endogenous
  }
  else return(list(error="La variabile endogena non è del formato atteso (numerico o carattere)"))
  

  # controllo nomi/indici variabili esogene e compatibilità con la variabile endogena
  if((is.numeric(covariates)|is.character(covariates)) & is.vector(covariates)){
    if(length(covariates)==0){
      model$kk <- 0
    }
    else{
      if(endogenous%in%covariates)
        return(list(error="La variabile endogena non può essere inclusa anche tra le esogene"))
      if(!is.null(var.names) & is.character(covariates)){
        if(sum(covariates%in%var.names)<length(covariates))
          return(list(error="Uno o più indici di variabile esogena non presenti nel vettore di variabili"))
      }
      if(!is.null(var.names) & is.numeric(covariates)){
        if(max(covariates)>length(var.names))
          return(list(error="Uno o più indici di variabile esogena non presenti nel vettore di variabili"))
        covariates <- var.names(covariates)
      }
      if(!is.null(nomi)){
        if(sum(covariates%in%nomi)<length(covariates))
          return(list(error="Uno o più indici di variabile esogena non compatibili con i dati passati"))
      }
      model$kk <- length(covariates)
      model$beta_coeffs <- rep(T, model$kk)
      names(model$beta_coeffs) <- covariates
    }
  }
  else return(list(error="Le variabili esogene non sono del formato atteso (vettore numerico o vettore carattere)"))

  # definizione della componente effetti fissi (con contestuale controllo del formato booleano)
  if(fixed.effects){
    if((!is.null(nomi)) & ("fixed_effects"%in%nomi)){
        model$fixed_effects <- TRUE
    }
    else if(is.null(nomi)){
      model$fixed_effects <- TRUE
    }
    else return(list(error="La matrice dei coefficienti non contiene gli effetti fissi, sebbene siano previsti nel modello"))
  }
  else{
    model$fixed_effects <- FALSE
    if((!is.null(model$coeffs)) & ("fixed_effects"%in%nomi))
        model$coeffs <- model$coeffs[,-c("fixed_effects")]
  }

  # definizione della componente effetti temporali (con contestuale controllo del formato)
  if(!is.null(nomi)){
    if(is.logical(time.effects) & length(time.effects)==1)
      model$time_effects <- time.effects
    else if(is.numeric(time.effects) & is.vector(time.effects))
      model$time_effects <- time.effects
    else return(list(error="Il parametro time.effects non è del formato previsto (valore booleano oppure vettore numerico)"))
  }
  else{
    if(is.logical(time.effects) & length(time.effects)==1)
      model$time_effects <- time.effects
    else return(list(error="Il parametro time.effects non è del formato previsto (valore booleano)"))
  }
  
  ## returning the built model
  model
}
