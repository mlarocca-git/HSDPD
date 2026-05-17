check.sdpd.model <- function(res.fit=NULL, model=NULL, ww.index=NULL, ww.values=NULL){
  ## questa funzione permette di controllare la validità del modello o del modello stimato
  ## valutando le condizioni di stazionarietà
  
  if(!is.null(res.fit)){
    coeffs <- res.fit$coeff.hat
    if(!is.null(model))
      cat("\nAttenzione....è stato valutato il modello ereditato da res.fit (l'oggetto model è stato quindi ignorato).")
    model <- res.fit$model
  }
  else coeffs <- NULL
  
  if(!is.null(model)){
    if(is.null(coeffs))
      coeffs <- model$coeffs
    px <- dimnames(coeffs)[[1]]
  }
  
  ## checking validity of the spatial matrix
  if(is.null(ww.index) | is.null(ww.values)){
    ww.index <- model$ww.index[px,]
    ww.values <- model$ww.values[px,]
  }
  if(is.null(ww.index) | is.null(ww.values))
    return(error="The spatial matrix components are missing")
  else{
    WW <- build.spatial.matrix(ww.index=ww.index, ww.values=ww.values)
  }
  
  if(!is.null(coeffs)){
    pp <- dim(coeffs)[1]
    ## calcoliamo autovalori delle componenti della matrice ridotta A
    ll0 <- ll1 <- ll2 <- matrix(0, ncol=pp, nrow=pp)
    vettore1 <- vettore2 <- eigenA <- rep(0, pp)
    if("lambda0"%in%dimnames(coeffs)[[2]]){
      ll0 <- diag(coeffs[,"lambda0"])%*%WW
      matrice1 <- solve(diag(rep(1, pp))-ll0)
      if(sum(is.na(matrice1))==0)
        vettore1 <- eigen(matrice1)$values
      else 
        vettore1 <- rep(NA, pp)
    }
    else matrice1 <- diag(vettore1)

    if("lambda1"%in%dimnames(coeffs)[[2]])
      ll1 <- diag(coeffs[,"lambda1"])

    if("lambda2"%in%dimnames(coeffs)[[2]])
      ll2 <- diag(coeffs[,"lambda2"])%*%WW

    matrice2 <- ll1+ll2
    if(sum(is.na(matrice2))==0)
      vettore2 <- eigen(matrice2)$values
    else
      vettore2 <- rep(NA, pp)
    if(sum(is.na(matrice1%*%matrice2))==0)
      eigenA <- eigen(matrice1%*%matrice2)$values
    else
      eigenA <- rep(NA, pp)
    diagnostics <- cbind(max.mod.eigen1=Mod(vettore1), max.mod.eigen2=Mod(vettore2), max.mod.eigenA=Mod(eigenA))
  }
  else cat("\nNo coefficients to evaluate....")
  
  ## output
  list(diagnostics=data.frame(diagnostics))
}
