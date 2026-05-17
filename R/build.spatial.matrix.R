build.spatial.matrix <- function(ww.index, ww.values){
  
  ## this function built the square spatial matrix WW
  pp <- dim(ww.index)[1]
  WW <- matrix(0, nrow=pp, ncol=pp)
  dimnames(WW)[[1]] <- dimnames(WW)[[2]] <- dimnames(ww.index)[[1]]
  for(ii in 1:pp){
    indici <- as.character(ww.index[ii,ww.index[ii,]>0])
    WW[ii, indici] <- ww.values[ii,ww.index[ii,]>0]
  }
  
  ## checking the consistency of the spatial matrix  
  if(sum(is.na(WW))>0)
    return(error="The spatial matrix cannot have NA values.")
  if(sum(diag(abs(WW)))>0)
    return(error="The spatial matrix has not zero diagonal.")
  indici.w <- apply(abs(WW), 1, sum)==0
  if(sum(indici.w)>0)
    return(error="The spatial matrix has one (or more) zero-row(s).")
  
  ## output
  WW
}