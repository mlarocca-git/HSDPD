fit.sdpd.covs <- function(series, X=NULL, kk=0, px.neighbors=NULL, nn, pp, NAcovs="pairwise.complete.obs"){
  if(is.null(px.neighbors)){
    seriesBoundary <- series
    index <- matrix(rep(dimnames(series)[[1]], pp), byrow = T, nrow = pp)
    dimnames(index)[[1]] <- dimnames(series)[[1]]
  }
  else{
    seriesBoundary <- rbind(series, px.neighbors$seriesBoundary)
    index <- px.neighbors$index
  }
  cov12 <- cov(t(series)[2:nn,],t(seriesBoundary)[-nn,], use=NAcovs)
  cov11 <- cov(t(series), t(seriesBoundary), use=NAcovs)
  if(kk==0 | is.null(kk))
    covX <- NULL
  else if(kk==1)
    covX <- cov(t(X)[2:nn,],t(seriesBoundary)[-nn,], use=NAcovs)
  else if(kk>1){
    covX <- array(0, dim=c(dim(X)[1], dim(X)[2], dim(seriesBoundary)[1]))
    dimnames(covX)[[1]] <- dimnames(X)[[1]]
    dimnames(covX)[[2]] <- dimnames(X)[[2]]
    dimnames(covX)[[3]] <- dimnames(seriesBoundary)[[1]]
    for(jj in 1:kk)
      covX[jj,,] <- cov(t(X[jj,,2:nn]),t(seriesBoundary)[-nn,], use=NAcovs)
  }
  COVs <- list(index=index, cov11=cov11, cov12=cov12, covX=covX)
}

