fit.sdpd.assemble <- function(obj.stime){
  if(is.null(obj.stime)){
    cat("\n A null object has been received.")
    res <- NULL
  }
  else if(is.null(obj.stime$error)){
    indici <- as.character(obj.stime$px)
    px <- obj.stime$px
    names(px) <- indici
    vec.mean <- apply(obj.stime$resid, 1, mean, na.rm=T)
    vec.sd <- apply(obj.stime$resid, 1, sd, na.rm=T)
    vec.LB <- apply(obj.stime$resid, 1, fun.LBtest)
    vec.LB <- p.adjust(vec.LB, method = "BY")
    vec.JB <- apply(obj.stime$resid, 1, fun.JBtest)
    vec.JB <- p.adjust(vec.JB, method = "BY")
    eigenA <- round(rep(obj.stime$diagnostics[1], dim(obj.stime$resid)[1]), 2)
    names(vec.mean) <- names(vec.sd) <- names(vec.LB) <- names(vec.JB) <- names(eigenA) <- dimnames(obj.stime$resid)[[1]]

    res <- list(px=px, lon=obj.stime$lon, lat=obj.stime$lat, group=data.frame(obj.stime$group)[indici,],
                  fitted=data.frame(obj.stime$fitted, check.names=F)[indici,], resid=data.frame(obj.stime$resid, check.names=F)[indici,],
                  coeff.hat=data.frame(obj.stime$coeff.hat, check.names=F)[indici,],
                  mean.resid=vec.mean[indici], sd.resid=vec.sd[indici], pvalue.LB=vec.LB[indici], pvalue.JB=vec.JB[indici],
                  max.eigenA=eigenA[indici],
                  mu.means=data.frame(obj.stime$mean.equation$mus, check.names=F)[indici,],
                  mu.tmeans=data.frame(obj.stime$mean.equation$mus_t_weighted, check.names=F)[indici,])
  }
  else{
    cat("\n There was an error in the estimation:", obj.stime$group[1,], obj.stime$error)
    res <- NULL
  }
  res
}
