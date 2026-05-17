plot.sdpd.boot.series <- function(bb, pp, nn, eigen, yy.star1, opts.boot, res.fit, group.index){
  if(is.null(opts.boot$ylimiti))
    opts.boot$ylimiti <- range(res.fit$data$series)
  cat("\n boot iterations:", bb-300+1, "-", bb, "...one of the bootstrapped time series is plotted for control")
  jji <- sample((1:pp)[group.index==1], size=1)
  titolo <- "True series (red) vs bootstrap series (black)"
  titoloB <- NULL
  titolo <- paste(titolo, " - location=", dimnames(res.fit$data$series)[[1]][jji], sep="")
  if(!is.null(opts.boot$label.index))
    titolo <- paste(titolo, " - group=", opts.boot$label.index[jji], sep="")
  
  titolo <- paste(titolo, "\n core.points=", sum(abs(res.fit$data$W[jji,])>0), ", n.points=", sum(res.fit$data$px.neighbors$index[jji,]>0, na.rm=T), sep="")
  titolo <- paste(titolo, "  (", sep="")
  l.sum <- 0
  if("lambda0"%in%dimnames(res.fit$coeff.hat)[[2]]){
    titolo <- paste(titolo, "l0=", round(res.fit$coeff.hat[jji,"lambda0"], digits=2), sep="")
    l.sum <- l.sum + res.fit$coeff.hat[jji,"lambda0"]
  }
  if("lambda1"%in%dimnames(res.fit$coeff.hat)[[2]]){
    titolo <- paste(titolo, ", l1=", round(res.fit$coeff.hat[jji,"lambda1"], digits=2), sep="")
    l.sum <- l.sum + res.fit$coeff.hat[jji,"lambda1"]
  }
  if("lambda2"%in%dimnames(res.fit$coeff.hat)[[2]]){
    titolo <- paste(titolo, ", l2=", round(res.fit$coeff.hat[jji,"lambda2"], digits=2), sep="")
    l.sum <- l.sum + res.fit$coeff.hat[jji,"lambda2"]
  }
  titolo <- paste(titolo, titoloB, ") eigen=", round(eigen[1], digits=2), sep="")
  if("lambda1"%in%dimnames(res.fit$coeff.hat)[[2]])
    titolo <- paste(titolo, ", #{|l1|>1}=", sum(ifelse(abs(res.fit$coeff.hat[,"lambda1"])>1, 1, 0)), sep="")
  titolo <- paste(titolo, ", #{|l0+l1+l2|>1}=", sum(ifelse(abs(l.sum)>1, 1, 0)), sep="")
  titolo <- paste(titolo, ", DEVS=(", round(mean(yy.star1[jji,]-res.fit$data$series[jji,]), 2), ";", round(sd(yy.star1[jji,]-res.fit$data$series[jji,]), 2), ")", sep="")
  nome.tsplot <- paste(opts.boot$cartel, opts.boot$label.index[jji], dimnames(res.fit$data$series)[[1]][jji], "bootplot_",bb , ".jpeg", sep="")
  sottotitolo <- paste(res.fit$model$name_endogenous, "~H-SDPD  (based on data from ",
                       dimnames(res.fit$data$series)[[2]][1], " to ", dimnames(res.fit$data$series)[[2]][2], ")", sep="")
  jpeg(file = nome.tsplot, width = 1100, height = 600)
  ts.plot(yy.star1[jji,], col=group.index[jji], ylim=opts.boot$ylimiti, ylab="temperatures")
  title(main=titolo, sub=sottotitolo, cex.main=1.5)
  lines(seq(1,nn), res.fit$data$series[jji,], col="red")
  abline(h=res.fit$coeff.hat[jji,"fixed_effects"])
  if(kk==1)
    lines(seq(1,nn), res.fit$data$X[jji,], col="blue", lty=2)
  else if(kk>1)
    lines(seq(1,nn), res.fit$data$X[1,jji,], col="blue", lty=2)
  dev.off()
}

