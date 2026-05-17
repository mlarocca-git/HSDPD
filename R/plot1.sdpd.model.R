plot1.sdpd.model <- function (res.fit, n.units="all", n.vars="all", t.axis=list(t.labels=NULL, t.points=NULL), xlimit=NULL, ylimit=NULL, max.col=5, col.punti=1) 
{
	## res.fit is an object with the following components:
	## "coeff.hat" "series"    "WW"        "XX"        "fitted"    "resid"     "errors"    "warnings"  "model" 
  
  ## definizioni variabili e dimensioni
  dserie <- res.fit$data$series
  XX <- res.fit$data$X
	nn <- dim(dserie)[2]
	pp <- dim(dserie)[1]
	alpha.hat <- res.fit$coeff.hat
	error <- NULL

	## definizione pannelli da rappresentare graficamente
	panels <- dimnames(res.fit$coeff.hat)[[2]]
	i.panels <- rep(FALSE, length(panels))
	if(is.logical(n.vars)){
    nn.p <- min(length(n.vars), length(i.panels))
	 	i.panels[1:nn.p] <- n.vars[1:nn.p]
	}
	else if(is.numeric(n.vars)&length(n.vars)==1){
	  nn.p <- min(length(i.panels), as.integer(n.vars))
	  i.panels[1:nn.p] <- TRUE
	}
	else if(is.numeric(n.vars)){
	  check <- as.integer(n.vars)%in%seq(1,length(i.panels))
	  n.vars <- as.integer(n.vars)[check]
	  nn.p <- min(length(n.vars), length(i.panels))
	  i.panels[n.vars[1:nn.p]] <- TRUE
	}
	else if(n.vars[1]=="all"){
	  n.vars <- dimnames(res.fit$coeff.hat)[[2]]
	  i.panels <- panels%in%n.vars
	}
	else if(is.character(n.vars)){
	  i.panels <- panels%in%n.vars
	}
	else return(error="Il valore passato all'argomento n.vars non segue il formato previsto")
	
	panels <- panels[i.panels]
	if(length(which(panels=="fixed_effects"))>0)
  	panels <- panels[-which(panels=="fixed_effects")]
	n.col <- min(max.col, length(panels))
  if(n.units[1]=="all")
    n.units <- seq(1, dim(res.fit$data$series)[1])
	if(is.character(n.units)|is.numeric(n.units)){
	  indici <- seq(1, dim(res.fit$data$series)[1])
	  names(indici) <- dimnames(res.fit$data$series)[[1]]
	  n.units <- indici[n.units]
	}
	else return(error="Il valore passato all'argomento n.units non segue il formato previsto")
	if(res.fit$model$fixed_effects)
	  fixed_effects.text <- "fixed effect (dashed grey) & "
	else
	  fixed_effects.text <- ""

		## ciclo principale per la costruzione dei grafici
	op <- par(no.readonly = TRUE)
	for(ii in n.units){
	  if(sum(is.na(res.fit$coeff.hat[ii,]))>0){
	    cat("\n There are missing estimated coefficients for this series ")
	    next
	  }
	  beta.cont <- 0
	  par(mfcol=c(2,n.col), mai=c(1, 0.8, 0.5, 0.5), las=2)
		for(panel in panels[1:n.col]){
		  xlimit <- ylimit <- NULL
		  if(panel=="lambda0") jj <- 1
		  else if(panel=="lambda1") jj <- 2
		  else if(panel=="lambda2") jj <- 3
		  else{
		    jj <- 4
		    beta.cont <- beta.cont+1
		    if(is.matrix(res.fit$data$X)){
		      titolo1 <- substitute(bold(str0*str1)*str2*hat(beta)[i] == y, list(str0="i=", str1=dimnames(dserie)[[1]][ii], str2=": ", y = round(alpha.hat[ii,panel], digits=6)))
		      XX <- res.fit$data$X
		      nomeX <- "X"
		    }
		    else if(is.array(res.fit$data$X)){
		      titolo1 <- substitute(bold(str0*str1)*str2*hat(beta)[i] == y, list(str0="i=", str1=dimnames(dserie)[[1]][ii], str2=": ", y = round(alpha.hat[ii,panel], digits=6)))
		      XX <- res.fit$data$X[panel,,]
		      nomeX <- panel
		    }
		  }
		  if(jj==4){
		    lag <- 0
		    Wdserie <- XX
			  yll <- expression(y[i*t])
			  xll <- substitute(str0[str1*str2], list(str0=nomeX, str1="i", str2="t"))
			  plot.title <- paste("Left axis: observed series (black) & ", nomeX, "'s signal (green)\n Right axis: ", fixed_effects.text, "residuals (solid grey)", sep="")
			  off.color <- 1
		  }
			else if(jj==3){
			  lag <- 1
				Wdserie <- res.fit$data$W%*%dserie
				yll <- expression(y[i*t])
				xll <- expression(bold(w)[i]*bold(y)[t-1])
				plot.title <- paste("Left axis: observed series (black) & spatial-dynamic signal (yellow)\n Right axis: ", fixed_effects.text, "residuals (solid grey)", sep="")
				titolo1 <- substitute(bold(str0*str1)*str2*hat(lambda)[x*i] == y, list(str0="i=", str1=dimnames(dserie)[[1]][ii], str2=": ", x = 2, y = round(alpha.hat[ii,panel], digits=6)))
				off.color <- 5
			}
			else if(jj==2){
			  lag <- 1
				Wdserie <- dserie
				yll <- expression(y[i*t])
				xll <- expression(y[i*(t-1)])
				plot.title <- paste("Left axis: observed series (black) & pure dynamic signal (blue)\n Right axis: ", fixed_effects.text, "residuals (solid grey)", sep="")
				titolo1 <- substitute(bold(str0*str1)*str2*hat(lambda)[x*i] == y, list(str0="i=", str1=dimnames(dserie)[[1]][ii], str2=": ", x = 1, y = round(alpha.hat[ii,panel], digits=6)))
				off.color <- 2
			}
			else{
			  lag <- 0
				Wdserie <- res.fit$data$W%*%dserie
				yll <- expression(y[i*t])
				xll <- expression(bold(w)[i]*bold(y)[t])
				plot.title <- paste("Left axis: observed series (black) & pure spatial signal (red)\n Right axis: ", fixed_effects.text, "residuals (solid grey)", sep="")
				titolo1 <- substitute(bold(str0*str1)*str2*hat(lambda)[x*i] == y, list(str0="i=", str1=dimnames(dserie)[[1]][ii], str2=": ", x = 0, y = round(alpha.hat[ii,panel], digits=6)))
				off.color <- 0
			}
			#### plot of spatial regression
		  if(is.null(xlimit))
				xlimit <- range(Wdserie[ii,], na.rm=T)
			if(is.null(ylimit))
				ylimit <- range(dserie[ii,], na.rm=T)
		  plot(Wdserie[ii,1:(nn-lag)], dserie[ii,(1+lag):nn], ylim=ylimit, xlim=xlimit, xlab="", ylab="", cex.lab=1.5, cex.axis=0.9, col=col.punti)
			title(xlab=xll, ylab=yll, cex.lab=1.5, cex.main=1.5, main=titolo1)
			lines(xlimit, xlimit*alpha.hat[ii,panel], col=2+off.color, lwd=1, lty=1)
			#### plot of spatial and residual estimated series
			fun.asse2 <- function(x, ylim1, ylim2){
			  y <- ylim1[1]+(x-ylim2[1])*range(ylim1)/range(ylim2)
			  y
			}
			ylimitb <- range(dserie[ii,], alpha.hat[ii,panel]*Wdserie[ii,1:(nn-lag)],na.rm=T)
			if(res.fit$model$fixed_effects)
			  ylimitc <- range(res.fit$resid[ii,], res.fit$coeff.hat[ii,"fixed_effects"], na.rm=T)
			else
			  ylimitc <- range(res.fit$resid[ii,], na.rm=T)
			dum <- 0
			if(abs(ylimitc[2]-ylimitb[2])<diff(ylimitb)*0.1 & abs(ylimitc[1]-ylimitb[1])<diff(ylimitb)*0.1)
			  ylimitb <- range(ylimitb, ylimitc)
			else if(ylimitc[2]<ylimitb[1]-diff(ylimitb)*0.1){
			  ylimitb[1] <- ylimitb[1]-diff(ylimitb)*0.1
			}
			else if(ylimitc[1]>ylimitb[2]+diff(ylimitb)*0.1){
			  ylimitb[2] <- ylimitb[2]+diff(ylimitb)*0.1
			}
			ts.plot(dserie[ii,(1+lag):nn], ylim=ylimitb, xlab="", ylab="", gpars = list(axes=F, cex.main=0.9), type="n", main=paste(plot.title))
			lines(seq(1, nn), fun.asse2(res.fit$resid[ii,], ylim1=ylimitc, ylim2=ylimitc), col="grey")
			abline(h=fun.asse2(0, ylim1=ylimitc, ylim2=ylimitc), col="grey")
			if(res.fit$model$fixed_effects)
			  abline(h=fun.asse2(res.fit$coeff.hat[ii,"fixed_effects"], ylim1=ylimitc, ylim2=ylimitc), col="grey", lty=2)
			lines(seq(1+lag, nn), dserie[ii,(1+lag):nn])
			lines(seq(1+lag, nn), alpha.hat[ii,panel]*Wdserie[ii,1:(nn-lag)], col=2+off.color, lwd=1)
			axis(2)
			#if(res.fit$model$fixed_effects)
			#  axis(4, at=c(seq(0, ylimitb[2], by=10), fun.asse2(res.fit$coeff.hat[ii,"fixed_effects"])), labels=c(seq(0, ylimitb[2], by=10), round(res.fit$coeff.hat[ii,"fixed_effects"], 2)), cex.axis=0.6, col.ticks="grey")
			#else
			#  axis(4, at=c(seq(0, ylimitb[2], by=10))-dum*ylimitc[1]+dum*ylimitb[1], labels=c(seq(0, ylimitb[2], by=10)), cex.axis=0.6, col.ticks="grey")
			if(!is.null(t.axis$t.points)){
				if(length(t.axis$t.points)!=length(t.axis$t.labels))
				error <- error + 1
				t.axis$t.points <- NULL
			}
			if(is.null(t.axis$t.points)){
			  t.axis$t.points <- seq(1, nn, by=nn%/%15)
			  t.axis$t.labels  <- dimnames(res.fit$data$series)[[2]][t.axis$t.points]
			}
			axis(1, at=t.axis$t.points, labels=t.axis$t.labels, cex.axis=0.8)
			box()
		}
	  cat("\n unit=", ii)
	  if(!ii==n.units[length(n.units)]){
	    cat("\n click on the plot window for the next panel....")
	    pos <- locator(1)
	  }
	  else
	    cat("\n panels are finished...")
	}
	if(length(error)>0)
  	list(error=error)
}
