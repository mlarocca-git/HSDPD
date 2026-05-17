test.sdpd.model <- function (res.fit, px, model=NULL, n.boot=399,
                             H0=c("zero", "constant", "grouped", "nospatial", "noautoregressive", "noX", "constrained")[1],
                             coeff.H0=0,
                             group.index=rep(1, dim(res.fit$fitted)[1]),
                             opts.boot=list(markovian=FALSE, resid=c("fitted", "normal")[1], sigma.resid=1, boot.plot=FALSE, cartel="", ylimiti=NULL, label.index=NULL))
{
  ## H0=("constant", "grouped", ....) denotes the kind of hypothesis tested under the null.
  if(is.null(model))
    model.obj <- res.fit$model
  else
    model.obj <- check.sdpd.model(model=model)$model

  ## checking stationarity conditions
  data <- check.sdpd.model(res.fit=res.fit)
  if(length(data$errors)>0){
    cat("\n There are some errors....\n", data$errors)
    return(list(errors=data$errors, warnings=data$warnings))
  }
  if(length(data$warnings)>0){
    cat("\n There are some warnings....\n", data$warnings)
  }
  
  ### sample distribution approximation based on residual bootstrap
  if(n.boot>20){
    if(opts.boot$resid=="fitted"){
      resid <- res.fit$resid[,-1]
    }
    else if(opts.boot$resid=="normal"){
      if(length(opts.boot$sigma.resid)==1)
        vec.sigma <- rep(opts.boot$sigma.resid, data$pp)
      else if(length(opts.boot$sigma.resid)>=data$pp)
        vec.sigma <- opts.boot$sigma.resid[1:data$pp]
      else vec.sigma <- rep(1, data$pp)
      resid <- t(apply(as.matrix(vec.sigma), 1, FUN=function(xx,nn){rnorm(n=nn, sd=xx)}, nn=data$nn))
    }
    i.resid	<- seq(1,dim(resid)[2])
    boot.ind 	<- matrix(sample(i.resid, size=(data$nn)*n.boot, replace = TRUE), nrow=n.boot)
    boot.rep	<- array(0, dim=c(n.boot, data$pp, dim(res.fit$coeff.hat)[2]))
    dimnames(boot.rep)[[3]] <- dimnames(res.fit$coeff.hat)[[2]]
    devs.tsboot <- array(0, dim=c(n.boot, 2, data$pp))
    ## bootstrap iterations...
    for(bb in 1:n.boot){
      resid.boot <- resid[,boot.ind[bb,]]
      yy.star1 <- res.fit$data$series
      for(step in 1:5)
        yy.star1 <- fit.sdpd.series(dseries=yy.star1, W=res.fit$data$W, X.centr=res.fit$data$X, model=model.obj, coeff.hat=res.fit$coeff.hat, resid=resid.boot, markovian=opts.boot$markovian)$series
      COVs <- fit.sdpd.covs(series=yy.star1, X=res.fit$data$X, px.neighbors=res.fit$data$px.neighbors, kk=data$kk, nn=data$nn, pp=data$pp)
      boot.rep[bb,,] <- fit.sdpd.coefficients(W=res.fit$data$W, COVs=COVs, mu=data$mu, model=model.obj)$coeff.hat
      devs.tsboot[bb,,] <- apply(yy.star1-res.fit$data$series, 1, FUN=function(x){c(mean=mean(x), sd=sd(x))})
      if(bb%%300==0 & opts.boot$boot.plot)
        plot.sdpd.boot.series(bb=bb, pp=data$pp, nn=data$nn, eigen=data$diagnostics$Mod.eigenA[1], yy.star1=yy.star1, opts.boot=opts.boot, res.fit=res.fit, group.index=group.index)
    }
    ## collecting results...
    results <- coeff.constr <- array(0, dim=c(data$pp, dim(boot.rep)[3]))
    dimnames(results)[[1]] <- dimnames(coeff.constr)[[1]] <- dimnames(boot.rep)[[2]] <- dimnames(res.fit$coeff.hat)[[1]]
    dimnames(results)[[2]] <- dimnames(coeff.constr)[[2]] <- dimnames(res.fit$coeff.hat)[[2]]

    ### definizione delle statistiche per le varie tipologie di ipotesi del test...
    if(H0=="zero"|H0==1){
      boot.constr <- boot.rep
      boot.constr[,,] <- 0
      coeff.constr[,] <- coeff.H0
    }
    else if(H0=="constant"|H0==2){
      boot.constr <- apply(boot.rep, c(1,3), FUN=function(x){rep(mean(x), length(x))})
      boot.constr <- apply(boot.constr, c(1,3), FUN=function(x){x})
      coeff.constr[,] <- apply(res.fit$coeff.hat, 2, FUN=function(x){rep(mean(x), length(x))})
    }
    else if(H0=="grouped"|H0==3){
      fun.group <- function(x, index){
        x.cons <- x
        for(jj in 1:max(index))
          x.cons[index==jj] <- rep(mean(x[index==jj], sum(index==jj)))
        x.cons
      }
      boot.constr <- apply(boot.rep, c(1,3), FUN=fun.group, index=group.index)
      boot.constr <- apply(boot.constr, c(1,3), FUN=function(x){x})
      coeff.constr[,] <- apply(res.fit$coeff.hat, 2, FUN=fun.group, index=group.index)
    }
    else if(H0=="nospatial"|H0==4){
      boot.constr <- boot.rep
      boot.constr[,,c("lambda0", "lambda2")] <- 0
      coeff.constr[,] <- res.fit$coeff.hat; coeff.constr[,c("lambda0", "lambda2")] <- 0
    }
    else if(H0=="noautoregressive"|H0==5){
      boot.constr <- boot.rep
      boot.constr[,,"lambda1"] <- 0
      coeff.constr[,] <- res.fit$coeff.hat; coeff.constr[,"lambda1"] <- 0
    }
    else if(H0=="noX"|H0==6){
      beta.names <- names(model.obj$beta_coeffs)[model.obj$beta_coeffs]
      boot.constr <- boot.rep
      boot.constr[,,beta.names] <- 0
      coeff.constr[,] <- res.fit$coeff.hat; coeff.constr[,beta.names] <- 0
    }
    else if(H0=="constrain"|H0==7){
      boot.constr <- boot.rep
      boot.constr[,,"lambda1"] <- -1*boot.rep[,,"lambda2"]
      boot.constr[,,"lambda2"] <- -1*boot.rep[,,"lambda1"]
      coeff.constr[,] <- 0
    }
    else return(error=paste("\n The null hypothesis", H0, "is not available..."))

    ### summary of bootstrap distribution
    statistica 	<- res.fit$coeff.hat - coeff.constr
    boot.rep		<- boot.rep-boot.constr
    fun.BOOT <- function(x){
      # statistiche di sintesi delle serie bootstrap
      nNAx <- sum(is.na(x))
      sdx <- sd(x, na.rm=T)
      mx <- mean(x, na.rm=T)
      if(nNAx==0)
        pvalue <- tseries::jarque.bera.test(x)$p.value
      else
        pvalue <- NA
      c(mx, sdx, pvalue, nNAx)
    }
    boot <- apply(boot.rep, c(2,3), fun.BOOT)
    dimnames(boot)[[1]] <- c("bias.boot", "sd.boot", "pvJBnormaltest.boot", "NAs.boot")
    devs.tsboot <- t(apply(devs.tsboot, c(2,3), mean))
    dimnames(devs.tsboot)[[2]] <- c("mean", "sd")
    ### test mediante il metodo del normal basic bootstrap
    # pvalue <- 2*pnorm(0, mean=abs(statistica/boot["sd.boot",,]), lower.tail = TRUE)
    pvalue <- 2*pnorm(abs(statistica/boot["sd.boot",,]), lower.tail = FALSE)
    ### risultati:
    dimnames(pvalue)[[1]] <- dimnames(devs.tsboot)[[1]] <- dimnames(boot)[[2]]
    dimnames(pvalue)[[2]] <- dimnames(boot)[[3]]
    boot["bias.boot",,] <- boot["bias.boot",,] - res.fit$coeff.hat
  }
  else return(error="n. repliche bootstrap insufficienti")

  #### restituzione risultati
  list(pvalue=pvalue, n.boot=n.boot, H0=H0, diagnostics.model=data$diagnostics, coeff.hat=res.fit$coeff.hat,
       diagnostics.coeff.boot=boot, diagnostics.sdevs.tsboot=devs.tsboot, warnings=data$warnings)
}

