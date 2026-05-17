plot.sdpd.estimates <- function(obj.results, item="coeff.hat", main=NULL, sub=NULL, limits=NULL, mid=0){
  res <- list()
  item.bis <- as.matrix(obj.results[[item]])
  nomi <- dimnames(item.bis)[[2]]
  if(length(nomi)==0 & dim(item.bis)[2]==1){
    nomi <- as.character(item)
  }
  if(length(mid)==1)
    mid <- rep(mid, length(nomi))
  for(ii in 1:length(nomi)){
    if(dim(item.bis)[2]==1)
      valori <- as.numeric(item.bis)
    else
      valori <- item.bis[,ii]
    dati <- data.frame(lon=obj.results$lon, lat=obj.results$lat, group=obj.results$group$LABEL, value=valori)
    res[[ii]] <- dati %>% ggplot(aes(x=lon, y=lat, group=group, colour = value)) +
      geom_point(size=size.point) +
      scale_colour_gradient2(high = "red", low = "blue", mid = "white", midpoint = mid[ii], limits=limits) +
      labs(x="", y="", title = paste(main, nomi[ii], sub, sep=""), colour="value")
  }
  names(res) <- nomi
  res
}

