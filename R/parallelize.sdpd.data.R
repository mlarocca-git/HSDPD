parallelize.sdpd.data <- function(obj.series, model)
{
  ## group pixels and data in obj.series by groups
  n.groups <- length(table(obj.series$group[,"COD"]))
  if(n.groups>1){
    px <- obj.series$px
    df.gruppi   <- tibble(gruppo=obj.series$group[as.character(px),], px=as.numeric(px)) %>%
      nest_by(.by=gruppo, .key="gruppo")
    
    ## For each group of pixels, prepare data and quantities for H-SDPD model estimation
    df.data <- df.gruppi$gruppo %>% 
      map(build.sdpd.series, rry=obj.series$series, rrxx=obj.series$X, lon=obj.series$lon, lat=obj.series$lat,
          rrgroups=obj.series$group, model=model, check=FALSE,
          ww.index=obj.series$ww.index, ww.values=obj.series$ww.values, px.neighbors=obj.series$px.neighbors)
  }
  else{
    cat("\nWarning: no parallelization has been made since there is only one group")
    df.data <- obj.series
  } 
  
  ## returning list of objects
  df.data
}


