parallelize.sdpd.data <- function(obj.series)
{
  ## group pixels by groups
  ## series, X, px, group, lon, lat, model, ww.index, ww.values, px.neighbors
  px <- obj.series$px
  df.gruppi   <- tibble(gruppo=obj.series$group[as.character(px),], px=as.numeric(px)) %>%
    nest_by(.by=gruppo, .key="gruppo")
  
  ## For each group of pixels, prepare data and quantities for H-SDPD model estimation
  df.data <- df.gruppi$gruppo %>% 
              map(build.sdpd.series, rry=obj.series$series, rrxx=obj.series$X, lon=obj.series$lon, lat=obj.series$lat,
                  rrgroups=obj.series$group, model=obj.series$model, 
              ww.index=obj.series$ww.index, ww.values=obj.series$ww.values, px.neighbors=obj.series$px.neighbors)
  
  ## returning list of objects
  df.data
}
