build.sdpd.series <- function(df.obj=NULL, px=NULL, rry, rrxx=NULL, lon=NULL, lat=NULL,
                              rrgroups=NULL, label_groups=NULL,
                              model,
                              type.w=c("distance", "correlations")[1],
                              ww.index=NULL, ww.values=NULL,
                              px.neighbors=NULL,
                              vec.options=NULL){
  if(is.null(px))
    px <- df.obj$px
  if(is.data.frame(rry) | is.matrix(rry))
    res <- read.sdpd.series.from.dataframe(px=px, rry=rry, rrXX=rrxx, model=model, longit=lon, latit=lat,
                rrgroups=rrgroups, ww.index=ww.index, ww.values=ww.values, px.neighbors=px.neighbors)
  else if(inherits(rry, "SpatRaster")){
    if(is.null(vec.options)){
      ## ATTENZIONE (da fare): qui bisognerà testare la consistenza di vec.options
      vec.options <- list(px.core=1, px.neighbors=8, t_frequency=1, na.rm=T, NAcovs="pairwise.complete.obs")
    }
    res <- read.sdpd.series.from.raster(px=df.obj$px, rry=rry, rrXX=rrxx, lon=lon, lat=lat,
                rrgroups=rrgroups, model=model, label_groups=label_groups, type.w=type.w, vec.options=vec.options)
  }
  else return(list(error="L'oggetto rry non è del formato atteso (dataframe, matrix oppure raster)"))
  res
}

