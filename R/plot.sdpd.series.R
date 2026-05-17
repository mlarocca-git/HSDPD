plot.sdpd.series <- function(rry, lat, lon, n.digits=15, main=NULL, sub=NULL, xlab=NULL, ylab=NULL){
  serie1 <- rry |>
    mutate(Latitude=round(latitude,n.digits)) |>
    mutate(Longitude=round(longitude,n.digits)) |>
    filter(Latitude==round(lat,n.digits),Longitude==round(lon,n.digits)) |>
    select(-c(Longitude, Latitude, longitude, latitude))
  tempo <- dimnames(serie1)[[2]]
  range.tempo <- as.Date(as.numeric(tempo))
  if(dim(serie1)[1]==0){
    return("These coordinates are not present in the database")
  }
  res <- data.frame(time=range.tempo, series=as.numeric(serie1)) %>%
    ggplot(aes(x=time)) +
    geom_line(aes(y=series, color="series")) +
    guides(x = guide_axis(angle = 0)) +
    scale_x_date(date_labels = "%Y", date_breaks="1 year") +
    theme(legend.position = "bottom") +
    scale_y_continuous(
      # Features of the first axis
      name = ylab) +
    labs(title = main, x=xlab, y=ylab, color="")
  res
}

