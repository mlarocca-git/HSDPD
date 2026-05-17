
plot.sdpd.series.FITs <- function(df.results, latitude, longitude, n.digits=1, main=NULL, sub=NULL, xlab=NULL, ylab=NULL){
  tempo <- dimnames(df.results$fitted)[[2]]
  if(is.Date(tempo))
    range.tempo <- as.Date(tempo)
  else if(is.numeric(as.numeric(tempo)))
    range.tempo <- as.Date(as.numeric(tempo))
  else
    range.tempo <- tempo
  if(is.null(sub))
    sub <- paste("\n(based on data from ", range(tempo)[1], " to ", range(tempo)[2], ")", sep="")
  if(is.null(main))
    main <- paste("Estimated model series (longitude=", longitude, " - latitude=", latitude, ")", sep="")
  serie1 <- df.results |>
    mutate(Latitude=round(lat,n.digits)) |>
    mutate(Longitude=round(lon,n.digits)) |>
    filter(Latitude==round(latitude,n.digits),Longitude==round(longitude,n.digits)) |>
    select(fitted, resid)
  if(dim(serie1)[1]==0){
    return("These coordinates are not present in the database")
  }
  offset.axis <- range(serie1$fitted, na.rm=T)[1]-0.1*diff(range(serie1$fitted, na.rm=T))
  res <- data.frame(time=range.tempo, fitted=as.numeric(serie1$fitted), resid=as.numeric(serie1$resid)) %>%
    mutate(observed=resid+fitted) %>%
    ggplot(aes(x=time)) +
    geom_line(aes(y=fitted, color="Fitted")) +
    geom_line(aes(y=observed, color="Observed")) +
    geom_line(aes(y=resid+offset.axis, color="Residuals")) +
    geom_hline(aes(yintercept = mean(resid, na.rm=T)+offset.axis)) +
    guides(x = guide_axis(angle = 0)) +
    scale_x_date(date_labels = "%Y", date_breaks="1 year") +
    theme(legend.position = "bottom") +
    scale_y_continuous(
      # Features of the first axis
      name = ylab,
      # Add a second axis and specify its features
      sec.axis = sec_axis(transform=~.-offset.axis, name="2° axis for residuals")) +
    labs(title = main, x=xlab, y=ylab, color="")
  res
}
