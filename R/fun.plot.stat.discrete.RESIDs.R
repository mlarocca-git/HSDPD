fun.plot.stat.discrete.RESIDs <- function(
    df.results,
    statistic = mean,
    main=NULL,
    significant.test = FALSE,
    BYadjusted = FALSE,
    alpha=NULL,
    mid_value=0,...) {
  # Compute summary statistics on residuals
  if (is.data.frame(df.results$resid)) {
    dati <- data.frame(
      lon = df.results$lon,
      lat = df.results$lat,
      value = apply(df.results$resid[, -1], 1, FUN = statistic, ...)
    )
  } else if (is.list(df.results$resid)) {
    dati <- data.frame(
      lon = df.results$lon,
      lat = df.results$lat,
      value = unlist(lapply(df.results$resid, FUN = statistic, ...))
    )
  } else {
    # Handle simple vector case
    dati <- data.frame(
      lon = df.results$lon,
      lat = df.results$lat,
      value = statistic(df.results$resid, ...)
    )
  }
  
  # Adjust p-values if requested
  if (BYadjusted) {
    dati$value <- p.adjust(dati$value, method = "BY")
    etichetta <- "\nGlobal\n"
  }
  
  # Convert to factor for significance test
  if (significant.test) {
    dati$value <- as.factor(ifelse(dati$value < alpha, "Significant", "Not significant"))
    # Create discrete plot
    res <- dati %>%
      ggplot(aes(x = lon, y = lat, colour = value)) +
      geom_point(size=size.point) +
      guides(fill = "none") +
      labs(
        title = main,
        x = "Longitude",
        y = "Latitude",
        colour = ""
      )
  }
  else{
    # Create plot
    res <- dati %>%
      ggplot(aes(x = lon, y = lat, colour = value)) +
      geom_point(size=size.point) +
      guides(fill = "none") +
      scale_colour_gradient2(high = "red", low = "blue", mid = "white", midpoint = mid_value) +
      labs(
        title = main,
        x = "Longitude",
        y = "Latitude",
        colour = ""
      )
    
  }
  res
}


