.new_sdpd_data <- function(series,
                           px,
                           lat = NULL,
                           lon = NULL,
                           rr_xx = NULL,
                           rr_groups = NULL) {
  structure(
    list(
      series = series,
      px = px,
      lat = lat,
      lon = lon,
      rr_xx = rr_xx,
      rr_groups = rr_groups
    ),
    class = "sdpd_data"
  )
}
