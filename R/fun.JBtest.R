fun.JBtest <- function(x, ...) {
  tryCatch(
    {
      if (length(x) > 7 && !all(is.na(x))) {
        # Simple Jarque-Bera test implementation
        n <- length(x)
        x_centered <- x - mean(x, na.rm = TRUE)
        s <- sqrt(sum(x_centered^2, na.rm = TRUE) / n)
        skewness_val <- sum(x_centered^3, na.rm = TRUE) / (n * s^3)
        kurtosis_val <- sum(x_centered^4, na.rm = TRUE) / (n * s^4)
        jb_stat <- n * (skewness_val^2 / 6 + (kurtosis_val - 3)^2 / 24)
        p_value <- 1 - pchisq(jb_stat, df = 2)
        return(p_value)
      } else {
        return(NA)
      }
    },
    error = function(e) {
      return(NA)
    }
  )
}
