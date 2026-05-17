fun.LBtest <- function(x, ...) {
  tryCatch(
    {
      if (length(x) > 10 && !all(is.na(x))) {
        test_result <- Box.test(x, lag = min(12, length(x) / 4), type = "Ljung-Box")
        return(test_result$p.value)
      } else {
        return(NA)
      }
    },
    error = function(e) {
      return(NA)
    }
  )
}

