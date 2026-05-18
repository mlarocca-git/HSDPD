#' Compute a Ljung-Box Test P-Value
#'
#' Computes the p-value of the Ljung-Box test for serial autocorrelation in a
#' numeric vector.
#'
#' The function is a safe wrapper around [stats::Box.test()]. It returns
#' `NA_real_` when the input is too short, entirely missing, has zero variance,
#' or produces an error.
#'
#' @param x Numeric vector. Input values to test for serial autocorrelation.
#' @param ... Additional arguments. Currently unused.
#'
#' @return Numeric scalar. The Ljung-Box test p-value, or `NA_real_` if the test
#'   cannot be computed.
#'
#' @details
#' The test is computed using `type = "Ljung-Box"` and lag equal to
#' `min(12, floor(length(x) / 4))`.
#'
#' Missing values are removed before computing the test.
#'
#' @examples
#' fun_lb_test(rnorm(100))
#'
#' fun_lb_test(c(1, 2, 3, NA, 5, 6, 7, 8, 9, 10, 11))
#'
#' @seealso [stats::Box.test()]
#'
#' @export
fun_lb_test <- function(x, ...) {
  tryCatch(
    {
      if (length(x) > 10 && !all(is.na(x))) {
        x <- x[!is.na(x)]

        if (length(x) <= 10) {
          return(NA_real_)
        }

        if (stats::sd(x) == 0) {
          return(NA_real_)
        }

        lag_value <- min(12, floor(length(x) / 4))

        test_result <- stats::Box.test(
          x,
          lag = lag_value,
          type = "Ljung-Box"
        )

        return(test_result$p.value)
      }

      NA_real_
    },
    error = function(error) {
      NA_real_
    }
  )
}
