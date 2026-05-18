#' Compute a Jarque-Bera Normality-Test P-Value
#'
#' Computes the p-value of the Jarque-Bera normality test for a numeric vector.
#'
#' The function implements a simple Jarque-Bera test directly from sample
#' skewness and kurtosis. It returns `NA_real_` when the input is too short,
#' entirely missing, has zero variance, or produces an error.
#'
#' @param x Numeric vector. Input values to test for normality.
#' @param ... Additional arguments. Currently unused.
#'
#' @return Numeric scalar. The Jarque-Bera test p-value, or `NA_real_` if the
#'   test cannot be computed.
#'
#' @details
#' The test statistic is computed as:
#' \deqn{
#'   JB = n \left(\frac{S^2}{6} + \frac{(K - 3)^2}{24}\right)
#' }
#' where `S` is sample skewness and `K` is sample kurtosis. Under the null
#' hypothesis of normality, the statistic is asymptotically distributed as a
#' chi-squared random variable with 2 degrees of freedom.
#'
#' @examples
#' fun_jb_test(rnorm(100))
#'
#' fun_jb_test(c(1, 2, 3, NA, 5, 6, 7, 8))
#'
#' @seealso [stats::pchisq()]
#'
#' @export
fun_jb_test <- function(x, ...) {
  tryCatch(
    {
      if (length(x) > 7 && !all(is.na(x))) {
        x <- x[!is.na(x)]

        n <- length(x)

        if (n <= 7) {
          return(NA_real_)
        }

        x_centered <- x - mean(x)

        standard_deviation <- sqrt(sum(x_centered^2) / n)

        if (is.na(standard_deviation) || standard_deviation == 0) {
          return(NA_real_)
        }

        skewness_value <- sum(x_centered^3) /
          (n * standard_deviation^3)

        kurtosis_value <- sum(x_centered^4) /
          (n * standard_deviation^4)

        jb_statistic <- n * (
          skewness_value^2 / 6 +
            (kurtosis_value - 3)^2 / 24
        )

        p_value <- 1 - stats::pchisq(jb_statistic, df = 2)

        return(p_value)
      }

      NA_real_
    },
    error = function(error) {
      NA_real_
    }
  )
}

