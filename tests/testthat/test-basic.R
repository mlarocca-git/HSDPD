test_that("package loads", {

  expect_true(requireNamespace("HSDPD", quietly = TRUE))

})

test_that("Jarque-Bera helper returns a valid p-value or NA", {
  x <- rnorm(100)
  p_value <- fun_jb_test(x)

  expect_true(is.numeric(p_value))
  expect_length(p_value, 1)
  expect_true(is.na(p_value) || (p_value >= 0 && p_value <= 1))
})

test_that("Ljung-Box helper returns a valid p-value or NA", {
  x <- rnorm(100)
  p_value <- fun_lb_test(x)

  expect_true(is.numeric(p_value))
  expect_length(p_value, 1)
  expect_true(is.na(p_value) || (p_value >= 0 && p_value <= 1))
})

