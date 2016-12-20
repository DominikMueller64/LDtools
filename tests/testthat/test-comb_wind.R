library('LDtools')
context("comb_wind")

set.seed(123L)

pos <- rnorm(100L)


test_that("'comb_wind' recognizes invalid arguments.", {
  expect_error(comb_wind(pos = 1, min_dist = 5, max_dist = -3))
  expect_error(comb_wind(pos = 1, min_dist = TRUE, max_dist = NA))
  expect_error(comb_wind(pos = 1, min_dist = c(8, 2), max_dist = c(1, NA)))
  expect_error(comb_wind(pos = 1, min_dist = 1, max_dist = c(2, 3)))
})
