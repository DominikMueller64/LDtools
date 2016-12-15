library('LDtools')
context("comb_wind")

set.seed(123L)

pos <- rnorm(100L)

out_exp <- structure(c(65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 
                       79, 80, 81, 82, 83, 84, 85, 85, 86, 87, 88, 88, 89, 89, 90, 90, 
                       91, 92, 93, 93, 93, 94, 94, 94, 95, 95, 95, 96, 96, 96, 96, 96, 
                       96, 96, 96, 96, 96, 97, 97, 97, 97, 97, 97, 97, 98, 98, 98, 98, 
                       98, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 100, 100, 100, 100, 
                       100, 100, 100, 100, 100, 100, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 
                       1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 3, 4, 
                       3, 4, 4, 4, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 8, 9, 10, 11, 
                       12, 13, 14, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 15, 
                       16, 17, 18, 19, 20, 21, 22, 23, 24, 15, 16, 17, 18, 19, 20, 21, 
                       22, 23, 24, 96, 97, 98, 99, 100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2), 
                     .Dim = c(87L, 3L), min_dist = c(-3,  4), max_dist = c(-2.7, 4.2))

test_that("'comb_wind' recognizes invalid arguments.", {
  expect_error(comb_wind(pos = pos, min_dist = min(pos), max_dist = max(pos)))
  expect_error(comb_wind(pos = 1, min_dist = 5, max_dist = -3))
  expect_error(comb_wind(pos = 1, min_dist = TRUE, max_dist = NA))
  expect_error(comb_wind(pos = 1, min_dist = c(8, 2), max_dist = c(1, NA)))
  expect_error(comb_wind(pos = 1, min_dist = 1, max_dist = c(2, 3)))
  expect_identical(out_exp,
      comb_wind(pos = sort(pos), min_dist = c(-3, 4), max_dist = c(-2.7,  4.2)))
})
