library('LDtools')
context("comb")

## set.seed(123L)
pos <- c(-0.560475646552213, -0.23017748948328, 0.070508391424576, 0.129287735160946,
         0.460916205989202, 1.55870831414912, 1.71506498688328)
pos2 <- c(-1.26506123460653, -0.686852851893526, -0.445661970099958,
          0.359813827057364, 1.22408179743946)
indices <- 1L:7L
indices2 <- 10L:14L

comb_wind_ref <- structure(c(1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 7, 6, 7, 6, 7, 6,
7, 1, 1, 1, 1, 1, 1, 1, 1, 1), .Dim = c(9L, 3L), min_dist = 1, max_dist = 2)
comb_adj_ref <- structure(c(1, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1,
1), .Dim = c(6L, 3L))
comb_all_ref <- structure(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4,
4, 5, 5, 6, 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6,
7, 6, 7, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1), .Dim = c(21L, 3L))
comb_flank_ref <- structure(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 11, 12, 12, 13,
12, 13, 12, 13, 13, 14, 14, 14, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
6, 7), .Dim = c(12L, 3L))
comb_nearest_k_ref <- structure(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 12, 11,
12, 11, 13, 12, 13, 12, 13, 14, 14, 13, 14, 13, 1, 1, 2, 2, 3,
3, 4, 4, 5, 5, 6, 6, 7, 7), .Dim = c(14L, 3L))
comb_sliding_ref <- structure(c(3, 3, 4, 3, 3, 4, 4, 4, 5, 5, 4, 5, 5, 5, 1, 1, 1,
2, 2, 2, 3), .Dim = c(7L, 3L), begin = c(0, 0.05, 0.1, 0.15,
0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7), end = c(1,
1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55,
1.6, 1.65, 1.7))
comb_wind_sets_ref <- structure(c(1, 1, 2, 3, 4, 5, 1, 2, 2, 3, 4, 5, 12, 13, 13, 13, 
13, 14, 13, 13, 14, 14, 14, 14, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 
2, 2), .Dim = c(12L, 3L), min_dist = c(0.1, 0.5), max_dist = c(1, 
1.5))
comb_all_sets_ref <- structure(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 
4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 10, 11, 
12, 13, 14, 10, 11, 12, 13, 14, 10, 11, 12, 13, 14, 10, 11, 12, 
13, 14, 10, 11, 12, 13, 14, 10, 11, 12, 13, 14, 10, 11, 12, 13, 
14, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), .Dim = c(35L, 
3L))



test_that("'checks'", {
  expect_error(LDtools:::check_pos(1L:2L))
  expect_true(LDtools:::check_pos(c(-0.1, 2)))
  expect_error(LDtools:::check_pos(c(0.11, 0.1)))
  expect_error(LDtools:::check_pos(list(0.1, 0.2)))

  expect_error(LDtools:::check_min_max_dist(5, -3))
  expect_error(LDtools:::check_min_max_dist(TRUE, NA))
  expect_error(LDtools:::check_min_max_dist(c(8, 2), c(1, NA)))
  expect_error(LDtools:::check_min_max_dist(1, c(2, 3)))
})

test_that("'comb'", {
  expect_equal(comb_wind(pos = pos, min_dist = 1, max_dist = 2), comb_wind_ref)
  expect_equal(comb_adj(pos), comb_adj_ref)
  expect_equal(comb_all(pos), comb_all_ref)
  expect_equal(comb_flank(indices, indices2, pos, pos2), comb_flank_ref)
  expect_equal(comb_nearest_k(indices, indices2, pos, pos2, k = 2L), comb_nearest_k_ref)
  expect_equal(comb_sliding(pos, start = 0.0, width = 1.0, advance = 0.05), comb_sliding_ref)
  expect_equal(comb_wind_sets(indices, indices2, pos, pos2, min_dist = c(0.1, 0.5),
                              max_dist = c(1, 1.5)), comb_wind_sets_ref)
  expect_equal(comb_all_sets(indices, indices2, pos, pos2), comb_all_sets_ref)
})
