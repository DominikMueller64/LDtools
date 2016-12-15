library('LDtools')
context("LD")

# set.seed(123L)

x <- genetics::genotype( c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
                           'C/A', 'C/C', NA_character_, 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
                           'C/A', 'C/A', 'C/A', 'A/A') )

xt <- drop(transform_geno(matrix(x), sep = '/')$geno)

y <- genetics::genotype( c('T/A', 'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A',
                           'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', NA_character_, 'T/T',
                           'T/A', 'T/A', 'T/A', 'T/T') )

yt <- drop(transform_geno(matrix(y), sep = '/')$geno)

microbenchmark::microbenchmark(
ret_gen <- genetics::LD(g1 = x, g2 = y),
ret_LDt <- LD(x = xt, y = yt, is_phased = FALSE, check = TRUE)
)

test_that("Output of LDtools::LD corresponds to genetics::LD", {
  expect_equal(ret_LDt$D, ret_gen$D, tolerance = 1e-3)
  expect_equal(ret_LDt$Dprime, ret_gen[["D'"]], tolerance = 1e-3)
  expect_equal(ret_LDt$r, ret_gen$r, tolerance = 1e-3)
  expect_equal(ret_LDt$chi2, ret_gen[["X^2"]], tolerance = 1e-3)
  expect_equal(ret_LDt$pval, ret_gen[["P-value"]], tolerance = 1e-3)
})
