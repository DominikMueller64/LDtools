library('LDtools')
context("LD")

# set.seed(123L)

LDr <- function(x, y) {

  geno_x <- transform_geno(x, sep = '/')
  geno_y <- transform_geno(y, sep = '/')

  x <- geno_x$geno
  pA <- geno_x$p
  y <- geno_y$geno
  pB <- geno_y$p

  pa <- 1.0 - pA
  pb <- 1.0 - pB

  Dmin <- max(-pA * pB, -pa * pb)
  pmin <- pA * pB + Dmin

  Dmax <- min(pA * pb, pB * pa)
  pmax <- pA * pB + Dmax

  counts <- table(x, y)

  n3x3 <- matrix(0, nrow=3, ncol=3)
  colnames(n3x3) <- rownames(n3x3) <- 0:2

  # ensure the matrix is 3x3, with highest frequency values in upper left
  for(i in rownames(counts))
    for(j in colnames(counts))
      n3x3[3-as.numeric(i),3-as.numeric(j)] <- counts[i,j]


  loglik <- function(pAB,...)
  {
    (2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
      (2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
      (2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
      (2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) +
      n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
  }

  solution <- optimize(
    loglik,
    lower=pmin+.Machine$double.eps,
    upper=pmax-.Machine$double.eps,
    maximum=TRUE
  )
  pAB <- solution$maximum

  estD <- pAB - pA*pB
  if (estD>0)
    estDp <- estD / Dmax
  else
    estDp <- estD / Dmin

  n <-  sum(n3x3)

  corr <- estD / sqrt( pA * pB * pa * pb )

  dchi <- (2*n*estD^2)/(pA * pa * pB* pb)
  dpval <- 1 - pchisq(dchi,1)

  retval <- list(
    call=match.call(),
    "D"=estD,
    "D'"=estDp,
    "r" = corr,
    "R^2" = corr^2,
    "n"=n,
    "X^2"=dchi,
    "P-value"=dpval
  )

  return(retval)
}

x <- matrix( c('C/A', 'C/A', 'C/C', 'C/A', 'C/C', 'C/A', 'C/A', 'C/A',
                           'C/A', 'C/C', NA_character_, 'A/A', 'C/A', 'A/A', 'C/A', 'C/C',
                           'C/A', 'C/A', 'C/A', 'A/A') )

xt <- drop(transform_geno(matrix(x), sep = '/')$geno)

y <- matrix( c('T/A', 'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', 'T/A',
                           'T/A', 'T/T', 'T/A', 'T/T', 'T/A', 'T/A', NA_character_, 'T/T',
                           'T/A', 'T/A', 'T/A', 'T/T') )

yt <- drop(transform_geno(matrix(y), sep = '/')$geno)


ret_LDr <- LDr(x = x, y = y)
ret_LDt <- LD(x = xt, y = yt, is_phased = FALSE, check = TRUE)

test_that("Output of LDtools::LD corresponds to LDretics::LD", {
  expect_equal(ret_LDt$D, ret_LDr$D, tolerance = 1e-3)
  expect_equal(ret_LDt$Dprime, ret_LDr[["D'"]], tolerance = 1e-3)
  expect_equal(ret_LDt$r, ret_LDr$r, tolerance = 1e-3)
  expect_equal(ret_LDt$chi2, ret_LDr[["X^2"]], tolerance = 1e-3)
  expect_equal(ret_LDt$pval, ret_LDr[["P-value"]], tolerance = 1e-3)
})
