#' Generate unphased genotypes
#' 
#' @description Generate a random sample of unphased genotypes.
#' 
#' @param nucleotides A character vector. The allowed nucleatides.
#' @param n An integer. The number of genotypes.
#' @param m An integer. The number of loci.
#' @param prob_na A numeric. The probability of a missing value.
#' @param shape1 A numeric. The first shape parameter of the beta distribution.
#' @param shape2 A numeric. The second shape parameter of the beta distribution.
#' @param sep A character. The seperator used between nucleotides.
#' 
#' @details This function can be used to simulate genotypes and is mainly for
#' testing purposes. 
#' #' Note: If \code{shape1 == shape2}, the beta distribution is symmetric.
#' 
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#' 
#' @examples
#' NULL
#'  
#' @export
gen_geno <- function(nucleotides = c("A", "C", "G", "T"),
                     n, m, prob_na = 0.2,
                     shape1 = 4, shape2 = 4, sep = '/') {  
  geno <- matrix(NA_character_, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    nuc <- sample(x = nucleotides, size = 2L, replace = FALSE)
    p <- stats::rbeta(n = 1L, shape1 = shape1, shape2 = shape1)
    a <- sample(x = nuc, size = n, replace = TRUE, prob = c(p, 1.0-p))
    b <- sample(x = nuc, size = n, replace = TRUE, prob = c(p, 1.0-p))
    geno[, j] <- paste(a, b, sep = sep)
  }
  geno[stats::runif(n = m * n) < prob_na] <- NA_character_
  geno
}