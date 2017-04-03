# LDtools: Tools for Computation of Linkage Disequilibrium.
#
# Copyright (C) 2016 Dominik Mueller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
#' testing purposes. The allele frequency of each locus is drawn from a beta with
#' parameters \code{shape1} and \code{shape2}.
#' Note: If \code{shape1 == shape2}, the beta distribution is symmetric.
#'
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#'
#' @examples
#' .gen_geno(n = 2, m = 5)
#'
#' @export
.gen_geno <- function(nucleotides = c("A", "C", "G", "T"),
                     n, m, prob_na = 0.2,
                     shape1 = 4, shape2 = 4, sep = '/') {

  if (!is.atomic(nucleotides) || !is.character(nucleotides))
    stop("'nucleotides' must be a character vector")

  if (!is.character(sep) | length(sep) != 1)
    stop("'sep' must be a single character.")

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
