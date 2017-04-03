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

#' Transform genotypes
#'
#' @description Transform genotypes from two-letter coding to numeric
#' coding.
#'
#' @param geno A character matrix containing the genotypes.
#' @param sep A character string use to split the genotypes in \code{geno}.
#' @param check Should checks be performed?
#'
#' @return A numeric matrix of the same dimensions as geno containing the
#' recoded genotypes.
#'
#' @details A genotype must consist of a sequence of three characters,
#' where the mittle character is the seperator as specified in
#' \code{sep}. Recoding is performed based on the major allele at the
#' respective locus, i.e., if for instance the major allele is "A",
#' the recoded genotypes of "A/A", "A/B", "B/A" and "B/B" will be,
#' 2, 1, 1, 0, given \code{sep = "/"}.
#'
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#'
#' @examples
#' geno <- structure(c("A/G", "G/A", "A/A", NA, "A/A", "G/G", NA, "G/G",
#'                     "G/G", "A/A", "G/G", NA, "C/C", "C/C", "C/A", NA, "C/C", "C/C"),
#'                  .Dim = c(6L, 3L))
#' transform_geno(geno, sep = '/')
#'
#' @export
transform_geno <- function(geno, sep = '/', check = TRUE){

  if (check) {
    if (!is.matrix(geno) || !is.character(geno))
      stop("'geno' must be a character matrix.")

    if (!is.character(sep) || length(sep) != 1L)
      stop("'sep' must be a string.")

    # Check if the structure of genotypes is correct.
    for(x in split(x = geno, f = col(geno))) {
      x <- x[!is.na(x)]
      uq <- character()
      for(y in strsplit(x = x, split = '', fixed = TRUE)) {
        if(length(y) != 3L || y[2L] != sep)
          stop(paste0('Genotypic data do not comply with the',
                      'required format'))
        uq <- unique(c(uq, y))
      }
      if (!(length(uq) %in% c(2L, 3L)))
        stop(paste0('There must be at least one and at maximum',
                    'two allels per locus.'))
    }
  }

  n <- ncol(geno)
  new_geno <- matrix(NA_integer_, ncol = n, nrow = nrow(geno),
                     dimnames = dimnames(geno))
  major_freq <- numeric(n)
  for (i in seq_len(n)) {
    l <- strsplit(x = geno[, i, drop = TRUE], split = sep, fixed = TRUE)
    tb <- table(unlist(l, recursive = FALSE, use.names = FALSE))
    major_allele <- names(which.max(tb))
    major_freq[i] <- tb[major_allele] / sum(tb)
    new_geno[, i] <- purrr::map_int(l, ~ sum(.x == major_allele))
  }

  return(list('geno' = new_geno, 'p' = major_freq))
}
