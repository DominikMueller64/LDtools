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

#' Compute linkage phase similarity
#'
#' @description Compute the linkage phase similary between two populations.
#'
#' @param X1 A numeric matrix. Data for the first population.
#' @param X2 A numeric matrix. Data for the second population.
#' @param min_dist A numeric vector. Minimum distances.
#' @param max_dist A numeric vector. Maximum distances.
#' @param method A string. Method use for computing linkag phase
#' similarities, one of 'correlation', 'cosine_similarity' or
#' 'sign_match'.
#' @param signif_level A double. If the chi-square test yields a p-value above
#' \code{signif_level} for any population, this pair of loci will be
#' excluded from the computation of linkage phases.
#'
#' @inheritParams LD_mult
#'
#' @return A data.frame with columns min_dist, max_dist, n (number of pairs),
#' and lps (linkage phase similarity coefficient).
#'
#' @details
#' The computation is faster if \code{signif_level = 1} and all pairs are included.
#'
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#'
#' @examples
#' data('population', package = 'LDtools')
#' pos <- map$pos
#' mid <- nrow(X) %/% 2L
#' X1 <- X[1L:mid, , drop = FALSE]
#' X2 <- X[(mid + 1L):nrow(X), , drop = FALSE]
#' min_dist <- seq(sqrt(.Machine$double.eps), floor(max(pos) / 2), by = 1)
#' max_dist <- min_dist + 1
#'
#' dat <- LP(X1, X2, pos, min_dist, max_dist, method = 'correlation')
#' with(dat,
#'  plot(x = (min_dist + max_dist) / 2, y = lps, type = 'o', ylim = c(0, 1))
#' )
#'
#' @export
LP <- function(X1, X2, pos, min_dist, max_dist,
               method = c('correlation', 'cosine_similarity', 'sign_match'),
               signif_level = 1,
               is_phased = TRUE, any_na = FALSE, check = TRUE) {


  method <- match.arg(arg = method)

  if (!is.numeric(signif_level) || length(signif_level) != 1L || signif_level <= 0)
    stop("'signif_level' must be numeric, positive and of length 1.")

  matr <- comb_wind(pos = pos, min_dist = min_dist, max_dist = max_dist)

  r_only <- TRUE
  if (signif_level < 1) r_only <- FALSE

  LD1 <- LD_mult(X = X1, matr = matr, pos = NULL, is_phased = is_phased,
                 any_na = any_na, r_only = r_only, check = check)

  LD2 <- LD_mult(X = X2, matr = matr, pos = NULL , is_phased = is_phased,
                 any_na = any_na, r_only = r_only, check = check)

  r1 <- LD1$r
  r2 <- LD2$r
  blocks <- LD1$block
  if (signif_level < 1) {
    ix <- (LD1$pval <= signif_level) & (LD2$pval <= signif_level)
    r1 <- r1[ix]
    r2 <- r2[ix]
    blocks <- blocks[ix]
  }
  rm('LD1', 'LD2')
  ix <- is.finite(r1) & is.finite(r2)
  r1 <- r1[ix]
  r2 <- r2[ix]
  blocks <- blocks[ix]

  lps <- numeric(length(min_dist))
  n <- integer(length(min_dist))

  for (block in seq_along(min_dist)) {
    ix <- blocks == block
    r1s <- r1[ix]
    r2s <- r2[ix]
    n[block] <- length(r1s)
    if (method == 'correlation') {
      lps[block] <- suppressWarnings(stats::cor(r1s, r2s, method = 'pearson'))
    } else if (method == 'cosine_similarity') {
      lps[block] <- drop(crossprod(r1s, r2s) / sqrt(crossprod(r1s) * crossprod(r2s)))
    } else if (method == 'sign_match') {
      lps[block] <- mean(sign(r1s) == sign(r2s))
    } else {
      stop('No matching method.')
    }
  }

  data.frame('min_dist' = min_dist,
             'max_dist' = max_dist,
             'n' = n, 'lps' = lps)
}
