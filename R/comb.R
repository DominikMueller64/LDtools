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

#' @title Combinations of loci.
#' 
#' @description Create a matrix of combinations between loci for which linkage
#' disequilibrium (LD) should be computed.
#' 
#' @param pos An increasingly sorted numeric vector. The positions of the loci. 
#' @param pos2 An increasingly numeric vector. The positions of a second set of loci. 
#' 
#' @param min_dist A numeric vector. The minimum distances between loci. Each entry must be positive.
#' @param max_dist A numeric vector. The maximum distances between loci.
#' 
#' @param indices An integer vector. A vector with the indices of the loci.
#' @param indices2 An integer vector. A vector with the indices of a second set
#' of loci.
#' 
#' @param start A double. Specifies where to start with the sliding window.
#' @param width A double. Width of the sliding window.
#' @param advance A double. Increment of the sliding window.
#' 
#' @param k An integer. The number of nearest loci to be considered.
#' 
#' @return A numeric matrix.
#' 
#' @details
#' The purpose of these functions is to provide a convenient way to obtain all
#' desired combinations of loci for computation of LD between them. Here is an
#' overview of their functionality:
#' 
#' \itemize{
#'   \item comb_all: Combinations of all loci.
#'   \item comb_adj: Combinations of adjacent (neighboring) loci.
#'   \item comb_nearest: Combinations of a set of loci and its closest pendant
#'   in another set.
#'   \item comb_nearest_k: Combinations of a set of loci and its closest
#'   k pendant in another set.
#'   \item comb_flank: Combinations between flanking loci.
#'   \item comb_wind: Combinations of loci with a given minimum and maximum distance.
#'   \item comb_sliding: Combinations of loci within an advancing sliding window.
#' }
#' 
#' In general, these functions all return a matrix with three columns.
#' The first two columns refer to the indices of the pair of loci, while the
#' third column indicates to which group of combinations the pair belongs.
#' Currently, the assigned groups only vary for \code{\link{comb_sliding}},
#' where the group indicates to which window the pair of loci belongs.
#' 
#' The returned matrix is suitable as input for the function
#' \code{\link{LD_mult}} for efficient computation of LD.
#' 
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#'
#' @seealso \code{\link{LD_mult}}
#' @name comb
#' 
#' @examples
#' pos <- sort(runif(10))
NULL

#' @rdname comb
#' @description Combinations among loci with minimum and maximum distance.
#' @examples
#' comb_wind(pos, min_dist = c(0.1, 0.3, 0.8), max_dist = c(0.2, 0.7, 0.9))
#' @export
comb_wind <- function(pos, min_dist, max_dist) {
    
  if (length(min_dist) != length(max_dist))
    stop("'min_dist' and 'max_dist' must be have the same length.")
  
  if (!is.atomic(pos) || !is.numeric(pos) || is.unsorted(pos))
    stop("'pos' must be an increasingly sorted vector.")
  
  if (!is.numeric(min_dist) || !is.numeric(max_dist) || any(min_dist >= max_dist))
    stop(paste("'min_dist' and 'max_dist' must be numeric with",
               "min_dist < max_dist for all elements."))

  if (any(min_dist <= 0))
    stop("'min_dist' cannot contain entries smaller than or equal to zero.")
  
  matr <- purrr::map2(min_dist, max_dist, ~ .comb_wind(pos, .x, .y))
  for (i in seq_along(matr)) matr[[i]][, 3L] <- i
  matr <- do.call(what = rbind, args = matr)
  
  attr(x = matr, which = 'min_dist') <- min_dist
  attr(x = matr, which = 'max_dist') <- max_dist
  matr
}

#' @rdname comb
#' @description Combinations among adjacent (neighboring) loci.
#' @examples
#' comb_adj(pos)
#' @export
comb_adj <- function(pos) {
  if (!is.atomic(pos) || !is.numeric(pos) || is.unsorted(pos))
    stop("'pos' must be an increasingly sorted vector.")
  .comb_adj(pos)
}


#' @rdname comb
#' @description Combinations among all loci.
#' @examples
#' comb_all(pos)
#' @export
comb_all <- function(pos) {
  if (!is.atomic(pos) || !is.numeric(pos) || is.unsorted(pos))
    stop("'pos' must be an increasingly sorted vector.")
  .comb_all(pos)
}


#' @rdname comb
#' @description Combinations among nearest loci.
#' @examples
#'  pos2 <- sort(runif(10L, min = 0, max = max(pos)))
#' indices1 <- which(sort(c(pos, pos2)) %in% pos)
#' indices2 <- which(sort(c(pos, pos2)) %in% pos2)
#' comb_nearest(indices1, indices2, pos, pos2)
#' @export
comb_nearest <- function(indices, indices2, pos, pos2) {
  comb_nearest_k(indices, indices2, pos, pos2, k = 1L)
}

#' @rdname comb
#' @description Combinations among nearest loci.
#' @examples
#' comb_nearest_k(indices1, indices2, pos, pos2, 3)
#' @export
comb_nearest_k <- function(indices, indices2, pos, pos2, k) {
  
  if (!is.atomic(indices)) 
    stop("'indices1' must be a vector.")
  
  if (!is.atomic(indices2)) 
    stop("'indices2' must be a vector.")
  
  if (!is.atomic(pos) || !is.numeric(pos) || is.unsorted(pos))
    stop("'pos1' must be an increasingly sorted vector.")
  
  if (!is.atomic(pos2) || !is.numeric(pos2) || is.unsorted(pos2))
    stop("'pos2' must be an increasingly sorted vector.")
  
  if (length(indices) != length(pos) || length(indices2) != length(pos2))  
    stop(paste("'indices' and 'pos' as well as 'indices2' and 'pos2' must have",
               "the same length."))
  if (!is.numeric(k) || k < 1L || k > length(pos2))
    stop(paste("'k' must be a natural number and not exceed",
               "the number of loci in the second set."))
  
  .comb_nearest_k(indices, indices2, pos, pos2, k)
}

#' @rdname comb
#' @description Combinations between flaning loci.
#' @examples
#' comb_flank(indices1, indices2, pos, pos2)
#' @export
comb_flank <- function(indices, indices2, pos, pos2) {
  
  if (!is.atomic(indices)) 
    stop("'indices1' must be a vector.")
  
  if (!is.atomic(indices2)) 
    stop("'indices2' must be a vector.")
  
  if (!is.atomic(pos) || !is.numeric(pos) || is.unsorted(pos))
    stop("'pos1' must be an increasingly sorted vector.")
  
  if (!is.atomic(pos2) || !is.numeric(pos2) || is.unsorted(pos2))
    stop("'pos2' must be an increasingly sorted vector.")
  
  if (length(indices) != length(pos) || length(indices2) != length(pos2))  
    stop(paste("'indices' and 'pos' as well as 'indices2' and 'pos2' must have",
               "the same length."))
  
  .comb_flank(indices, indices2, pos, pos2)
}




#' @rdname comb
#' @description Combinations among loci in a sliding window.
#' @examples
#' comb_sliding(pos, 0, 0.2, 0.1)
#' @export
comb_sliding <- function(pos, start, width, advance) {
  if (!is.atomic(pos) || !is.numeric(pos) || is.unsorted(pos))
    stop("'pos' must be an increasingly sorted vector.")
  
  if (!is.numeric(start))
    stop("'start' must be numeric.")
  
  if (!is.numeric(width) || width <= 0 ||
      !is.numeric(advance) || advance <= 0)
    stop("'width' and 'advance' must be numeric and > 0.")
  
  matr <- .comb_sliding(pos, start, width, advance)
  eps <- sqrt(.Machine$double.eps)
  n_wind <- ceiling((pos[length(pos)] - start - width - eps) / advance)
  begin <- start + advance * seq(from = 0, to = n_wind - 1)
  attr(x = matr, which = 'begin') <- begin
  attr(x = matr, which = 'end') <- begin + width
  matr
}
