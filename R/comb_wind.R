##' @title Combinations among loci with minimum and maximum distance.
##' @description Combinations among loci with minimum and maximum distance.
##' @param pos A numeric vector. Increasingly sorted, contains the positions of loci.
##' @param min_dist A double. Minimum distance.
##' @param max_dist A double. Maximum distance.
##' @return A numeric matrix.
##' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
##' @examples
##' comb_wind(pos, min_dist = c(0.1, 0.3, 0.8), max_dist = c(0.2, 0.7, 0.9))
##' @export
comb_wind <- function(pos, min_dist, max_dist) {

  if (length(min_dist) != length(max_dist))
    stop("'min_dist' and 'max_dist' must be have the same length.")

  if (!is.atomic(pos) || !is.numeric(pos) || is.unsorted(pos))
    stop("'pos' must be an increasingly sorted vector.")

  if (!is.numeric(min_dist) || !is.numeric(max_dist) || any(min_dist >= max_dist))
    stop(paste("'min_dist' and 'max_dist' must be numeric with",
               "min_dist < max_dist for all elements."))

  matr <- purrr::map2(min_dist, max_dist, ~ .comb_wind(pos, .x, .y))
  for (i in seq_along(matr)) matr[[i]][, 3L] <- i
  matr <- do.call(what = rbind, args = matr)

  attr(x = matr, which = 'min_dist') <- min_dist
  attr(x = matr, which = 'max_dist') <- max_dist
  matr
}

##' @title Combinations among different sets of loci with minimum and maximum distance.
##' @description Combinations among different sets of loci with minimum and maximum distance.
##' @param indices_a A integer vector. Indices of the first set of loci.
##' @param indices_b A integer vector. Indices of the second set of loci.
##' @param pos_a A numeric vector. Increasingly sorted, contains the positions of the first set of loci.
##' @param pos_b A numeric vector. Increasingly sorted, contains the positions of the second set of loci.
##' @param min_dist A double. Minimum distance.
##' @param max_dist A double. Maximum distance.
##' @return A numeric matrix.
##' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
##' @examples
##' NULL
##' @export
comb_wind_sets <- function(indices_a, indices_b,
                           pos_a, pos_b,
                           min_dist, max_dist) {

  matr <- vector(mode = 'list', length = length(pos_a))
  for(i in seq_along(pos_a)) {
    matr[[i]] <- .comb_wind_sets(indices_a = indices_a, indices_b = indices_b,
                                 pos_a = pos_a, pos_b = pos_b,
                                 min_dist = min_dist[i], max_dist = max_dist[i])
    matr[[i]][, 3L] <- i
  }
  matr <- do.call(what = rbind, args = matr)

  attr(x = matr, which = 'min_dist') <- min_dist
  attr(x = matr, which = 'max_dist') <- max_dist
  matr
}

