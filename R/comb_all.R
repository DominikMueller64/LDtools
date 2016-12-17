##' @title Combinations among all loci.
##' @description Combinations among all loci.
##' @param pos A numeric vector. Increasingly sorted, contains the positions of loci.
##' @return A numeric matrix.
##' @examples
##' NULL 
##' @author dominik
comb_all <- function(pos) {
  if (!is.atomic(pos) || !is.numeric(pos) || is.unsorted(pos))
    stop("'pos' must be an increasingly sorted vector.")
  .comb_all(pos)
}

##' @title Combinations among all loci of two sets.
##' @description Combinations among all loci of two sets.
##' @inheritParams comb_wind_sets
##' @return A numeric matrix.
##' @examples
##' NULL 
##' @author dominik
comb_all_sets <- function(indices_a, indices_b, pos_a, pos_b) {
  .comb_all_sets(indices_a = indices_a, indices_b = indices_b,
                 pos_a = pos_a, pos_b = pos_b)
}


