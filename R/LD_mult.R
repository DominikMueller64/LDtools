#' Compute linkage disequilibrium (LD) among several loci
#' 
#' @description Compute the linkage disequilibrium (LD) among several
#' combinations of loci.
#' 
#' @param X A numeric matrix. The genotypes of the loci. If the genotypes are
#' phased, these must only contain 0 and 1. If they are unphased, they must
#' only contain 0, 1, or 2, where 1 codes for the heterozygot.
#' 
#' @param matr A numeric Matrix with three columns.
#' The first two columns specify the indices of
#' the loci in \code{X} for which (LD) is computed. 
#' 
#' @param pos An increasingly sorted numeric vector. The positions of the loci. 
#' 
#' @inheritParams LD
#' 
#' @details For the documentation of the calculated LD statistics, see
#' \code{\link{LD}}.
#' 
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#' @seealso{ \code{\link{transform_geno}}, \code{\link{LD}},
#' \code{\link{comb_all}}, \code{\link{comb_adj}}, \code{\link{comb_nearest}},
#' \code{\link{comb_wind}}, \code{\link{comb_sliding}}. }
#'        
#' @examples
#' \dontrun{ 
#' data('population', package = 'LDtools')
#' pos <- map$pos
#' 
#' # Adjacent loci
#' mean(LD_mult(X, comb_adj(pos), check = FALSE)$r2, na.rm = TRUE)
#' 
#' # Nearest loci
#' i <- sort(sample(seq_along(pos), size = 1000L))
#' i2 <- setdiff(seq_along(pos), i)
#' mean(LD_mult(X, comb_nearest(i, i2, pos[i], pos[i2]), check = FALSE)$r2, na.rm = TRUE)
#' mean(LD_mult(X, comb_nearest_k(i, i2, pos[i], pos[i2], 2), check = FALSE)$r2, na.rm = TRUE)
#' 
#' # Flanking loci
#' mean(LD_mult(X, comb_flank(i, i2, pos[i], pos[i2]), check = FALSE)$r2, na.rm = TRUE)
#' 
#' # Sliding window
#' matr <- comb_sliding(pos, start = 0, width = 5, advance = 0.25)
#' dat <- LD_mult(X, matr, pos, check = FALSE)
#' 
#' dat <- by(data = dat, INDICES = factor(dat$block),
#'           FUN = function(x) {
#'           data.frame(pos = mean(c(x$pos1, x$pos2)),
#'                      r = mean(x$r, na.rm = TRUE))
#'           },
#'           simplify = FALSE)
#' dat <- do.call(rbind, dat)
#' with(dat, plot(pos, r, type = 'l'))
#' }
#'  
#' @export
LD_mult <- function(X, matr, pos = NULL, is_phased = TRUE,
                    any_na = FALSE, r_only = FALSE, check = TRUE){
  
  m <- ncol(X)
  n <- nrow(X)
  if (check) {
    
    if (!is.matrix(X) || !is.numeric(X))
      stop("'X' must be a numeric matrix.")
    
    if (!is.matrix(matr) || !is.numeric(matr))
      stop("'matr' must be a numeric matrix")
    
    if (any((tmp <- as.vector(matr[, c(1L, 2L)])) > ncol(X)) || any(tmp < 1))
      stop("'matr' must specify valid indices of loci.")
    
    if (!is.logical(is_phased) || length(is_phased) != 1L)
      stop("'is_phased' must be logical and of length one.")
    
    if (!is.logical(any_na) || length(any_na) != 1L)
      stop("'any_na' must be logical and of length one.")
    
    if (!is.logical(r_only) || length(r_only) != 1L)
      stop("'r_only' must be logical and of length one.")
    
    if (n < 2L)
      stop("'X' has less than 2 rows.")
  
    uq <- stats::na.exclude(unique(c(X)))
    
    if (is_phased) {
      if (!all(uq %in% c(0, 1)))
        stop("'X' is not a \"phased\" matrix (all element 0 or 1).")
    } else {
      if (!all(uq %in% c(0, 1, 2)))
        stop("'X' is not a \"unphased\" matrix (all element 0 or 1 or 2).")
    }
  }
  
  p <- .colMeans(X, n, m, na.rm = TRUE)
  if (is_phased) {
    X <- scale(X, center = p, scale = FALSE) # replace this by efficient version.
  } else {
    p <- p / 2.0
  }
  
  if (r_only) {
    ret <- .LD_mult_r(X = X, p = p, matr = matr, is_phased = is_phased, any_na = any_na) 
  } else {
    ret <- .LD_mult(X = X, p = p, matr = matr, is_phased = is_phased, any_na = any_na) 
  }
  if (!is.null(pos)) {
    pos1 <- pos[matr[, 1L, drop = TRUE]]
    pos2 <- pos[matr[, 2L, drop = TRUE]] 
    ret$pos1 <- pos1
    ret$pos2 <- pos2
    ret$dist <- pos2 - pos1
  }
  ret
}
