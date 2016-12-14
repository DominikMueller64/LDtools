#' Transform genotypes
#' 
#' @description Transform genotypes from two-letter coding to numeric
#' coding.
#' 
#' @param geno A character matrix containing the genotypes.
#' @param sep A character string use to split the genotypes in \code{geno}.
#'
#' @return A numeric matrix of the same dimensions as geno containing the
#' recoded genotypes.
#' 
#' @details Recoding is performed based on the major allele at the respective
#' locus, i.e., if for instance the major allele is "A", the recoded
#' genotypes of "A/A", "A/B", "B/A" and "B/B" will be 2, 1, 1, 0, given
#' \code{sep = "/"}.
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
transform_geno <- function (geno, sep = '/'){
  
  if (!is.matrix(geno) || !is.character(geno)) 
    stop("'geno' must be a character matrix.")
  
  if (!is.character(sep) || length(sep) != 1L)
    stop("'sep' must be a string.")
    
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

