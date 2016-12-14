#' @name population
#' @title A population for examples in LD computations.
#' @description The data set contains a simulated population
#' (one chromosome) for illustration LD computations.
#' @docType data
#' @format The matrix \code{X} contains haplotype data for 100 individuals
#' at 2000 simulated SNPs and 200 QTL. The data.frame \code{map} contains a
#' genetic map of the chromosome.
NULL

#' @name X
#' @title Haplotype data 
#' @description The haplotype data of the simulated population 
#' @docType data
#' @format A matrix with haplotype data for 100 diploid individuals at
#' 2000 simulated SNPs and 200 QTL. The alleles are coded with 0 and 1.
NULL

#' @name map
#' @title Genetic map
#' @description The genetic map of the chromsome.
#' @docType data
#' @format A data.frame with a the columns: \code{id} gives the unique
#' identifier of the locus, \code{pos} gives the genetic position and 
#' \code{is_qtl} indicates if the locus is a QTL.
NULL