#' Computation of linkage disequilibrium from biallelic markers.
#' 
#' This provides efficient functions for computing linkage
#' disequilibrium (LD) between biallelic marker loci from either phased
#' (haplotype) or unphased (genotype) data.
#' 
#' All performance-critical parts are implemented in C++.
#' 
#' @section Functions
#'
#' @docType package
#' @name LDtools
#' @useDynLib LDtools
#' @importFrom Rcpp sourceCpp
NULL