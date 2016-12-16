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


#' Compute linkage disequilibrium (LD)
#'
#' @description Compute the linkage disequilibrium (LD) between a pair of loci.
#' 
#' @param x A numeric vector. The haplotype/genotype data of the first locus. If they are phased,
#' there must only contain 0 and 1. If they are unphased, they must only contain
#' 0, 1, or 2, where 1 codes for the heterozygot.
#' @param y A numeric vector. The haplotype/genotype data of the second locus.
#' 
#' @param is_phased A logical. Are the data phased?
#' @param any_na A logical. May some genotypes contain missing values? If not,
#' computations are more efficient for phased genotypes.
#' @param r_only Should only the r-statistic be computed?
#' (Produces minimal output for saving memory and computation time.)
#' @param check A logical. Should checks be performed?
#' 
#' @details For computing LD among a large number of loci, please use 
#' \code{\link{LD_mult}}, which is much more efficient for that.
#' 
#' @details
#' Linkage disequilibrium (LD) is the non-random association of
#' marker alleles and can arise from marker proximity or from selection
#' bias.
#'  
#' Three estimators of LD are computed:
#'  \itemize{
#'    
#'    \item{D}{ raw difference in frequency between the
#'       observed number of AB pairs and the expected number:
#'         
#'         \deqn{%
#'           D = p_{AB} - p_A p_B %
#'         }{%
#'           D = p(AB) - p(A)*p(B) %
#'         }
#'       
#'     }
#'     \item{D'}{ scaled D spanning the range [-1,1] 
#'       
#'      
#'      \deqn{D' = \frac{D}{D_{max} } }{D' = D / Dmax}
#'        
#'        where, if D > 0:
#'        \deqn{%
#'         D_{max} = \min( p_A p_b, p_a p_B )  %
#'         }{%
#'         Dmax = min( p(A)p(b), p(a)p(B) )   %
#'         } 
#'         or if D < 0:
#'         \deqn{%
#'         D_{max} = \max{ -p_A p_B, -p_a p_b }  %
#'         }{%
#'         Dmax = max( -p(A)p(B), -p(a)p(b) )  %
#'         }
#'     }
#'        
#'      \item{r}{ correlation coefficient between the markers
#'      
#'      \deqn{%
#'      r = \frac{-D}{\sqrt( p_A * p_a * p_B * p_b  )} %
#'       }{%
#'       r = -D / sqrt( p(A) * p(a) * p(B) * p(b) ) %
#'       }
#'       }
#'   }
#'         
#'    where
#'    \itemize{
#'    \item{-}{ \eqn{p_A}{p(A)} is defined as the observed probability of
#'    allele 'A' for marker 1, }
#'    \item{-}{ \eqn{p_a=1-p_A}{p(a) = 1-p(A)} is defined as the observed probability of
#'   allele 'a' for marker 1, }
#'   \item{-}{\eqn{p_B}{p(B)} is defined as the observed probability of
#'   allele 'B' for marker 2, and }
#'   \item{-}{\eqn{p_b=1-p_B}{p(b) = 1- p(B)} is defined as the observed probability of
#'   allele 'b' for marker 2, and }
#'    \item{-}{\eqn{p_{AB}}{p(AB)} is defined as the probability of
#'    the marker allele pair 'AB'. }
#'    }
#'    
#'    For genotype data, AB/ab cannot be distinguished from
#'    aB/Ab. Consequently, we estimate \eqn{p_{AB}}{p(AB)} using maximum
#'    likelihood and use this value in the computations.
#'  
#' @return
#' \code{LD} returns a list with the following components:
#' \itemize{
#'   \item{D}{ Linkage disequilibrium estimate}
#'   \item{Dprime}{ Scaled linkage disequilibrium estimate}
#'   \item{r}{ Correlation coefficient} 
#'   \item{r2}{ Squared correlation coefficient} 
#'   \item{n}{ Number of observations}
#'   \item{chis2}{ Chi-square statistic for linkage
#'     equilibrium (i.e., D = Dprime = r = r2 = 0)}
#'   \item{pval}{ Chi-square p-value for marker independence}
#' }
#' 
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#' The documentation is adapted from the \link[genetics]{LD} function of the 
#' \href{https://cran.r-project.org/web/packages/genetics/index.html}{genetics} 
#' package by \href{https://sites.google.com/a/warnes.net/resume/}{Gregory Warnes}.
#' @seealso{ \code{\link{transform_geno}}, \code{\link{LD_mult}}  }
#' 
#' @examples
#' # phased data
#' data('population', package = 'LDtools')
#' LD(x = X[, 1], y = X[, 2], is_phased = TRUE, any_na = FALSE)
#' 
#' # unphased data
#' tmp <- transform_geno(gen_geno(n = 20L, m = 2L))
#' LD(x = tmp$geno[, 1], y = tmp$geno[, 2], is_phased = FALSE,
#' any_na = TRUE)
#' 
#' @export
LD <- function(x, y, is_phased = TRUE, any_na = TRUE, 
               r_only = FALSE, check = TRUE) {
  
  if (check) {
    
    if (!is.atomic(x) || !is.numeric(x))
      stop("'x' must be a numeric vector.")
    
    if (!is.atomic(y) || !is.numeric(y))
      stop("'y' must be a numeric vector.")
    
    if ((n <- length(x)) != length(y))
      stop("'x' and 'y' must have the same length.")
    
    if (!is.logical(is_phased) || length(is_phased) != 1L)
      stop("'is_phased' must be logical and of length one.")
    
    if (!is.logical(any_na) || length(any_na) != 1L)
      stop("'any_na' must be logical and of length one.")
    
    if (!is.logical(r_only) || length(r_only) != 1L)
      stop("'r_only' must be logical and of length one.")
    
    uqx <- unique(x)
    uqy <- unique(y)
    if (is_phased) {
      
      if (length(uqx) + length(uqy) < 2L)
        stop("Data must be polymorphic")
      
      if (!all(union(uqx, uqy) %in% c(0, 1, NA_real_, NA_integer_)))
        stop("Phased data are only allowed to contain 0 and 1.")
      
    } else {
      
      if (all(uqx == 0) || all(uqx == 2) || all(uqy == 0) || all(uqy == 2))
        stop("Data must be polymorphic")
      
      if (!all(union(uqx, uqy) %in% c(0, 1, 2, NA_real_, NA_integer_)))
        stop("Unphased data are only allowed to contain 0, 1 and 2.")
    }
  }
  
  px <- mean(x, na.rm = TRUE)
  py <- mean(y, na.rm = TRUE)
  if (is_phased) {
    x <- x - px
    y <- y - px
  } else {
    px <- px / 2.0
    py <- py / 2.0
  }
  
  if(is.na(px) || is.na(py))
    stop("All gentoypes of either 'x' or 'y' are missing.")
  
  .LD(x = x, y = y, px = px, py = py, is_phased = is_phased,
      any_na = any_na, r_only = r_only)
}


