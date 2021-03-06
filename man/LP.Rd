% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LP.R
\name{LP}
\alias{LP}
\title{Compute linkage phase similarity}
\usage{
LP(X1, X2, pos, min_dist, max_dist, method = c("correlation",
  "cosine_similarity", "sign_match"), signif_level = 1, is_phased = TRUE,
  any_na = FALSE, check = TRUE)
}
\arguments{
\item{X1}{A numeric matrix. Data for the first population.}

\item{X2}{A numeric matrix. Data for the second population.}

\item{pos}{An increasingly sorted numeric vector. The positions of the loci.}

\item{min_dist}{A numeric vector. Minimum distances.}

\item{max_dist}{A numeric vector. Maximum distances.}

\item{method}{A string. Method use for computing linkag phase
similarities, one of 'correlation', 'cosine_similarity' or
'sign_match'.}

\item{signif_level}{A double. If the chi-square test yields a p-value above
\code{signif_level} for any population, this pair of loci will be
excluded from the computation of linkage phases.}

\item{is_phased}{A logical. Are the data phased?}

\item{any_na}{A logical. May some genotypes contain missing values? If not,
computations are more efficient for phased genotypes.}

\item{check}{A logical. Should checks be performed?}
}
\value{
A data.frame with columns min_dist, max_dist, n (number of pairs),
and lps (linkage phase similarity coefficient).
}
\description{
Compute the linkage phase similary between two populations.
}
\details{
The computation is faster if \code{signif_level = 1} and all pairs are included.
}
\examples{
data('population', package = 'LDtools')
pos <- map$pos
mid <- nrow(X) \%/\% 2L
X1 <- X[1L:mid, , drop = FALSE]
X2 <- X[(mid + 1L):nrow(X), , drop = FALSE]
min_dist <- seq(sqrt(.Machine$double.eps), floor(max(pos) / 2), by = 1)
max_dist <- min_dist + 1

dat <- LP(X1, X2, pos, min_dist, max_dist, method = 'correlation')
with(dat,
 plot(x = (min_dist + max_dist) / 2, y = lps, type = 'o', ylim = c(0, 1))
)

}
\author{
Dominik Mueller (\email{dominikmueller64@yahoo.de})
}
