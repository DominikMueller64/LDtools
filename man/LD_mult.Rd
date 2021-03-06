% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD_mult.R
\name{LD_mult}
\alias{LD_mult}
\title{Compute linkage disequilibrium (LD) among several loci}
\usage{
LD_mult(X, matr, pos = NULL, is_phased = TRUE, any_na = FALSE,
  r_only = FALSE, check = TRUE)
}
\arguments{
\item{X}{A numeric matrix. The genotypes of the loci. If the genotypes are
phased, these must only contain 0 and 1. If they are unphased, they must
only contain 0, 1, or 2, where 1 codes for the heterozygot.}

\item{matr}{A numeric Matrix with three columns.
The first two columns specify the indices of
the loci in \code{X} for which (LD) is computed.}

\item{pos}{An increasingly sorted numeric vector. The positions of the loci.}

\item{is_phased}{A logical. Are the data phased?}

\item{any_na}{A logical. May some genotypes contain missing values? If not,
computations are more efficient for phased genotypes.}

\item{r_only}{Should only the r-statistic be computed?
(Produces minimal output for saving memory and computation time.)}

\item{check}{A logical. Should checks be performed?}
}
\description{
Compute the linkage disequilibrium (LD) among several
combinations of loci.
}
\details{
For the documentation of the calculated LD statistics, see
\code{\link{LD}}.
}
\examples{
\dontrun{
data('population', package = 'LDtools')
pos <- map$pos

# Adjacent loci
mean(LD_mult(X, comb_adj(pos), check = FALSE)$r2, na.rm = TRUE)

# Nearest loci
i <- sort(sample(seq_along(pos), size = 1000L))
i2 <- setdiff(seq_along(pos), i)
mean(LD_mult(X, comb_nearest(i, i2, pos[i], pos[i2]), check = FALSE)$r2, na.rm = TRUE)
mean(LD_mult(X, comb_nearest_k(i, i2, pos[i], pos[i2], 2), check = FALSE)$r2, na.rm = TRUE)

# Flanking loci
mean(LD_mult(X, comb_flank(i, i2, pos[i], pos[i2]), check = FALSE)$r2, na.rm = TRUE)

# Sliding window
matr <- comb_sliding(pos, start = 0, width = 5, advance = 0.25)
dat <- LD_mult(X, matr, pos, check = FALSE)

dat <- by(data = dat, INDICES = factor(dat$block),
          FUN = function(x) {
          data.frame(pos = mean(c(x$pos1, x$pos2)),
                     r = mean(x$r, na.rm = TRUE))
          },
          simplify = FALSE)
dat <- do.call(rbind, dat)
with(dat, plot(pos, r, type = 'l'))
}

}
\seealso{
{ \code{\link{transform_geno}}, \code{\link{LD}},
\code{\link{comb_all}}, \code{\link{comb_adj}}, \code{\link{comb_nearest}},
\code{\link{comb_wind}}, \code{\link{comb_sliding}}. }
}
\author{
Dominik Mueller (\email{dominikmueller64@yahoo.de})
}
