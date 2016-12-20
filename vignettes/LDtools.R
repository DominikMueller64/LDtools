## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- missing=FALSE, warning=FALSE---------------------------------------
library('magrittr') # %>% pipe operator
data('population', package = 'LDtools', envir = .GlobalEnv)

## ---- missing=FALSE, warning=FALSE---------------------------------------
library('LDtools')
LD(x = X[, 1], y = X[, 2], is_phased = TRUE, any_na = FALSE)

## ------------------------------------------------------------------------
tmp <- .gen_geno(n = 20L, m = 2L) %>% transform_geno()
LD(x = tmp$geno[, 1], y = tmp$geno[, 2], is_phased = FALSE, any_na = TRUE)

## ------------------------------------------------------------------------
pos <- map$pos
matr <- comb_wind(pos = pos, min_dist = 5, max_dist = 6)
dat <- LD_mult(X, matr, pos, is_phased = TRUE)

## ---- fig.height=4, fig.width=6------------------------------------------
matr <- comb_sliding(pos, start = 0, width = 5, advance = 1)
dat <- LD_mult(X, matr, pos)
# Summarize data for plotting.
dat <- by(data = dat, INDICES = factor(dat$block),
          FUN = function(x) {
          data.frame(pos = mean(c(x$pos1, x$pos2)),
                     r2 = mean(x$r2, na.rm = TRUE))
          },
          simplify = FALSE) %>% do.call(what = rbind)

with(dat, plot(pos, r2, type = 'l'))

## ---- fig.height=4, fig.width=6------------------------------------------
mid <- nrow(X) %/% 2L
X1 <- X[1L:mid, , drop = FALSE]
X2 <- X[(mid + 1L):nrow(X), , drop = FALSE]
min_dist <- seq(sqrt(.Machine$double.eps), floor(max(pos) / 2), by = 1)
max_dist <- min_dist + 1

dat <- LP(X1, X2, pos, min_dist, max_dist, method = 'correlation')
with(dat,
 plot(x = (min_dist + max_dist) / 2, y = lps, type = 'o', ylim = c(0, 1))
)

