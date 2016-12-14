#' Compute summary statistics on a sliding window
#'
#' @description Compute summary statistics on a sliding window along a
#' vector containing data where the positions of the data points are
#' stored in another vector.
#' 
#' @param x A numeric vector containing the values.
#' @param pos A incresingly sorted vector containing the corresponding
#' positions.
#' @param start A double. The starting point of the sliding window.
#' @param width A double. The width of the sliding window.
#' @param advance A double. The step size of the sliding window.
#' @param stat A string. The summary statistic to be calculated.
#' 
#' @return A \code{data.frame} with columns \code{start} (start of the
#' window), \code{end} (end of the window), \code{stat} (the computed
#' summary statistic) and \code{n} (number of data points in the window).
#'
#' @details Note that only windows that fully fit into the range of \code{pos} are
#' considered.
#' 
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de}).
#' 
#' @examples
#' set.seed(123L)
#' n <- 1000L
#' x <- arima.sim(n = n, list(ar = 0.99))
#' pos <- sort(runif(n))
#' plot(x = pos, y = x)
#' advance <- 0.01
#' width <- c(0.02, 0.1, 0.2)
#' colors <- c('red', 'green', 'blue') 
#' for (i in seq_along(width)) {
#'   df <- sliding_window(x, pos, 0.0, width[i], advance, "mean")
#'   points(x = df$begin + width[i] / 2, df$stat, col = colors[i], pch = 19, lwd = 2.0, type = 'l')
#' }
#' @export
sliding_window <- function(x, pos, start,  width, advance,
                           stat = c('mean', 'median', 'min', 'max', 'sd')) {
  
  if (!is.atomic(x) || !is.numeric(x))
    stop("'x' must be a numeric vector.")
  
  if (!is.atomic(pos) || !is.numeric(pos) || is.unsorted(pos, strictly = TRUE))
    stop("'pos' must be an increasingly sorted numeric vector.")
  
  if (!is.numeric(start) || length(start) != 1L)
    stop("'start' must a double.")
  
  if (!is.numeric(width) || length(width) != 1L || width <= 0.0)
    stop("'width' must a positive double.")
  
  if (!is.numeric(advance) || length(advance) != 1L || advance <= 0.0)
    stop("'advance' must a positive double.")
  
  stat <- match.arg(arg = stat)
  ix <- !is.na(x)
  x <- x[ix]
  pos <- pos[ix]
  .sliding_window(x, pos, start, width, advance, stat)
}