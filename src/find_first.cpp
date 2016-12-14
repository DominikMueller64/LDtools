// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace arma;

// Find index.
// 
// @description Find the index of the first element of a increasingly sorted
// vector that is larger than or equal to a given value.
//
// @param x A Vector. Increasingly sorted.
// @param y A Double.
// 
// @return A integer. The index of the first element in \code{x} which is 
// larger than or equal to \code{y}. If \code{y} is larger than all the
// elements of \code{x}, it returns \code{length(x) + 1}.
//
// [[Rcpp::export(".index_geq")]]
int index_geq(const arma::vec& x, const double y) {
  auto it = std::lower_bound(x.begin(), x.end(), y);
  return it - x.begin();
}

// [[Rcpp::export(".index_greater")]]
int index_greater(const arma::vec& x, const double y) {
  auto it = std::upper_bound(x.begin(), x.end(), y);
  return it - x.begin();
}

// [[Rcpp::export(".find_closest")]]
int find_closest(const arma::vec& x, const double y) {
  auto it = std::lower_bound(x.begin(), x.end(), y);
  if (it == x.begin())
    return 0;
  return (std::abs(*(it - 1) - y) < std::abs(*it - y) ? (it - 1) : it) - x.begin();
}
