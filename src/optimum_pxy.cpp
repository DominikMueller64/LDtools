// LDtools: Tools for Computation of Linkage Disequilibrium.
//
// Copyright (C) 2016 Dominik Mueller
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <utility>
#include <boost/math/tools/minima.hpp>

// [[Rcpp::export]]
arma::imat get_counts(const arma::vec& x,
                      const arma::vec& y){
  int N = x.size();
  arma::imat counts(3, 3, arma::fill::zeros);

  for (int i = 0; i != N; i++) {
    if (arma::is_finite(x(i)) && arma::is_finite(y(i))) {
      counts(2 - x(i), 2 - y(i)) += 1;
    }
  }
  return counts;
}

class objective_fun
{
public:

  double px, py;
  arma::imat n3x3;

  double operator()(double pxy) {
    double like =
      (2*n3x3(0,0) + n3x3(0,1) + n3x3(1,0)) * log(pxy) +
      (2*n3x3(0,2) + n3x3(0,1) + n3x3(1,2)) * log(px - pxy) +
      (2*n3x3(2,0) + n3x3(1,0) + n3x3(2,1)) * log(py - pxy) +
      (2*n3x3(2,2) + n3x3(2,1) + n3x3(1,2)) * log(1 - px - py + pxy) +
      n3x3(1,1) * log(pxy * (1 - px - py + pxy) + (px - pxy) * (py - pxy)) ;

    return -like;
  }
};

//' @export
// [[Rcpp::export(".optimum_pxy")]]
double optimum_pxy(const arma::imat& n3x3,
                   const double px,
                   const double py,
                   const double lower,
                   const double upper)
{

  objective_fun f;
  // double minloc; // location of the minimum
  // double minval; // value at the minimum
  f.px = px;
  f.py = py;
  f.n3x3 = n3x3;
  int bits_prec = 25;

  std::pair<double, double> r;
  r = boost::math::tools::brent_find_minima(f, lower, upper, bits_prec);
  return r.first;
}
