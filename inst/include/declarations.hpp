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

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

R_xlen_t index_greater(const arma::vec& x, const double y);

R_xlen_t index_geq(const arma::vec& x, const double y);

R_xlen_t find_closest(const arma::vec& x, const double y);

double optimum_pxy(const arma::imat& n3x3,
                   const double px,
                   const double py,
                  const double lower,
                   const double upper);

arma::imat get_counts(const arma::vec& x,
                      const arma::vec& y);
