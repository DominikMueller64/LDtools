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
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include "headers.h"

// Convert std::list to Rcpp::List
// @description Convert a std::list with arbitrary content to
// an Rcpp::List.
// @param stdlist A std::list.
// @return A Rcpp::List.
// @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
template<typename T>
Rcpp::List stdlist2List(T& stdlist)
{
  const int n = stdlist.size();
  Rcpp::List list(n);
  typename T::iterator iter;
  int i;
  for (iter = stdlist.begin(), i = 0; iter != stdlist.end(); ++iter, ++i) {
    list[i] = *iter;
  }
  return list;
}

// Convert std::list to arma::mat
// @description Convert a std::list containing matrices (arma::mat)
// to a single arma::mat.
// @param stdlist A std::list.
// @return A arma::mat.
// @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
arma::mat stdlist2mat(const std::list< arma::mat >& stdlist)
{
  std::list< arma::mat >::const_iterator it;
  int len = 0;
  arma::mat tmp;
  for (it = stdlist.begin(); it != stdlist.end(); ++it) {
    tmp = *it;
    len += tmp.n_rows;
  }
  arma::mat matr(len, tmp.n_cols);
  
  int ct = 0;
  for (it = stdlist.begin(); it != stdlist.end(); ++it) {
    tmp = *it;
    int n = tmp.n_rows;
    for (int i = 0; i < n; ++i, ++ct) {
      matr.row(ct) = tmp.row(i);
    }
  }
  return matr;
}

// Rcpp::NumericMatrix stdlist2Matrix(const std::list< Rcpp::NumericMatrix >& stdlist)
// {
//   // const int n = stdlist.size();
//   std::list< Rcpp::NumericMatrix >::const_iterator it;
//   int len = 0;
//   Rcpp::NumericMatrix tmp;
//   for (it = stdlist.begin(); it != stdlist.end(); ++it) {
//     tmp = *it;
//     len += tmp.nrow();
//   }
//   Rcpp::NumericMatrix matr(len, tmp.ncol());
//   
//   int ct = 0;
//   for (it = stdlist.begin(); it != stdlist.end(); ++it) {
//     tmp = *it;
//     int n = tmp.nrow();
//     for (int i = 0; i < n; ++i, ++ct) {
//       matr(ct, _) = tmp(i, _);
//     }
//   }
//   return matr;
// }


// Compute all combinations under spatial constrictions.
// 
// @description Compute all combinations of loci with a given minimum
// and maximum distance.
// 
// @param pos A incresingly sorted numeric vector with positions.
// @param min_dist A double giving the minimum distance between two loci.
// @param max_dist A double giving the maximim distance between two loci.
// 
// @return A list. Each element is a matrix that refers to a single locus and
// specifies all combinations of this locus with other loci fulfilling the
// criteria.
// 
// @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
// 
// [[Rcpp::export(".comb_wind")]]
arma::mat comb_wind(const arma::vec& pos,
                    const double min_dist,
                    const double max_dist)
{
  const int M = pos.n_elem;
  std::list < arma::mat > list;
  int len = 0;
  for (int i = 0; i < M; ++i) {
    int ixB1 = index_geq(pos, pos(i) + min_dist);
    int ixB2 = index_geq(pos, pos(i) + max_dist) - 1;
    if (ixB1 <= ixB2) {
      int n = ixB2 - ixB1 + 1;
      len += n;
      arma::mat matr(n, 3);
      int ct;
      for (int j = ixB1, ct = 0; j < ixB2 + 1; ++j, ++ct) {
        matr(ct, 0) = i + 1;
        matr(ct, 1) = j + 1;
        matr(ct, 2) = 1;
      }
      list.push_back(matr);
    }
  }
  return stdlist2mat(list);
}

// Compute all combinations of adjacent loci.
// 
// @description Compute all combination of loci with adjacent genetic
// positions.
// 
// @param n A integer. The number of loci.
// 
// @return A List. The list contains as its single element a matrix with all
// the combinations of adjacent loci.
// 
// @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
//
// [[Rcpp::export(".comb_adj")]]
arma::mat comb_adj(const arma::vec& pos)
{
  int n = pos.n_elem;
  arma::mat matr(n - 1, 3);
  for (int i = 0; i < n - 1; ++i) {
    matr(i, 0) = i + 1;
    matr(i, 1) = i + 2;
    matr(i, 2) = 1;
  }
  return matr;
}


// Compute all combinations of loci.
// 
// @description Compute all combination of loci.
// 
// @param n A integer. The number of loci.
// 
// @return A List. The list contains as its single element a matrix with all
// the combinations of loci.
// 
// @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
// 
// [[Rcpp::export(".comb_all")]]
arma::mat comb_all(const arma::vec& pos)
{
  int n = pos.n_elem;
  int len = (int) n * (n - 1) / 2;
  arma::mat matr(len, 3);
  int ct = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j, ++ct) {
      matr(ct, 0) = i + 1;
      matr(ct, 1) = j + 1;
      matr(ct, 2) = 1;
    }
  }
  return matr;
}


// // ' @examples
// // '
// // 'n1 <- 100
// // 'n2 <- 1000
// // 'pos1 <- sort(runif(n1))
// // 'pos2 <- sort(runif(n2))
// // 'all(dplyr::near(
// // '  apply(abs(outer(X = pos1, Y = pos2, FUN = `-`)), 1L, which.min),
// // '  .comb_nearest(1:n1, 1:n2, pos1, pos2)[[1L]][, 2]
// // '))
// // //' @export
// // // [[Rcpp::export(".comb_nearest")]]
// // arma::mat comb_nearest(const arma::ivec& indices,
// //                              const arma::ivec& indices2,
// //                              const arma::vec& pos,
// //                              const arma::vec& pos2)
// // {
// //   // std::sort(pos1.begin(), pos1.end(), std::greater<double>());
// //   int n = pos.n_elem;
// //   arma::mat matr(n, 3);
// //   for (int i = 0; i < n; ++i) {
// //     int j = find_closest(pos2, pos(i));
// //     // Rcout << pos1(i) << " " << j << std::endl;
// //     matr(i, 0) = indices(i);
// //     matr(i, 1) = indices2(j);
// //     matr(i, 2) = 1;
// //   }
// //   return matr;
// // }

// [[Rcpp::export(".comb_flank")]]
arma::mat comb_flank(const arma::ivec& indices,
                     const arma::ivec& indices2,
                     const arma::vec& pos,
                     const arma::vec& pos2)
  
{
  int n = pos.n_elem;
  int n2 = pos2.n_elem;
  std::list< arma::mat > stdlist(n);
  for (int i = 0; i < n; ++i) {
    double p = pos(i);  
    int ixr = index_geq(pos2, p);
    int ixl = ixr - 1;
    if (ixl > 0) {
      arma::mat matr(1, 3);
      matr(0, 0) = indices(i);
      matr(0, 1) = indices2(ixl);
      matr(0, 2) = i + 1;
      stdlist.push_back(matr);
    }
    if (ixr < n2) {
      arma::mat matr(1, 3);
      matr(0, 0) = indices(i);
      matr(0, 1) = indices2(ixr);
      matr(0, 2) = i + 1;
      stdlist.push_back(matr);
    }
  }
  return stdlist2mat(stdlist);
}

// Compute all combinations between closest loci.
// 
// @description Compute all combination between a set of loci and its closest
// pendants in another set.
// 
// @param indices1 A integer vector. A vector with the indices of the first
// set of loci.
// @param indices2 A integer vector. A vector with the indices of the second
// set of loci.
// @param pos1 A numeric vector. The positions of the first set of loci.
// @param pos2 A numeric vector. The positions of the second set of loci.
// 
// @return A List. Each row of the contained matrix refers to a locus of the first set of
// loci and specifies the pair of this locus and the locus
// closest to it from the second set of loci.
// 
// @details Vectors \code{pos1} and \code{pos2} must be increasingly sorted.
// This function is usefull if for a known set of QTL the LD between QTL and
// their closest markers should be computed.
// 
// @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
//
// [[Rcpp::export(".comb_nearest_k")]]
arma::mat comb_nearest_k(const arma::ivec& indices,
                         const arma::ivec& indices2,
                         const arma::vec& pos,
                         const arma::vec& pos2,
                         const int k)
{
  int n = pos.n_elem;
  int n2 = pos2.n_elem;
  arma::mat matr(k * n, 3);
  int ct = 0;
  for (int i = 0; i < n; ++i) {
    double p = pos(i);  
    int ixr = index_geq(pos2, p);
    int ixl = ixr - 1;
    for (int j = 0; j < k; ++j, ++ct) {
      matr(ct, 0) = indices(i);
      if (ixl < 0) {
        matr(ct, 1) = indices2(ixr);
        ++ixr;
      } else if (ixr > n2 - 1) {
        matr(ct, 1) = indices2(ixl);
        --ixl;
      }
      else {
        if (std::abs(pos2(ixl) - p) < std::abs(pos2(ixr) - p)) {
          matr(ct, 1) = indices2(ixl);
          --ixl;
        } else {
          matr(ct, 1) = indices2(ixr);
          ++ixr;
        }
      }
      matr(ct, 2) = i + 1; 
    }
  }
  return matr;
}


// Compute combination of loci in a sliding window.
// 
// @description Compute all combinations of loci in a moving sliding
// window.
// 
// @param pos An increasingly sorted numeric vector containing the genetic
// position of the loci.
//
// @param start A double, specifyin where to start with the sliding window.
// @param width A double, the width of the sliding window.
// @param advance A double, the increment of the sliding window.
// 
// @return A List. Contained matrices refer to individual windows where all
// possible combination of loci falling into this window are specified.
//
// @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
//
// [[Rcpp::export(".comb_sliding")]]
arma::mat comb_sliding(const arma::vec& pos,
                       const double start,
                       const double width,
                       const double advance)
{
  
  double eps = std::numeric_limits<double>::epsilon();
  int n = pos.n_elem;
  int n_wind = ceil((*(pos.end() - 1) - start - width - eps) / advance); 
  std::list< arma::mat > stdlist(n_wind);
  
  int ix1 = 0;
  int ix2;
  double begin = start;
  double end = begin + width;
  int wind = 0;
  int len = 0;
  while (wind < n_wind && end <= pos(n - 1)) {
    ix1 = index_geq(pos, begin);
    ix2 = index_greater(pos, end) - 1;
    // Rcout << begin << " " << end << " " << ix1 << " " << ix2 << std::endl;
    if (ix1 < ix2) {
      int m = ix2 - ix1 + 1;
      arma::mat matr((int) m * (m - 1) / 2, 3);
      int ct = 0;
      for (int i = ix1; i < ix2; ++i) {
        for (int j = i + 1; j < ix2 + 1; ++j, ++ct, ++len) {
          // Rcout << i << " " << j << std::endl;
          matr(ct, 0) = i + 1;
          matr(ct, 1) = j + 1;
          matr(ct, 2) = wind + 1;
        }
      }
      stdlist.push_back(matr);
      ++wind;
    } 
    begin += advance;
    end = begin + width;
  }
  return stdlist2mat(stdlist);
}
