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
  const int n = pos.size();
  std::list < arma::mat > list;
  int len = 0;
  for (int i = 0; i < n; ++i) {
    int ix1 = index_geq(pos, pos(i) + min_dist);
    int ix2 = index_geq(pos, pos(i) + max_dist) - 1;
    if (ix1 <= ix2) {
      arma::mat matr(ix2 - ix1 + 1, 3);
      int ct = 0;
      for (int j = ix1; j < ix2 + 1; ++j, ++ct) {
        matr(ct, 0) = i + 1;
        matr(ct, 1) = j + 1;
        matr(ct, 2) = 1;
      }
      len += ct;
      list.push_back(matr);
    }
  }
  return stdlist2mat(list);
}

// [[Rcpp::export(".comb_wind_sets")]]
arma::mat comb_wind_sets(const arma::ivec& indices_a,
                         const arma::ivec& indices_b,
                         const arma::vec& pos_a,
                         const arma::vec& pos_b,
                         const double min_dist,
                         const double max_dist)
{
  Rcpp::Rcout << "motherfucker" << std::endl;
  const int n_a = pos_a.size();
  std::list < arma::mat > list;
  int len = 0;
  for (int i = 0; i < n_a; ++i) {
    int ix1 = index_geq(pos_b, pos_a(i) + min_dist);
    int ix2 = index_geq(pos_b, pos_a(i) + max_dist) - 1;
    if (ix1 <= ix2) {
      arma::mat matr(ix2 - ix1 + 1, 3);
      int ct = 0;
      for (int j = ix1; j < ix2 + 1; ++j, ++ct) {
        matr(ct, 0) = indices_a(i + 1);
        matr(ct, 1) = indices_b(j + 1);
        matr(ct, 2) = 1;
      }
      len += ct;
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

arma::mat comb_all_sets(const arma::ivec& indices_a,
                        const arma::ivec& indices_b,
                        const arma::vec& pos_a,
                        const arma::vec& pos_b)
{
  int n_a = pos_a.n_elem;
  int n_b = pos_b.n_elem;
  arma::mat matr(n_a * n_b, 3);
  int ct = 0;
  for (int i = 0; i < n_a - 1; ++i) {
    for (int j = 0; j < n_b; ++j, ++ct) {
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
// // 'pos_b <- sort(runif(n2))
// // 'all(dplyr::near(
// // '  apply(abs(outer(X = pos1, Y = pos2, FUN = `-`)), 1L, which.min),
// // '  .comb_nearest(1:n1, 1:n2, pos1, pos2)[[1L]][, 2]
// // '))
// // //' @export
// // // [[Rcpp::export(".comb_nearest")]]
// // arma::mat comb_nearest(const arma::ivec& indices_a,
// //                              const arma::ivec& indices_b,
// //                              const arma::vec& pos,
// //                              const arma::vec& pos2)
// // {
// //   // std::sort(pos1.begin(), pos1.end(), std::greater<double>());
// //   int n = pos.n_elem;
// //   arma::mat matr(n, 3);
// //   for (int i = 0; i < n; ++i) {
// //     int j = find_closest(pos2, pos(i));
// //     // Rcout << pos1(i) << " " << j << std::endl;
// //     matr(i, 0) = indices_a(i);
// //     matr(i, 1) = indices_b(j);
// //     matr(i, 2) = 1;
// //   }
// //   return matr;
// // }

// [[Rcpp::export(".comb_flank_sets")]]
arma::mat comb_flank_sets(const arma::ivec& indices_a,
                          const arma::ivec& indices_b,
                          const arma::vec& pos_a,
                          const arma::vec& pos_b)
  
{
  int n_a = pos_a.n_elem;
  int n_b = pos_b.n_elem;
  std::list< arma::mat > stdlist(n_a);
  for (int i = 0; i < n_a; ++i) {
    int ixr = index_geq(pos_b, pos_a(i));
    int ixl = ixr - 1;
    if (ixl > 0) {
      arma::mat matr(1, 3);
      matr(0, 0) = indices_a(i);
      matr(0, 1) = indices_b(ixl);
      matr(0, 2) = i + 1;
      stdlist.push_back(matr);
    }
    if (ixr < n_b) {
      arma::mat matr(1, 3);
      matr(0, 0) = indices_a(i);
      matr(0, 1) = indices_b(ixr);
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
// @param indices_a A integer vector. A vector with the indices of the first
// set of loci.
// @param indices_b A integer vector. A vector with the indices of the second
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
// [[Rcpp::export(".comb_nearest_k_sets")]]
arma::mat comb_nearest_k_sets(const arma::ivec& indices_a,
                              const arma::ivec& indices_b,
                              const arma::vec& pos_a,
                              const arma::vec& pos_b,
                              const int k)
{
  int n = pos_a.n_elem;
  int n2 = pos_b.n_elem;
  arma::mat matr(k * n, 3);
  int ct = 0;
  for (int i = 0; i < n; ++i) {
    double p = pos_a(i);  
    int ixr = index_geq(pos_b, p);
    int ixl = ixr - 1;
    for (int j = 0; j < k; ++j, ++ct) {
      matr(ct, 0) = indices_a(i);
      if (ixl < 0) {
        matr(ct, 1) = indices_b(ixr);
        ++ixr;
      } else if (ixr > n2 - 1) {
        matr(ct, 1) = indices_b(ixl);
        --ixl;
      }
      else {
        if (std::abs(pos_b(ixl) - p) < std::abs(pos_b(ixr) - p)) {
          matr(ct, 1) = indices_b(ixl);
          --ixl;
        } else {
          matr(ct, 1) = indices_b(ixr);
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
  //int n = pos.n_elem;
  double max_pos = *(pos.end() - 1);
  int n_wind = ceil((max_pos - start - width - eps) / advance); 
  std::list< arma::mat > stdlist(n_wind);
  
  int ix1 = 0;
  int ix2;
  double begin = start;
  double end = begin + width;
  int wind = 0;
  int len = 0;
  while (wind < n_wind && end <= max_pos) {
    ix1 = index_geq(pos, begin);
    ix2 = index_greater(pos, end) - 1;
    if (ix1 < ix2) {
      int m = ix2 - ix1 + 1;
      arma::mat matr((int) m * (m - 1) / 2, 3);
      int ct = 0;
      for (int i = ix1; i < ix2; ++i) {
        for (int j = i + 1; j < ix2 + 1; ++j, ++ct, ++len) {
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


// [[Rcpp::export(".comb_sliding_set")]]
arma::mat comb_sliding_sets(const arma::ivec& indices_a,
                            const arma::ivec& indices_b,
                            const arma::vec& pos_a,
                            const arma::vec& pos_b,
                            const double start,
                            const double width,
                            const double advance)
{
  
  double eps = std::numeric_limits<double>::epsilon();
  double last_pos_a = *(pos_a.end() - 1);
  double last_pos_b = *(pos_b.end() - 1);
  double min_pos = std::min(last_pos_a, last_pos_b);
  int n_wind = ceil((min_pos - start - width - eps) / advance); 
  std::list< arma::mat > stdlist(n_wind);
  
  int ix1_a = 0;
  int ix1_b = 0;
  int ix2_a, ix2_b;
  double begin = start;
  double end = begin + width;
  int wind = 0;
  int len = 0;

  while (wind < n_wind && end <= min_pos) {
    ix1_a = index_geq(pos_a, begin);
    ix2_a = index_greater(pos_a, end) - 1;
    ix1_b = index_geq(pos_b, begin);
    ix2_b = index_greater(pos_b, end) - 1;
    if (ix1_a <= ix2_a && ix1_b <= ix2_b) {
      int m_a = ix2_a - ix1_a + 1;
      int m_b = ix2_b - ix1_b + 1;
      arma::mat matr(m_a * m_b, 3);
      int ct = 0;
      for (int i = ix1_a; i < ix2_a; ++i) {
        for (int j = ix1_b; i < ix2_b; ++j, ++ct, ++len) {
          matr(ct, 0) = indices_a(i + 1);
          matr(ct, 1) = indices_b(j + 1);
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
