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
#include <cmath>
#include <limits>
#include <algorithm>
#include "../inst/include/declarations.hpp"
using namespace Rcpp;

// [[Rcpp::export(".sliding_window")]]
Rcpp::DataFrame sliding_window(arma::vec& x,
                               arma::vec& pos,
                               double start,
                               double width,
                               double advance,
                               std::string stat)
{
  double eps = std::numeric_limits<double>::epsilon();
  int n = pos.size();
  int n_wind = std::ceil((pos(n - 1) - start - width - eps) / advance);
  // Rcout << eps;
  // int n_wind = (int) (ceil( (arma::max(pos) - arma::min(pos) - width) - start) / advance);
  // Rcout << n << "  " << pos(n-1) << "  " << start << "  " << width << " " << advance;

  Rcpp::NumericVector begin_ret(n_wind, NA_REAL), end_ret(n_wind, NA_REAL), stat_ret(n_wind, NA_REAL), n_ret(n_wind, NA_REAL);
  // std::fill(stat_ret.begin(), stat_ret.end(), NA_REAL);
  // std::fill(begin_ret.begin(), begin_ret.end(), NA_REAL);
  // std::fill(end_ret.begin(), end_ret.end(), NA_REAL);
  // std::fill(n_ret.begin(), n_ret.end(), NA_REAL);

  int ixa = 0;
  int ixb;
  double begin = start;
  double end = begin + width;
  int bin = 0;
  while (bin < n_wind && (end = begin + width) <= pos[n - 1]) {

    ixa = index_geq(pos, begin);
    ixb = index_greater(pos, end) - 1;

    if (ixb < ixa) {
      stat_ret[bin] = NA_REAL;
    } else {
      arma::vec x_sub = x.subvec(ixa, ixb);
      if(stat == "mean"){
        stat_ret[bin] = arma::mean(x_sub);
      } else if(stat == "median"){
        stat_ret[bin] = arma::median(x_sub);
      } else if(stat == "min"){
        stat_ret[bin] = arma::min(x_sub);
      } else if(stat == "max"){
        stat_ret[bin] = arma::max(x_sub);
      } else if(stat == "sd"){
        stat_ret[bin] = arma::stddev(x_sub);
      }
    }
    begin_ret[bin] = begin;
    end_ret[bin] = end;
    n_ret[bin] = ixb - ixa + 1;
    begin = begin + advance;
    ++bin;
  }

  Rcpp::DataFrame ret = Rcpp::DataFrame::create(Rcpp::Named("begin") = begin_ret,
                                                Rcpp::Named("end") = end_ret,
                                                Rcpp::Named("stat") = stat_ret,
                                                Rcpp::Named("n") = n_ret);
  ret.attr("start") = start;
  ret.attr("width") = width;
  ret.attr("advance") = advance;
  return ret;
}

// //[[Rcpp::export(".sliding_window)]]
// NumericMatrix sliding_window(arma::vec & x,
//                              arma::vec & pos,
//                              double width,
//                              double advance,
//                              std::string stat)
// {
//   int n = pos.n_elem;
//   //int n_wind = (int) ceil( (pos(n-1) - pos(0) - width) / advance);
//   int n_wind = (int) ceil( (arma::max(pos) - arma::min(pos) - width) / advance);
//   //Rcout << "n_wind" << n_wind << std::endl;
//   arma::mat out(n_wind, 4);
//   std::fill(out.begin(), out.end(), NA_REAL);
//
//   double start = pos(0);
//   double end = start + width;
//
//   // Progress p(n_wind, true); // create the progress monitor
//
//   for (int bin = 0; bin != n_wind; bin++) {
//
// //     if (Progress::check_abort() ){
// //       stop("User aborted!\n");
// //     }
// //     p.increment(); // update progress
//
//     // arma::uvec q(1);
//     // q = arma::find(pos > start, 1, "first"); ixA = q(0);
//     // q = arma::find(pos < end, 1, "last"); ixB = q(0);
//     int ixA = std::lower_bound(pos.begin(), pos.end(), start) - pos.begin();
//     int ixB = std::upper_bound(pos.begin(), pos.end(), end) - pos.begin() - 1;
//
//     out(bin, 0) = start; // start of window
//     // only go on if at least one datum in window
//     if (ixA <= ixB) {
//
//     //Rcout << "bin= " << bin << "\tixA= " << ixA + 1 << "\tixB= " << ixB << "\n";
//
//       arma::running_stat<double> statVal;
//       // arma::running_stat<double> statPos;
//
//       for(int i = ixA; i <= ixB; i++){
//         double x = val(i);
//               statVal(x);
//               double p = pos(i);
//               statPos(p);
//       }
//
//       if(statistic == "mean"){
//         out(bin,1) = statVal.mean();
//       } else if(statistic == "median"){
//         Rcout << "Not implemented by Armadillo's running_stat class.\n";
//       } else if(statistic == "min"){
//         out(bin,1) = statVal.min();
//       } else if(statistic == "max"){
//         out(bin,1) = statVal.max();
//       } else if(statistic == "sd"){
//         out(bin,1) = statVal.stddev();
//       }
//
//       if(position == "mean"){
//         out(bin,2) = statPos.mean(); // average of marker positions
//       } else if (position == "min"){
//         out(bin,2) = statPos.min(); // positions of first marker in window
//       }
//
//       out(bin,3) = statVal.count();
//
//     }
//
//     start += advance;
//     end = start + width;
//   }
//   return wrap(out);
// }
