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
// #include <Rcpp.h>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <algorithm>
// using namespace Rcpp;
using namespace Rcpp;
using namespace arma;


// The return value of Minimize is the minimum of the function f.
// The location where f takes its minimum is returned in the variable minLoc.
// Notation and implementation based on Chapter 5 of Richard Brent's book
// "Algorithms for Minimization Without Derivatives".

template <class TFunction>
double Minimize
  (
      TFunction& f,  	// [in] objective function to minimize
      double leftEnd,     // [in] smaller value of bracketing interval
      double rightEnd,    // [in] larger value of bracketing interval
      double epsilon,     // [in] stopping tolerance
      double& minLoc      // [out] location of minimum
  )
{
  double d, e, m, p, q, r, tol, t2, u, v, w, fu, fv, fw, fx;
  static const double c = 0.5*(3.0 - sqrt(5.0));
  static const double SQRT_DBL_EPSILON = sqrt(DBL_EPSILON);
  
  double& a = leftEnd; double& b = rightEnd; double& x = minLoc;
  
  v = w = x = a + c*(b - a); d = e = 0.0;
  fv = fw = fx = f(x);
  int counter = 0;
  loop:
    counter++;
  m = 0.5*(a + b);
  tol = SQRT_DBL_EPSILON*fabs(x) + epsilon; t2 = 2.0*tol;
  // Check stopping criteria
  if (fabs(x - m) > t2 - 0.5*(b - a))
  {
    p = q = r = 0.0;
    if (fabs(e) > tol)
    {
      // fit parabola
      r = (x - w)*(fx - fv);
      q = (x - v)*(fx - fw);
      p = (x - v)*q - (x - w)*r;
      q = 2.0*(q - r);
      (q > 0.0) ? p = -p : q = -q;
      r = e; e = d;
    }
    if (fabs(p) < fabs(0.5*q*r) && p < q*(a - x) && p < q*(b - x))
    {
      // A parabolic interpolation step
      d = p/q;
      u = x + d;
      // f must not be evaluated too close to a or b
      if (u - a < t2 || b - u < t2)
        d = (x < m) ? tol : -tol;
    }
    else
    {
      // A golden section step
      e = (x < m) ? b : a;
      e -= x;
      d = c*e;
    }
    // f must not be evaluated too close to x
    if (fabs(d) >= tol)
      u = x + d;
    else if (d > 0.0)
      u = x + tol;
    else
      u = x - tol;
    fu = f(u);
    // Update a, b, v, w, and x
    if (fu <= fx)
    {
      (u < x) ? b = x : a = x;
      v = w; fv = fw; 
      w = x; fw = fx; 
      x = u; fx = fu;
    }
    else
    {
      (u < x) ? a = u : b = u;
      if (fu <= fw || w == x)
      {
        v = w; fv = fw; 
        w = u; fw = fu;
      }
      else if (fu <= fv || v == x || v == w)
      {
        v = u; fv = fu;
      }
    }
    goto loop;  // Yes, the dreaded goto statement. But the code here is faithful to Brent's orginal pseudocode.
  }
  return  fx;
}

// //' @export
// // [[Rcpp::export(".optimum_pxy")]] 
// double optimum_pxy(const IntegerMatrix n3x3, 
//                    const double px, 
//                    const double py,
//                    const double lower,
//                    const double upper)
// {
//   objective_fun f;
//   double minloc; // location of the minimum
//   double minval; // value at the minimum
//   f.px = px;
//   f.py = py;
//   f.n3x3 = n3x3;
//   minval = Minimize<objective_fun>(f, lower, upper, 1e-6, minloc);
//   return minloc;
// }

// //' @export
// // [[Rcpp::export(".get_counts")]] 
// IntegerMatrix get_counts(const arma::vec& x,
//                          const arma::vec& y){
//   int N = x.n_elem;
//   IntegerMatrix counts(3, 3);
//   std::fill(counts.begin(), counts.end(), 0);
//   
//   for (int i = 0; i != N; i++) {
//     if (is_finite(x(i)) && is_finite(y(i))) {
//       counts(2 - x(i), 2 - y(i)) += 1;
//     }
//   }
//   return(counts);
// }


// // [[Rcpp::export]]
// List estimate_LD(const NumericVector& a,
//                 const NumericVector& b,
//                 double px, double py) {
//   
//   // compute minor allele frequencies
//   double px = 1.0 - px;
//   double qy = 1.0 - py;
//   
//   // compute stistics for calculating Dprime
//   double Dmin = std::max(-px * py, -px * qy);
//   double pmin = px * py + Dmin;
//   
//   double Dmax = std::min(px * qy, py * qx);
//   double pmax = px * py + Dmax;
//   
//   // compute the count (n3x3) matrix
//   IntegerMatrix n3x3 = getCounts(a, b);
//   
//   // find maximum likelihood value for pxy
//   double pxy = optimum_pxy(n3x3, px, py, pmin, pmax);
//   
//   double estD = pxy - px * py;
//   double estDp = 0.0;
//   if(estD > 0.0){
//     estDp = estD / Dmax;
//   } else {
//     estDp = estD / Dmin;
//   }
//   
//   int n = Rcpp::sum(n3x3);
//   double cor = estD / sqrt( px * py * qx * qy );
//   double dchi = (2 * n * pow(estD, 2.0)) / (px * qx * py * qy);
//   double dpval = 1.0 - pchisq(NumericVector::create(dchi), 1.0, true, false)[0];
//   
//   return List::create(Named("D") = estD,
//                       Named("Dprime") = estDp,
//                       Named("r") = cor,
//                       Named("r2") = pow(cor, 2.0),
//                       Named("n") = n,
//                       Named("dchi") = dchi,
//                       Named("dpval") = dpval);
// }



// // [[Rcpp::export]]
// List estimate_LD(const arma::vec& a,
//                 const arma::vec& b,
//                 double px, double py) {
//   
//   // compute minor allele frequencies
//   double qx = 1.0 - px;
//   double qy = 1.0 - py;
//   
//   // compute stistics for calculating Dprime
//   double Dmin = std::max(-px * py, -qx * qy);
//   double pmin = px * py + Dmin;
//   
//   double Dmax = std::min(px * qy, py * qx);
//   double pmax = px * py + Dmax;
//   
//   // compute the count (n3x3) matrix
//   IntegerMatrix n3x3 = getCounts(a, b);
//   
//   // find maximum likelihood value for pxy
//   double pxy = optimum_pxy(n3x3, px, py, pmin, pmax);
//   
//   double estD = pxy - px * py;
//   double estDp = 0.0;
//   if(estD > 0.0){
//     estDp = estD / Dmax;
//   } else {
//     estDp = estD / Dmin;
//   }
//   
//   int n = Rcpp::sum(n3x3);
//   double cor = estD / sqrt( px * py * qx * qy );
//   double dchi = (2 * n * pow(estD, 2.0)) / (px * qx * py * qy);
//   double dpval = 1.0 - pchisq(NumericVector::create(dchi), 1.0, true, false)[0];
//   
//   return List::create(Named("D") = estD,
//                       Named("Dprime") = estDp,
//                       Named("r") = cor,
//                       Named("r2") = pow(cor, 2.0),
//                       Named("n") = n,
//                       Named("dchi") = dchi,
//                       Named("dpval") = dpval);
// }

// //' @export
// // [[Rcpp::export]] 
// DataFrame LDwithMLE(NumericMatrix cMajMat,
//                     NumericVector majFreqs){
//   
//   
//   int M = cMajMat.ncol();
//   int len = M * (M - 1) / 2;
//   
//   IntegerVector ix1(len);
//   IntegerVector ix2(len);
//   NumericVector D(len);
//   NumericVector Dprime(len);
//   NumericVector r(len);
//   NumericVector r2(len);
//   NumericVector dchi(len);
//   NumericVector dpval(len);
//   IntegerVector n(len);
//   
//   // NumericMatrix out(M*(M-1)/2, 9);
//   // std::fill(out.begin(), out.end(), NA_REAL);
//   
//   int ct = 0;
//   for(int i = 0; i != M-1; i++){
//     for(int j = i + 1; j != M; j++){
//       
//       // extract allele frequncies of major allele
//       double px = majFreqs(i);
//       double py = majFreqs(j);
//       
//       List tmp = estimate_LD(cMajMat(_,i), cMajMat(_,j), px, py);
//       ix1[ct] = i + 1;
//       ix2[ct] = j + 1;
//       D[ct] = tmp["D"];
//       Dprime[ct] = tmp["Dprime"];
//       r[ct] = tmp["r"];
//       r2[ct] = tmp["r2"];
//       n[ct] = tmp["n"];
//       dchi[ct] = tmp["dchi"];
//       dpval[ct] = tmp["dpval"];
      
      // // compute minor allele frequencies
      // double qx = 1.0 - px;
      // double qy = 1.0 - py;
      // 
      // // compute stistics for calculating Dprime
      // double Dmin = std::max(-px*py, -qx*qy);
      // double pmin = px*py + Dmin;
      // 
      // double Dmax = std::min(px*qy, py*qx);
      // double pmax = px*py + Dmax;
      // 
      // 
      // // compute the count (n3x3) matrix
      // IntegerMatrix n3x3 = getCounts(cMajMat(_,i), cMajMat(_,j));
      // 
      // // find maximum likelihood value for pxy
      // double pxy = optimum_pxy(n3x3, px, py, pmin, pmax);
      // 
      // double estD = pxy - px * py;
      // double estDp = 0.0;
      // if(estD > 0.0){
      //   estDp = estD / Dmax;
      // } else {
      //   estDp = estD / Dmin;
      // }
      // 
      // int n = Rcpp::sum(n3x3);
      // double corr = estD / sqrt( px * py * qx * qy );
      // double dchi = (2 * n * pow(estD, 2.0)) / (px * qx * py * qy);
      // double dpval = 1.0 - pchisq(NumericVector::create(dchi), 1.0, true, false)[0];
      // 
      // 
      // 
      // out(ct,0) = i; // marker 1
      // out(ct,1) = j; // marker 2
      // out(ct,2) = estD; // D
      // out(ct,3) = estDp; // Dprime
      // out(ct,4) = corr; // correlation
      // out(ct,5) = pow(corr, 2.0); // LD
      // out(ct,6) = n; // n
      // out(ct,7) = dchi; // chi-square test statistic
      // out(ct,8) = dpval; // p-value
//         
//       ct++;
//     }
//   }
//   // return(out);
//   
//   return DataFrame::create(Named("ix1") = ix1, Named("ix2") = ix2,
//                            Named("D") = D, Named("Dprime") = Dprime,
//                            Named("r") = r, Named("r2") = r2,
//                            Named("n") = n, Named("dchi") = dchi,
//                            Named("dpval") = dpval);
// }
