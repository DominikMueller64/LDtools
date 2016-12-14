// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include "headers.h"

// Center matrix.
//
// @description Center a matrix in-place.
// @param x A matrix of type double.
// @param center A vector of type double for centering
// @return Nothing.
// @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
// @examples
// x <- matrix(1:20 * 1.0, nrow = 4)
// p <- colMeans(x)
// y <- LDtools:::center_matrix(x, center = p)
//
// [[Rcpp::export(".center_matrix")]]
void center_matrix(arma::mat& X,
                   const arma::vec& center) {
  int n = X.n_rows;
  int m = X.n_cols;
  for (int j = 0; j < m; ++j) {
    double c = center(j);
    for (int i = 0; i < n; ++i) {
      X(i, j) = X(i, j) - c;
    }
  }
}


// Compute the inner producte.
// 
// @description This Function compute the inner product between two
// vector, catering for missing values.
// 
// @param x A numeric vector.
// @param y Another numeric vector.
// @param n Length of either vectors.
// 
// @details This function excludes all entries where either the element of
// \code{x} or the element of \code{y} contains a missing value. It
// works like \code{sum(x*y, na.rm = TRUE)}.
// 
// @return A List. The inner product and the number of finite products from
// which it was computed. If there are less than 2 finite finite,
// it returns \code{NA}.
// 
// @details This auxiliary function is intended only for internal usage.
// 
// @author Dominik Müller (\email{dominikmueller64@@yahoo.de})
struct inner_out
{
  double result;
  int n;
};

inner_out inner(const arma::vec& x,
                const arma::vec& y) {
  
  double result = 0;
  int n_fin = 0;
  auto itx = x.begin();
  auto ity = y.begin();
  while (itx != x.end()) {
    double tmp = *itx++ * *ity++;
    if (arma::is_finite(tmp)) {
      result += tmp;
      ++n_fin;
    }
  }
  if (n_fin < 2) {
    result = NA_REAL;
  }
  inner_out out = {result, n_fin};
  return out;
}

// Structure for returning multiple values in function LD.
struct LD_out {
  int   n;
  double D;
  double Dprime;
  double r;
  double r2;
  double chi2;
  double pval;
};

LD_out LD_struct(const arma::vec& x,
                 const arma::vec& y,
                 const double px,
                 const double py,
                 const bool any_na,
                 const bool is_phased)
{
  LD_out ret;
  double D, sum_prod, Dprime, r;
  int n = x.n_elem;
  double qx = 1 - px;
  double qy = 1 - py;
  
  double Dmin = std::max(-px * py, -qx * qy);
  double Dmax = std::min(px * qy, qx * py);
  
  if (is_phased) {
    if (any_na) {
      inner_out tmp = inner(x, y);
      n = tmp.n;
      sum_prod = tmp.result;
    } else {
      sum_prod = arma::dot(x, y);
    }
    D = sum_prod / n;
  } else {
    // pA = px, pB = py, pa = qx, pb = qy 
    arma::imat n3x3 = get_counts(x, y); // compute the count (n3x3) matrix
    double pxy = optimum_pxy(n3x3, px, py, px * py + Dmin,
                             px * py + Dmax); // max. likel. for pAB
    D = pxy - px * py;
    n = arma::accu(n3x3);
  }
  Dprime = D / ((D > 0) ? Dmax : Dmin);
  
  double prod = px * qx * py * qy;
  r = D / sqrt(prod);
  double chi2 = (2 * n * pow(D, 2)) / prod;
  double pval = 1 - R::pchisq(chi2, 1, true, false);
  
  ret.n = n;
  ret.D = D;
  ret.Dprime = Dprime;
  ret.r = r;
  ret.r2 = pow(r, 2);
  ret.chi2 = chi2;
  ret.pval = pval;
  return ret;
}

double LD_r(const arma::vec& x,
                   const arma::vec& y,
                   const double px,
                   const double py,
                   const bool any_na,
                   const bool is_phased)
{
  int n = x.n_elem;
  double D, sum_prod;
  double qx = 1 - px;
  double qy = 1 - py;
  double Dmin = std::max(-px * py, -qx * qy);
  double Dmax = std::min(px * qy, qx * py);
  
  if (is_phased) {
    if (any_na) {
      inner_out tmp = inner(x, y);
      sum_prod = tmp.result;
      n = tmp.n;
    } else {
      sum_prod = arma::dot(x, y);
    }
    D = sum_prod / n;
  } else {
    arma::imat n3x3 = get_counts(x, y);
    double pxy = optimum_pxy(n3x3, px, py, px * py + Dmin,
                             px * py + Dmax);
    D = pxy - px * py;
  }
  return D / sqrt(px * qx * py * qy);
}

// Compute LD between two loci.
// 
// @description Compute different LD statistics between two loci.
//
// @param x A numeric Vector. The genotypes at the first locus.
// @param y A numeric Vector. The genotypes at the second locus.
// @param px A double. The allele frequency at the first locus.
// @param py A double. The allele frequency at the second locus.
// @param any_na A logical. May the data contain any missing values?
// @param is_phased A logical. Are the genotypes referring to phased
// or to unphased data?
//
// @return A List containing the calculated LD statistics.
// 
// @details If the genotypes are phased, they MUST be centered by the 
// allele frequency.
// 
// @author Dominik Müller (\email{dominikmueller64@@yahoo.de})
//
//' @export
// [[Rcpp::export(".LD")]]
Rcpp::List LD(const arma::vec& x,
              const arma::vec& y,
              const double px,
              const double py,
              const bool any_na,
              const bool is_phased,
              const bool r_only)
{
  if (r_only) {
    return Rcpp::List::create(Rcpp::Named("r") = LD_r(x, y, px, py, any_na, is_phased));
    
  } else {
  LD_out t = LD_struct(x, y, px, py, any_na, is_phased);
  return Rcpp::List::create(Rcpp::Named("n") = t.n,
                            Rcpp::Named("D") = t.D,
                            Rcpp::Named("Dprime") = t.Dprime,
                            Rcpp::Named("r") = t.r,
                            Rcpp::Named("r2") = pow(t.r, 2),
                            Rcpp::Named("chi2") = t.chi2,
                            Rcpp::Named("pval") = t.pval);
  }
}


//' @export
// [[Rcpp::export(".LD_mult")]]
DataFrame LD_mult(arma::mat& X,
                  const arma::vec& p,
                  const arma::mat& matr,
                  const bool is_phased,
                  const bool any_na)
{
  // const int N = X.n_rows;
  // const int M = X.n_cols;
  int len = matr.n_rows;
  
  std::map<std::pair<int, int>, int> table;
  // Initialize vectors.
  IntegerVector ix1v(len);
  IntegerVector ix2v(len);
  IntegerVector blockv(len);
  IntegerVector nv(len);
  NumericVector Dv(len);
  NumericVector Dprimev(len);
  NumericVector rv(len);
  NumericVector r2v(len);
  NumericVector chi2v(len);
  NumericVector pvalv(len);
  
  int ct = 0;
  for (int j = 0; j < len; ++j, ++ct) {
    int ix1 = (ix1v(ct) = matr(j, 0)) - 1; // Implicit convertion to int!
    int ix2 = (ix2v(ct) = matr(j, 1)) - 1;
    blockv(ct) = matr(j, 2);
    
    std::map<std::pair<int, int>, int>::iterator it;
    it = table.find(std::make_pair(ix1, ix2));
    
    if (it == table.end()) {
      // Not found, compute and insert.
      // 3.5x speedup! -> Rcpp::List is slow, lots of overhead.
      LD_out ret = LD_struct(X.col(ix1), X.col(ix2), p(ix1), p(ix2), any_na, is_phased);
      nv(ct) = ret.n;
      Dv(ct) = ret.D;
      Dprimev(ct) = ret.Dprime;
      rv(ct) = ret.r;
      r2v(ct) = ret.r2;
      chi2v(ct) = ret.chi2;
      pvalv(ct) = ret.pval;
      table.insert(std::make_pair(std::make_pair(ix1, ix2), ct));
    } else {
      // Found, copy data.
      int ix = it->second;
      nv(ct) = nv(ix);
      Dv(ct) = Dv(ix);
      Dprimev(ct) = Dprimev(ix);
      rv(ct) = rv(ix);
      r2v(ct) = r2v(ix);
      chi2v(ct) = chi2v(ix);
      pvalv(ct) = pvalv(ct);
    }
  }
  return DataFrame::create(Named("i") = ix1v,
                           Named("j") = ix2v,
                           Named("block") = blockv,
                           Named("n") = nv,
                           Named("D") = Dv,
                           Named("Dprime") = Dprimev,
                           Named("r") = rv,
                           Named("r2") = r2v,
                           Named("chi2") = chi2v,
                           Named("pval") = pvalv);
}

// Switch for enabeling caching, as this wastes time if not necessary!
//' @export
// [[Rcpp::export(".LD_mult_r_dev")]]
Rcpp::DataFrame LD_mult_r_dev(arma::mat& X,
                    const arma::vec& p,
                    const arma::mat& matr,
                    const bool is_phased,
                    const bool any_na,
                    const bool cache)
{
  // const int N = X.n_rows;
  // const int M = X.n_cols;
  int len = matr.n_rows;
  
  Rcpp::IntegerVector blockv(len);
  Rcpp::NumericVector rv(len);
  
  if (cache) {
    std::map<std::pair<int, int>, int> table;
    int ct = 0;
    for (int j = 0; j < len; ++j, ++ct) {
      int ix1 = matr(j, 0) - 1; // Implicit convertion to int!
      int ix2 = matr(j, 1) - 1;
      blockv(ct) = matr(j, 2);
      std::map<std::pair<int, int>, int>::iterator it;
      it = table.find(std::make_pair(ix1, ix2));
      if (it == table.end()) {
        rv(ct) = LD_r(X.col(ix1), X.col(ix2), p(ix1), p(ix2), any_na, is_phased);
        table.insert(std::make_pair(std::make_pair(ix1, ix2), ct));
      } else {
        int ix = it->second;
        rv(ct) = rv(ix);
      }
    }
  } // if (cache)
  else {
    int ct = 0;
    for (int j = 0; j < len; ++j, ++ct) {
      int ix1 = matr(j, 0) - 1; // Implicit convertion to int!
      int ix2 = matr(j, 1) - 1;
      blockv(ct) = matr(j, 2);
      rv(ct) = LD_r(X.col(ix1), X.col(ix2), p(ix1), p(ix2), any_na, is_phased);
    }
  }
  
  return Rcpp::DataFrame::create(Rcpp::Named("block") = blockv,   
                                 Rcpp::Named("r") = rv);
}

//' @export
// [[Rcpp::export(".LD_mult_r")]]
Rcpp::DataFrame LD_mult_r(arma::mat& X,
                          const arma::vec& p,
                          const arma::mat& matr,
                          const bool is_phased,
                          const bool any_na)
{
  int len = matr.n_rows;
  
  std::map<std::pair<int, int>, int> table;
  Rcpp::IntegerVector blockv(len);
  Rcpp::NumericVector rv(len);
  
  int ct = 0;
  for (int j = 0; j < len; ++j, ++ct) {
    int ix1 = matr(j, 0) - 1; // Implicit convertion to int!
    int ix2 = matr(j, 1) - 1;
    blockv(ct) = matr(j, 2);
    
    std::map<std::pair<int, int>, int>::iterator it;
    it = table.find(std::make_pair(ix1, ix2));
    
    if (it == table.end()) {
      rv(ct) = LD_r(X.col(ix1), X.col(ix2), p(ix1), p(ix2), any_na, is_phased);
      table.insert(std::make_pair(std::make_pair(ix1, ix2), ct));
    } else {
      int ix = it->second;
      rv(ct) = rv(ix);
    }
  }
  return Rcpp::DataFrame::create(Rcpp::Named("block") = blockv,   
                                 Rcpp::Named("r") = rv);
}


// [[Rcpp::export]]
arma::mat test(const arma::mat& X) {
  int n = X.n_rows;
  int m = X.n_cols;

  mat Xs(X.memptr(), n, m);
  vec pv = conv_to< vec >::from(mean(Xs));
  for (int j = 0; j < m; j++) {
    double p = pv(j);
    double s = sqrt(p * (1 - p));
    for (int i = 0; i < n; i++) {
      // Rcout << Xs(i, j) << std::endl;
      Xs(i, j) = (Xs(i, j) - p) / s;
    }
  }
  // return Xs;

  mat out(m, m);
  for (int i = 0; i < m - 1; i++) {
    for (int j = i + 1; j < m; j++) {
      // mat tmp = arma::cor(X.col(i), X.col(j));
      // out(i, j) = tmp(0, 0);
      // out(j, i) = tmp(0, 0);
      // std::inner_product(Xs.begin_col(i), Xs.end_col(i), Xs.begin_col(j), 0.0);
      double tmp = pow(arma::dot(Xs.col(i), Xs.col(j)) / n, 2);
      out(i, j) = tmp;
      out(j, i) = tmp;
    }
    // mat tmp = arma::cor(X.col(i));
    out(i, i) = pow(arma::dot(Xs.col(i), Xs.col(i)) / n, 2);
  }
    out(m - 1, m - 1) = pow(arma::dot(Xs.col(m - 1), Xs.col(m - 1)) / n, 2);
  // mat tmp = arma::cor(X.col(m - 1));
  // out(m - 1, m - 1) = tmp(0, 0);
  return out;
}

