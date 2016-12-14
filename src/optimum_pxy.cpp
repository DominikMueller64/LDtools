// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <boost/math/tools/minima.hpp>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::imat get_counts(const arma::vec& x,
                      const arma::vec& y){
  int N = x.n_elem;
  arma::imat counts(3, 3);
  std::fill(counts.begin(), counts.end(), 0);
  
  for (int i = 0; i != N; i++) {
    if (arma::is_finite(x(i)) && arma::is_finite(y(i))) {
      counts(2 - x(i), 2 - y(i)) += 1;
    }
  }
  return(counts);
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
    
    return(-like);
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


