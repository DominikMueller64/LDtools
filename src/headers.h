int index_greater(const arma::vec& x, const double y);
int index_geq(const arma::vec& x, const double y);

int find_closest(const arma::vec& x, const double y);

double optimum_pxy(const arma::imat& n3x3, 
                   const double px, 
                   const double py,
                   const double lower,
                   const double upper);
  


arma::imat get_counts(const arma::vec& x,
                      const arma::vec& y);