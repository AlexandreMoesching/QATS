#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp17")]]

struct par {
  int m;            // Cardinality of the state space
  arma::vec logPi;  // Log-initial distribution
  arma::mat qq;     // Log-transition probabilities
  arma::mat GG;     // Cumulative log-emission densities
};

struct kih {
  kih() : k_star(2), i_star(3) {}
  arma::ivec k_star;
  arma::ivec i_star;
  double h_star = R_NegInf;
};

struct lrx0 {
  int l;
  int r;
  int x0;
};

double G0(const arma::vec& xx, int n, const par& par);
double G1(                int i1,                 const lrx0& lrx0,
                                                  const par& par);
double G2(int k1,         int i1, int i2,         const lrx0& lrx0,
                                                  const par& par);
double G3(int k1, int k2, int i1, int i2, int i3, const lrx0& lrx0,
                                                  const par& par);

double H1(                const lrx0& lrx0, const par& par);
double H2(int k1,         const lrx0& lrx0, const par& par);
double H3(int k1, int k2, const lrx0& lrx0, const par& par);

kih argH1(                const lrx0& lrx0, const par& par);
kih argH2(int k1,         const lrx0& lrx0, const par& par);
kih argH3(int k1, int k2, const lrx0& lrx0, const par& par);

void argH1_ref(kih& res,  const lrx0& lrx0, const par& par);
void argH2_ref(kih& res,  const lrx0& lrx0, const par& par);
void argH3_ref(kih& res,  const lrx0& lrx0, const par& par);

double    H1_dbl_cpp(int l, int r, int x0, int m,
                     const arma::vec& logPi,
                     const arma::mat& qq,
                     const arma::mat& GG);
arma::vec H2_vec_cpp(int l, int r, int x0, int m,
                     const arma::vec& logPi,
                     const arma::mat& qq,
                     const arma::mat& GG);
arma::mat H3_mat_cpp(int l, int r, int x0, int m,
                     const arma::vec& logPi,
                     const arma::mat& qq,
                     const arma::mat& GG);
