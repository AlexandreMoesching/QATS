#pragma once
#include <RcppArmadillo.h>
#include <array>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp17")]]

struct par {
  int m;            // Cardinality of the state space
  arma::vec logPi;  // Log-initial distribution
  arma::mat qq;     // Log-transition probabilities
  arma::mat GG;     // Cumulative log-emission densities
};

// Result struct for argH* functions: change points, states, and best score.
// Uses fixed-size stack arrays instead of heap-allocated arma::ivec to avoid
// allocations in the tight inner loops.
struct kih {
  kih() : h_star(R_NegInf) { k_star.fill(-1); i_star.fill(-1); }
  std::array<int, 2> k_star;  // Up to 2 change points
  std::array<int, 3> i_star;  // Up to 3 segment states
  double h_star;
};

struct lrx0 {
  int l;
  int r;
  int x0 = 0;
};

double G0(const arma::vec& xx, int n, const par& par);
double G1(                int i1,                 lrx0 win, const par& par);
double G2(int k1,         int i1, int i2,         lrx0 win, const par& par);
double G3(int k1, int k2, int i1, int i2, int i3, lrx0 win, const par& par);

double H1(                lrx0 win, const par& par);
double H2(int k1,         lrx0 win, const par& par);
double H3(int k1, int k2, lrx0 win, const par& par);

kih argH1(                lrx0 win, const par& par);
kih argH2(int k1,         lrx0 win, const par& par);
kih argH3(int k1, int k2, lrx0 win, const par& par);

void argH1_ref(kih& res,  lrx0 win, const par& par);
void argH2_ref(kih& res,  lrx0 win, const par& par);
void argH3_ref(kih& res,  lrx0 win, const par& par);

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
