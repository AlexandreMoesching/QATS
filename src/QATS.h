#pragma once
#include "log_prob.h"
#include <chrono>

struct opts {
  int d0;
  int n_seeds;
  bool rotate;
};

Rcpp::List QATS_timer_cpp(int d0, int n_seeds, bool rotate, int n_rep,
                    int n, int m,
                    const arma::vec& logPi, const arma::mat& qq,
                    const arma::mat& GG, arma::ivec SS, int UU);

void QATS_cpp(arma::ivec& xx, arma::ivec& zz, arma::ivec& SS, int& UU,
              int n, const par& par, const opts& opts);

void         OSH2_cpp(kih& res,
                      const lrx0& lrx0, const par& par, const opts& opts);
void         OSH3_cpp(kih& res,
                      const lrx0& lrx0, const par& par, const opts& opts);
void      OSH3_k0_cpp(kih& tmp,
                      const lrx0& lrx0, const par& par, const opts& opts);
void  OSH3_k0_dir_cpp(kih& tmp, int tau,
                      const lrx0& lrx0, const par& par, const opts& opts);
void OSH3_k0_diag_cpp(kih& tmp,
                      const lrx0& lrx0, const par& par, const opts& opts);

void rep_each_cpp(arma::ivec& xx, const arma::ivec& zz,
                  const arma::ivec& SS, int UU, int n);
