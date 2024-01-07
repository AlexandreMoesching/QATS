#include "sample.h"
#include "QATS.h"
#include "competitors.h"

arma::mat QATS_nseeds_norm_cpp(int n, int m,
                               const arma::vec& Pi, const arma::mat& pp,
                               double sigma,
                               int d0, arma::ivec n_seeds, bool rotate, int n_rep,
                               int n_sim);
List QATS_vs_Viterbi_norm_cpp(int n, int m,
                              const arma::vec& Pi, const arma::mat& pp,
                              double sigma,
                              int d0, int n_seeds, bool rotate, int n_rep,
                              int n_sim);
List compare_norm_cpp(int n, int m,
                      const arma::vec& Pi, const arma::mat& pp,
                      double sigma,
                      int d0, int n_seeds, bool rotate, int n_rep,
                      int n_sim);
double lp_norm_cpp(arma::ivec& xx_0, arma::ivec& xx_1, int n, int p);
