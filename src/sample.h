#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp17")]]

struct par0 {
  int n;              // Length of the sequence
  int m;              // Cardinality of the state space
  arma::vec Pi;       // Initial distribution
  arma::vec logPi;    // Log-initial distribution
  arma::mat pp;       // Transition probabilities
  arma::mat qq;       // Log-transition probabilities
  arma::ivec xx;      // Unobserved hidden states
  arma::vec yy;       // Observed sequence
  arma::mat f_mseq;   // Emission densities
  arma::mat g_mseq;   // Log-emission densities
  arma::mat GG;       // Cumulative log-emission densities
};

int SampleFromCumsumProbVec(int m, const arma::rowvec& CumsumProbVec);

par0 sample_norm_HMM_cpp(       int n, int m, const arma::vec& Pi,
                                const arma::mat& pp, double sigma);
List sample_norm_HMM_export_cpp(int n, int m, const arma::vec& Pi,
                                const arma::mat& pp, double sigma);
