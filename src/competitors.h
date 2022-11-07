#include <RcppArmadillo.h>
#include <chrono>
using namespace Rcpp;
using namespace std::chrono;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp17")]]

List Viterbi_timer_cpp(double n_rep,
                       int n,
                       int m,
                       const arma::vec& logPi,
                       const arma::mat& qq,
                       const arma::mat& g_mseq);

void Viterbi_cpp(arma::ivec& xx,
                 arma::imat& zeta,
                 arma::mat& rho,
                 const int& n,
                 const int& m,
                 const arma::vec& logPi,
                 const arma::mat& qq,
                 const arma::mat& g_mseq);

arma::vec G_classifier_cpp(double C1,
                           double C2,
                           double C3,
                           double C4,
                           int n,
                           int m,
                           const arma::vec& Pi,
                           const arma::vec& logPi,
                           const arma::mat& pp,
                           const arma::mat& qq,
                           const arma::mat& f_mseq,
                           const arma::mat& g_mseq);

arma::mat K_segmentation_cpp(int K_max,
                             int n,
                             int m,
                             const arma::vec& logPi,
                             const arma::mat& qq,
                             const arma::mat& g_mseq);
