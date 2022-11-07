#include "competitors.h"

//' Viterbi decoder, timing in C++
//'
//' @param n Length of the observation sequence
//' @param n_rep Number of repetition (for timing)
//' @param m Cardinality of the state space
//' @param logPi Initial log-distribution
//' @param qq Log-transition matrix
//' @param g_mseq Log-densities of the emission distributions
//'
//' @return Estimated sequence
//' @keywords internal
// [[Rcpp::export]]
List Viterbi_timer_cpp(double n_rep,
                       int n, int m,
                       const arma::vec& logPi,
                       const arma::mat& qq,
                       const arma::mat& g_mseq) {
  // Preallocation
  arma::ivec xx(n, arma::fill::zeros);
  arma::imat zeta(m, n);
  arma::mat rho(m, n);
  // Get starting time point
  auto start = std::chrono::high_resolution_clock::now();
  // Do stuff
  for (int i = 0; i < n_rep; i++) {
    Viterbi_cpp(xx, zeta, rho, n, m, logPi, qq, g_mseq);
  }
  // Get ending time point
  auto stop = std::chrono::high_resolution_clock::now();
  // Compute difference
  std::chrono::duration<double> time_s = stop - start;
  std::chrono::duration<double, std::micro> time_ms = stop - start;
  // Return
  return List::create(Named("xx") = xx + 1,
                      Named("time") = time_s.count()/n_rep,
                      Named("time_ms") = time_ms.count()/n_rep);
}

//' Viterbi decoder
//'
//' @param xx Pre-declared sequence to be updated
//' @param n Length of the observation sequence
//' @param m Cardinality of the state space
//' @param logPi Initial log-distribution
//' @param qq Log-transition matrix
//' @param g_mseq Log-densities of the emission distributions
//'
//' @return Estimated sequence
//' @keywords internal
// [[Rcpp::export]]
void Viterbi_cpp(arma::ivec& xx, arma::imat& zeta, arma::mat& rho,
                 const int& n, const int& m, const arma::vec& logPi,
                 const arma::mat& qq, const arma::mat& g_mseq) {
  // Declaration
  arma::vec tmp(m);
  // Forward pass
  rho.col(0) = logPi + g_mseq.col(0);
  if (n > 1) {
    for (int k = 1; k < n; k++) {
      for (int i = 0; i < m; i++) {
        tmp = rho.col(k-1) + qq.col(i);
        zeta(i, k-1) = tmp.index_max();
        rho(i, k) = max(tmp) + g_mseq(i, k);
      }
    }
  }
  // Backward pass
  xx[n-1] = rho.col(n-1).index_max();
  if (n > 1) {
    for (int k = n-2; k >= 0; k--) {
      xx[k] = zeta(xx[k+1], k);
    }
  }
}

//' G-classifier
//'
//' @param C1 Constant 1
//' @param C2 Constant 2
//' @param C3 Constant 3
//' @param C4 Constant 4
//' @param n Length of the observation sequence
//' @param m Cardinality of the state space
//' @param Pi Initial distribution
//' @param logPi Initial log-distribution
//' @param pp Transition matrix
//' @param qq Log-transition matrix
//' @param f_mseq Densities of the emission distributions
//' @param g_mseq Log-densities of the emission distributions
//'
//' @return Estimated sequence
//' @keywords internal
// [[Rcpp::export]]
arma::vec G_classifier_cpp(double C1, double C2, double C3, double C4,
                           int n, int m,
                           const arma::vec& Pi, const arma::vec& logPi,
                           const arma::mat& pp, const arma::mat& qq,
                           const arma::mat& f_mseq, const arma::mat& g_mseq) {
  // Declaration
  double C24 = C2 + C4, C234 = C2 + C3 + C4;
  arma::vec prob_x = Pi, tmp(m), xx(n);
  arma::mat logprob_x(m, n), ppT = pp.t(), zeta(m, n), hh(m, n), rho(m, n);
  // Pre-computing some useful quantities
  logprob_x.col(0) = logPi;
  if (n > 1) {
    for (int k = 1; k < n; k++) {
      prob_x = ppT * prob_x;
      logprob_x.col(k) = log(prob_x);
    }
  }
  if (C1 > 0) {
    // Declaration
    arma::mat alpha(m, n), alpha_bar(m, n);
    arma::vec cc(n), beta_bar(m, arma::fill::ones);
    // First pass forward
    tmp = f_mseq.col(0) % Pi;
    cc[0] = sum(tmp);
    alpha.col(0) = tmp/cc[0];
    if (n > 1) {
      for (int k = 1; k < n; k++) {
        for (int i = 0; i < m; i++) {
          tmp[i] = f_mseq(i, k) * sum(pp.col(i) % alpha.col(k-1));
        }
        cc[k] = sum(tmp);
        alpha.col(k) = tmp/cc[k];
      }
    }
    // First pass backward
    alpha_bar.col(n-1) = alpha.col(n-1);
    if (n > 1) {
      for (int k = n-2; k >= 0; k--) {
        for (int i = 0; i < m; i++) {
          tmp[i] = sum((arma::vectorise(pp.row(i)) %
            arma::vectorise(f_mseq.col(k))) %
            beta_bar);
        }
        beta_bar = tmp/cc[k+1];
        alpha_bar.col(k) = beta_bar % alpha.col(k);
      }
    }
    hh         = C1 * log(alpha_bar)        + C2 * g_mseq        + C3   * logprob_x;
    rho.col(0) = C1 * log(alpha_bar.col(0)) + C2 * g_mseq.col(0) + C234 * logprob_x.col(0);
  } else {
    hh         = C2 * g_mseq        + C3   * logprob_x;
    rho.col(0) = C2 * g_mseq.col(0) + C234 * logprob_x.col(0);
  }
  if (n > 1) {
    for (int k = 1; k < n; k++) {
      for (int i = 0; i < m; i++) {
        tmp = rho.col(k-1) + C24 * qq.col(i);
        zeta(i, k-1) = tmp.index_max();
        rho(i, k) = max(tmp) + hh(i, k);
      }
    }
  }
  // Last pass backward
  xx[n-1] = rho.col(n-1).index_max();
  if (n > 1) {
    for (int k = n-2; k >= 0; k--) {
      xx[k] = zeta(xx[k+1], k);
    }
  }
  // Return result
  return xx + 1;
}

//' K-segmentation
//'
//' @param K_max Maximum number of constant pieces
//' @param n Length of the observation sequence
//' @param m Cardinality of the state space
//' @param logPi Initial log-distribution
//' @param qq Log-transition matrix
//' @param g_mseq Log-densities of the emission distributions
//'
//' @return Estimated sequences
//' @keywords internal
// [[Rcpp::export]]
arma::mat K_segmentation_cpp(int K_max, int n, int m,
                             const arma::vec& logPi,
                             const arma::mat& qq,
                             const arma::mat& g_mseq) {
  // Declaration
  int i_max, s_max, t;
  double curr_max, tmp;
  arma::mat xx(K_max, n), ss(K_max, n);
  arma::cube gamma(m, n, K_max, arma::fill::value(R_NegInf));
  arma::cube delta_x(m, n, K_max), delta_s(m, n, K_max);
  // Forward pass
  gamma.slice(0).col(0) = g_mseq.col(0) + logPi;
  if (n > 1) {
    for (int k = 1; k < n; k++) {
      for (int i = 0; i < m; i++) {
        for (int s = 0; s < K_max; s++) {
          curr_max = R_NegInf;
          i_max = m;
          s_max = K_max;
          tmp = R_NegInf;
          for (int j = 0; j < m; j++) {
            if (j == i) {
              t = s;
            } else {
              t = s - 1;
            }
            if (t >= 0) {
              tmp = gamma(j, k-1, t) + qq(j, i);
            } else {
              tmp = R_NegInf;
            }
            if (tmp > curr_max) {
              curr_max = tmp;
              i_max = j;
              s_max = t;
            }
          }
          gamma(i, k, s) = curr_max + g_mseq(i, k);
          delta_x(i, k, s) = i_max;
          delta_s(i, k, s) = s_max;
        }
      }
    }
  }
  // Backward pass
  for (int K = 0; K < K_max; K++) {
    xx(K, n-1) = gamma.slice(K).col(n-1).index_max();
    ss(K, n-1) = K;
  }
  for (int k = n-2; k >= 0; k--) {
    for (int K = 0; K < K_max; K++) {
      xx(K, k) = delta_x(xx(K, k+1), k+1, ss(K, k+1));
      ss(K, k) = delta_s(xx(K, k+1), k+1, ss(K, k+1));
    }
  }
  // Return result
  return xx + 1;
}
