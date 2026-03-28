#include "competitors.h"
using namespace Rcpp;
using namespace std::chrono;

//' Viterbi decoder, timing in C++
//'
//' @param n_rep Number of repetitions (for timing)
//' @param n Length of the observation sequence
//' @param m Cardinality of the state space
//' @param logPi Initial log-distribution
//' @param qq Log-transition matrix
//' @param g_mseq Log-densities of the emission distributions
//'
//' @return Estimated sequence
//' @keywords internal
// [[Rcpp::export]]
List Viterbi_timer_cpp(int n_rep,
                       int n, int m,
                       const arma::vec& logPi,
                       const arma::mat& qq,
                       const arma::mat& g_mseq) {
  arma::ivec xx(n, arma::fill::zeros);
  arma::imat zeta(m, n);
  arma::mat  rho(m, n);
  auto start = high_resolution_clock::now();
  for (int i = 0; i < n_rep; i++) {
    Viterbi_cpp(xx, zeta, rho, n, m, logPi, qq, g_mseq);
  }
  auto stop = high_resolution_clock::now();
  duration<double>             time_s  = stop - start;
  duration<double, std::micro> time_ms = stop - start;
  return List::create(Named("xx")      = xx + 1,
                      Named("time")    = time_s.count()  / n_rep,
                      Named("time_ms") = time_ms.count() / n_rep);
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
void Viterbi_cpp(arma::ivec& xx,
                 arma::imat& zeta,
                 arma::mat& rho,
                 int n, int m,
                 const arma::vec& logPi,
                 const arma::mat& qq,
                 const arma::mat& g_mseq) {
  // Forward pass
  for (int j = 0; j < m; j++) rho(j, 0) = logPi(j) + g_mseq(j, 0);
  if (n > 1) {
    for (int k = 1; k < n; k++) {
      for (int i = 0; i < m; i++) {
        // Single pass: find argmax_j and max_j of rho(j, k-1) + qq(j, i)
        int    best_j   = 0;
        double best_val = rho(0, k - 1) + qq(0, i);
        for (int j = 1; j < m; j++) {
          const double val = rho(j, k - 1) + qq(j, i);
          if (val > best_val) { best_val = val; best_j = j; }
        }
        zeta(i, k - 1) = best_j;
        rho(i, k)      = best_val + g_mseq(i, k);
      }
    }
  }
  // Backward pass
  xx[n - 1] = rho.col(n - 1).index_max();
  if (n > 1) {
    for (int k = n - 2; k >= 0; k--) {
      xx[k] = zeta(xx[k + 1], k);
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
  double C24 = C2 + C4, C234 = C2 + C3 + C4;
  arma::vec prob_x = Pi, tmp(m), xx(n);
  arma::mat logprob_x(m, n), ppT = pp.t(), zeta(m, n), hh(m, n), rho(m, n);
  // Marginal log-probabilities of each state at each time
  logprob_x.col(0) = logPi;
  if (n > 1) {
    for (int k = 1; k < n; k++) {
      prob_x = ppT * prob_x;
      logprob_x.col(k) = log(prob_x);
    }
  }
  if (C1 > 0) {
    arma::mat alpha(m, n), alpha_bar(m, n);
    arma::vec cc(n), beta_bar(m, arma::fill::ones);
    // Forward pass
    for (int j = 0; j < m; j++) tmp[j] = f_mseq(j, 0) * Pi(j);
    cc[0] = 0.0; for (int j = 0; j < m; j++) cc[0] += tmp[j];
    for (int j = 0; j < m; j++) alpha(j, 0) = tmp[j] / cc[0];
    if (n > 1) {
      for (int k = 1; k < n; k++) {
        for (int i = 0; i < m; i++) {
          double s = 0.0;
          for (int j = 0; j < m; j++) s += pp(j, i) * alpha(j, k - 1);
          tmp[i] = f_mseq(i, k) * s;
        }
        cc[k] = 0.0; for (int j = 0; j < m; j++) cc[k] += tmp[j];
        for (int j = 0; j < m; j++) alpha(j, k) = tmp[j] / cc[k];
      }
    }
    // Backward pass
    for (int j = 0; j < m; j++) alpha_bar(j, n - 1) = alpha(j, n - 1);
    if (n > 1) {
      for (int k = n - 2; k >= 0; k--) {
        for (int i = 0; i < m; i++) {
          double s = 0.0;
          for (int j = 0; j < m; j++) s += pp(i, j) * f_mseq(j, k) * beta_bar(j);
          tmp[i] = s;
        }
        for (int j = 0; j < m; j++) beta_bar(j) = tmp[j] / cc[k + 1];
        for (int j = 0; j < m; j++) alpha_bar(j, k) = beta_bar(j) * alpha(j, k);
      }
    }
    hh         = C1 * log(alpha_bar)        + C2 * g_mseq        + C3   * logprob_x;
    rho.col(0) = C1 * log(alpha_bar.col(0)) + C2 * g_mseq.col(0) + C234 * logprob_x.col(0);
  } else {
    hh         = C2 * g_mseq        + C3   * logprob_x;
    rho.col(0) = C2 * g_mseq.col(0) + C234 * logprob_x.col(0);
  }
  // Last forward pass with combined score
  if (n > 1) {
    for (int k = 1; k < n; k++) {
      for (int i = 0; i < m; i++) {
        int    best_j   = 0;
        double best_val = rho(0, k - 1) + C24 * qq(0, i);
        for (int j = 1; j < m; j++) {
          const double val = rho(j, k - 1) + C24 * qq(j, i);
          if (val > best_val) { best_val = val; best_j = j; }
        }
        zeta(i, k - 1) = best_j;
        rho(i, k)      = best_val + hh(i, k);
      }
    }
  }
  // Backward pass
  xx[n - 1] = rho.col(n - 1).index_max();
  if (n > 1) {
    for (int k = n - 2; k >= 0; k--) {
      xx[k] = zeta(xx[k + 1], k);
    }
  }
  return xx + 1;
}

//' PMAP-classifier, timing in C++
//'
//' @param n_rep Number of repetitions (for timing)
//' @param n Length of the observation sequence
//' @param m Cardinality of the state space
//' @param Pi Initial distribution
//' @param pp Transition matrix
//' @param f_mseq Densities of the emission distributions
//'
//' @return Estimated sequence
//' @keywords internal
// [[Rcpp::export]]
List PMAP_timer_cpp(int n_rep,
                    int n, int m,
                    const arma::vec& Pi,
                    const arma::mat& pp,
                    const arma::mat& f_mseq) {
  arma::ivec xx(n, arma::fill::zeros);
  arma::mat alpha_hat(m, n), alpha_bar(m, n);
  arma::mat  beta_hat(m, n),  beta_bar(m, n, arma::fill::ones);
  arma::vec cc_inv(n);
  auto start = high_resolution_clock::now();
  for (int i = 0; i < n_rep; i++) {
    PMAP_cpp(xx, alpha_hat, alpha_bar, beta_hat, beta_bar, cc_inv,
             n, m, Pi, pp, f_mseq);
  }
  auto stop = high_resolution_clock::now();
  duration<double>             time_s  = stop - start;
  duration<double, std::micro> time_ms = stop - start;
  return List::create(Named("xx")      = xx + 1,
                      Named("time")    = time_s.count()  / n_rep,
                      Named("time_ms") = time_ms.count() / n_rep);
}

//' PMAP-classifier
//'
//' @param n Length of the observation sequence
//' @param m Cardinality of the state space
//' @param Pi Initial distribution
//' @param pp Transition matrix
//' @param f_mseq Densities of the emission distributions
//'
//' @return Estimated sequence
//' @keywords internal
// [[Rcpp::export]]
void PMAP_cpp(arma::ivec& xx,
              arma::mat& alpha_hat,
              arma::mat& alpha_bar,
              arma::mat& beta_hat,
              arma::mat& beta_bar,
              arma::vec& cc_inv,
              int n,
              int m,
              const arma::vec& Pi,
              const arma::mat& pp,
              const arma::mat& f_mseq) {
  // Forward pass (scaled)
  double c;
  for (int j = 0; j < m; j++) alpha_bar(j, 0) = f_mseq(j, 0) * Pi(j);
  c = 0.0; for (int j = 0; j < m; j++) c += alpha_bar(j, 0);
  cc_inv[0] = c;
  for (int j = 0; j < m; j++) alpha_hat(j, 0) = alpha_bar(j, 0) / c;
  if (n > 1) {
    for (int k = 1; k < n; k++) {
      c = 0.0;
      for (int i = 0; i < m; i++) {
        double s = 0.0;
        for (int j = 0; j < m; j++) s += pp(j, i) * alpha_hat(j, k - 1);
        alpha_bar(i, k) = f_mseq(i, k) * s;
        c += alpha_bar(i, k);
      }
      cc_inv[k] = c;
      for (int i = 0; i < m; i++) alpha_hat(i, k) = alpha_bar(i, k) / c;
    }
  }
  // Backward pass (scaled); beta_bar[:, n-1] = 1 on entry
  for (int j = 0; j < m; j++) beta_hat(j, n - 1) = beta_bar(j, n - 1) / cc_inv[n - 1];
  if (n > 1) {
    for (int k = n - 2; k >= 0; k--) {
      for (int i = 0; i < m; i++) {
        double s = 0.0;
        for (int j = 0; j < m; j++) s += pp(i, j) * f_mseq(j, k + 1) * beta_hat(j, k + 1);
        beta_bar(i, k) = s;
      }
      for (int i = 0; i < m; i++) beta_hat(i, k) = beta_bar(i, k) / cc_inv[k];
    }
  }
  // Pointwise MAP (single-pass argmax, no temporary vector)
  for (int k = 0; k < n; k++) {
    int    best_j   = 0;
    double best_val = alpha_hat(0, k) * beta_hat(0, k);
    for (int j = 1; j < m; j++) {
      const double val = alpha_hat(j, k) * beta_hat(j, k);
      if (val > best_val) { best_val = val; best_j = j; }
    }
    xx[k] = best_j;
  }
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
          curr_max = R_NegInf; i_max = m; s_max = K_max;
          for (int j = 0; j < m; j++) {
            t = (j == i) ? s : s - 1;
            tmp = (t >= 0) ? gamma(j, k-1, t) + qq(j, i) : R_NegInf;
            if (tmp > curr_max) { curr_max = tmp; i_max = j; s_max = t; }
          }
          gamma(i, k, s)   = curr_max + g_mseq(i, k);
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
  return xx + 1;
}
