#include "simulation.h"
using namespace Rcpp;
using namespace std::chrono;

//' Generate a sample and estimate paths using QATS and different seeds
//'
//' @param n Length of the sequence
//' @param m Cardinality of the state space
//' @param Pi Initial state distribution
//' @param pp Transition matrix
//' @param mu Means of the normal emission distributions
//' @param sigma Standard deviations of the normal emission distributions
//' @param d0 Smallest search interval
//' @param n_seeds Number of seeds for the optimistic search
//' @param rotate Indicates whether or not the gain functions have to be rotated
//' @param n_rep Number of repetitions (for timing)
//' @param n_sim Number of simulations
//'
//' @return A matrix containing estimation times and errors for QATS
//'
//' @export
// [[Rcpp::export(name = "QATS_nseeds_norm")]]
arma::mat QATS_nseeds_norm_cpp(int n, int m,
                               const arma::vec& Pi, const arma::mat& pp,
                               const arma::vec& mu, const arma::vec& sigma,
                               int d0, arma::ivec n_seeds, bool rotate, int n_rep,
                               int n_sim) {
  par0 smp;
  par  params;
  arma::ivec xx(n), zz(n), SS(n);
  int UU;
  auto start = high_resolution_clock::now();
  auto stop  = high_resolution_clock::now();
  duration<double> time = stop - start;

  const int num_n_seeds = n_seeds.n_elem;
  arma::mat res(n_sim, 4 * num_n_seeds);

  for (int j = 0; j < n_sim; j++) {
    smp    = sample_norm_HMM_cpp(n, m, Pi, pp, mu, sigma);
    params = {m, smp.logPi, smp.qq, smp.GG};

    for (int s = 0; s < num_n_seeds; s++) {
      const opts options = {d0, n_seeds[s], rotate};

      UU    = 1;
      start = high_resolution_clock::now();
      for (int i = 0; i < n_rep; i++) {
        QATS_cpp(xx, zz, SS, UU, n, params, options);
      }
      stop = high_resolution_clock::now();
      time = stop - start;

      res.row(j).subvec(s*4, s*4 + 3) = {
        time.count() / n_rep,
        lp_norm_cpp(smp.xx, xx, n, 0),
        lp_norm_cpp(smp.xx, xx, n, 1),
        lp_norm_cpp(smp.xx, xx, n, 2)
      };
    }
  }
  return res;
}

//' Generate a sample and estimate paths using QATS and Viterbi
//'
//' @param n Length of the sequence
//' @param m Cardinality of the state space
//' @param Pi Initial state distribution
//' @param pp Transition matrix
//' @param mu Means of the normal emission distributions
//' @param sigma Standard deviations of the normal emission distributions
//' @param d0 Smallest search interval
//' @param n_seeds Number of seeds for the optimistic search
//' @param rotate Indicates whether or not the gain functions have to be rotated
//' @param n_rep Number of repetitions (for timing)
//' @param n_sim Number of simulations
//'
//' @return A list containing estimation times and errors for both QATS and
//' Viterbi
//'
//' @export
// [[Rcpp::export(name = "QATS_vs_Viterbi_norm")]]
List QATS_vs_Viterbi_norm_cpp(int n, int m,
                              const arma::vec& Pi, const arma::mat& pp,
                              const arma::vec& mu, const arma::vec& sigma,
                              int d0, int n_seeds, bool rotate, int n_rep,
                              int n_sim) {
  const opts options = {d0, n_seeds, rotate};
  par0 smp;
  par  params;
  arma::ivec xx(n), zz(n), SS(n);
  arma::imat zeta(m, n);
  arma::mat  rho(m, n);
  int UU;
  auto start = high_resolution_clock::now();
  auto stop  = high_resolution_clock::now();
  duration<double> time = stop - start;

  arma::mat res_Vit(n_sim, 4), res_QATS(n_sim, 4);

  for (int j = 0; j < n_sim; j++) {
    smp    = sample_norm_HMM_cpp(n, m, Pi, pp, mu, sigma);
    params = {m, smp.logPi, smp.qq, smp.GG};

    // Viterbi
    start = high_resolution_clock::now();
    for (int i = 0; i < n_rep; i++) {
      Viterbi_cpp(xx, zeta, rho, n, m, smp.logPi, smp.qq, smp.g_mseq);
    }
    stop = high_resolution_clock::now();
    time = stop - start;
    res_Vit.row(j) = {time.count() / n_rep,
                      lp_norm_cpp(smp.xx, xx, n, 0),
                      lp_norm_cpp(smp.xx, xx, n, 1),
                      lp_norm_cpp(smp.xx, xx, n, 2)};

    // QATS
    UU    = 1;
    start = high_resolution_clock::now();
    for (int i = 0; i < n_rep; i++) {
      QATS_cpp(xx, zz, SS, UU, n, params, options);
    }
    stop = high_resolution_clock::now();
    time = stop - start;
    res_QATS.row(j) = {time.count() / n_rep,
                       lp_norm_cpp(smp.xx, xx, n, 0),
                       lp_norm_cpp(smp.xx, xx, n, 1),
                       lp_norm_cpp(smp.xx, xx, n, 2)};
  }

  return List::create(Named("res_Vit")  = res_Vit,
                      Named("res_QATS") = res_QATS);
}

//' Generate a sample and estimate paths using QATS, PMAP and Viterbi
//'
//' @param n Length of the sequence
//' @param m Cardinality of the state space
//' @param Pi Initial state distribution
//' @param pp Transition matrix
//' @param mu Means of the normal emission distributions
//' @param sigma Standard deviations of the normal emission distributions
//' @param d0 Smallest search interval
//' @param n_seeds Number of seeds for the optimistic search
//' @param rotate Indicates whether or not the gain functions have to be rotated
//' @param n_rep Number of repetitions (for timing)
//' @param n_sim Number of simulations
//'
//' @return A list containing estimation times and errors for QATS, Viterbi,
//' and PMAP
//'
//' @export
// [[Rcpp::export(name = "compare_norm")]]
List compare_norm_cpp(int n, int m,
                      const arma::vec& Pi, const arma::mat& pp,
                      const arma::vec& mu, const arma::vec& sigma,
                      int d0, int n_seeds, bool rotate, int n_rep,
                      int n_sim) {
  const opts options = {d0, n_seeds, rotate};
  par0 smp;
  par  params;
  arma::ivec xx(n), zz(n), SS(n);
  arma::imat zeta(m, n);
  arma::mat  rho(m, n);
  arma::mat  alpha_hat(m, n), alpha_bar(m, n);
  arma::mat   beta_hat(m, n),  beta_bar(m, n, arma::fill::ones);
  arma::vec  cc_inv(n);
  int UU;
  auto start = high_resolution_clock::now();
  auto stop  = high_resolution_clock::now();
  duration<double> time = stop - start;

  arma::mat res_Vit(n_sim, 4), res_PMAP(n_sim, 4), res_QATS(n_sim, 4);

  for (int j = 0; j < n_sim; j++) {
    smp    = sample_norm_HMM_cpp(n, m, Pi, pp, mu, sigma);
    params = {m, smp.logPi, smp.qq, smp.GG};

    // Viterbi
    start = high_resolution_clock::now();
    for (int i = 0; i < n_rep; i++) {
      Viterbi_cpp(xx, zeta, rho, n, m, smp.logPi, smp.qq, smp.g_mseq);
    }
    stop = high_resolution_clock::now();
    time = stop - start;
    res_Vit.row(j) = {time.count() / n_rep,
                      lp_norm_cpp(smp.xx, xx, n, 0),
                      lp_norm_cpp(smp.xx, xx, n, 1),
                      lp_norm_cpp(smp.xx, xx, n, 2)};

    // PMAP — reset beta_bar to ones before each call
    beta_bar.ones();
    start = high_resolution_clock::now();
    for (int i = 0; i < n_rep; i++) {
      beta_bar.col(n-1).ones();   // backward initialisation
      PMAP_cpp(xx, alpha_hat, alpha_bar, beta_hat, beta_bar, cc_inv,
               n, m, smp.Pi, smp.pp, smp.f_mseq);
    }
    stop = high_resolution_clock::now();
    time = stop - start;
    res_PMAP.row(j) = {time.count() / n_rep,
                       lp_norm_cpp(smp.xx, xx, n, 0),
                       lp_norm_cpp(smp.xx, xx, n, 1),
                       lp_norm_cpp(smp.xx, xx, n, 2)};

    // QATS
    UU    = 1;
    start = high_resolution_clock::now();
    for (int i = 0; i < n_rep; i++) {
      QATS_cpp(xx, zz, SS, UU, n, params, options);
    }
    stop = high_resolution_clock::now();
    time = stop - start;
    res_QATS.row(j) = {time.count() / n_rep,
                       lp_norm_cpp(smp.xx, xx, n, 0),
                       lp_norm_cpp(smp.xx, xx, n, 1),
                       lp_norm_cpp(smp.xx, xx, n, 2)};
  }

  return List::create(Named("res_Vit")  = res_Vit,
                      Named("res_PMAP") = res_PMAP,
                      Named("res_QATS") = res_QATS);
}

double lp_norm_cpp(const arma::ivec& xx_0, const arma::ivec& xx_1, int n, int p) {
  if (p > 0) {
    return pow(sum(pow(abs(arma::conv_to<arma::vec>::from(xx_1 - xx_0)), p)), 1.0/p) / n;
  } else {
    return (double)arma::accu(xx_1 != xx_0) / n;
  }
}
