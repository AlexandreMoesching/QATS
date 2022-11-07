#include "simulation.h"

//' Generate a sample and estimate paths using QATS and different seeds
//'
//' @param n Length of the sequence
//' @param m Cardinality of the state space
//' @param Pi Initial state distribution
//' @param pp Transition matrix
//' @param sigma Standard deviation of the normal emission distribution
//' @param d0 Smallest search interval
//' @param n_seeds Number of seeds for the optimistic search
//' @param rotate Indicates whether or not the gain functions have to be rotated
//' @param n_rep Number of repetition (for timing)
//' @param n_sim Number of simulations
//'
//' @return A matrix containing estimation times and errors for QATS
//'
//' @export
// [[Rcpp::export(name = "QATS_nseeds_norm")]]
arma::mat QATS_nseeds_norm_cpp(int n, int m,
                               const arma::vec& Pi, const arma::mat& pp,
                               double sigma,
                               int d0, arma::ivec n_seeds, bool rotate, int n_rep,
                               int n_sim) {
  // Declare sample variables
  par0 par0;
  par par;

  // Declare estimation variables
  arma::ivec xx(n), zz(n), SS(n);
  int UU;

  // Declare timers
  auto start = std::chrono::high_resolution_clock::now();
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time = stop - start;

  // Declare output variables
  int num_n_seeds = n_seeds.n_elem;
  arma::mat res(n_sim, 4 * num_n_seeds);

  // Loop over the number of simulations
  for (int j = 0; j < n_sim; j++) {
    // Generate a sample
    par0 = sample_norm_HMM_cpp(n, m, Pi, pp, sigma);
    par = {m, par0.logPi, par0.qq, par0.GG};

    // Loop on the number of seeds
    for (int s = 0; s < num_n_seeds; s++) {
      // Declare options
      const opts opts = {d0, n_seeds[s], rotate};

      // Estimate with QATS
      UU = 1;
      start = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < n_rep; i++) {
        QATS_cpp(xx, zz, SS, UU, n, par, opts);
      }
      stop = std::chrono::high_resolution_clock::now();
      time = stop - start;

      // Return results
      res.row(j).subvec(s*4, s*4 + 3)
        = {time.count()/n_rep,
           lp_norm_cpp(par0.xx, xx, n, 0),
           lp_norm_cpp(par0.xx, xx, n, 1),
           lp_norm_cpp(par0.xx, xx, n, 2)};
    }
  }

  // Return
  return res;
}

//' Generate a sample and estimate paths using QATS and Viterbi
//'
//' @param n Length of the sequence
//' @param m Cardinality of the state space
//' @param Pi Initial state distribution
//' @param pp Transition matrix
//' @param sigma Standard deviation of the normal emission distribution
//' @param d0 Smallest search interval
//' @param n_seeds Number of seeds for the optimistic search
//' @param rotate Indicates whether or not the gain functions have to be rotated
//' @param n_rep Number of repetition (for timing)
//' @param n_sim Number of simulations
//'
//' @return A list containing estimation times and errors for both QATS and
//' Viterbi
//'
//' @export
// [[Rcpp::export(name = "QATS_vs_Viterbi_norm")]]
List QATS_vs_Viterbi_norm_cpp(int n, int m,
                              const arma::vec& Pi, const arma::mat& pp,
                              double sigma,
                              int d0, int n_seeds, bool rotate, int n_rep,
                              int n_sim) {
  // Declare options
  const opts opts = {d0, n_seeds, rotate};

  // Declare sample variables
  par0 par0;
  par par;

  // Declare estimation variables
  arma::ivec xx(n), zz(n), SS(n);
  arma::imat zeta(m, n);
  arma::mat rho(m, n);
  int UU;

  // Declare timers
  auto start = std::chrono::high_resolution_clock::now();
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time = stop - start;

  // Declare output variables
  arma::mat res_Vit(n_sim, 4), res_QATS(n_sim, 4);

  // Loop over the number of simulations
  for (int j = 0; j < n_sim; j++) {
    // Generate a sample
    par0 = sample_norm_HMM_cpp(n, m, Pi, pp, sigma);
    par = {m, par0.logPi, par0.qq, par0.GG};

    // Estimate with Viterbi
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_rep; i++) {
      Viterbi_cpp(xx, zeta, rho, n, m, par0.logPi, par0.qq, par0.g_mseq);
    }
    stop = std::chrono::high_resolution_clock::now();
    time = stop - start;
    res_Vit.row(j) = {time.count()/n_rep,
                      lp_norm_cpp(par0.xx, xx, n, 0),
                      lp_norm_cpp(par0.xx, xx, n, 1),
                      lp_norm_cpp(par0.xx, xx, n, 2)};

    // Estimate with QATS
    UU = 1;
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_rep; i++) {
      QATS_cpp(xx, zz, SS, UU, n, par, opts);
    }
    stop = std::chrono::high_resolution_clock::now();
    time = stop - start;
    res_QATS.row(j) = {time.count()/n_rep,
                       lp_norm_cpp(par0.xx, xx, n, 0),
                       lp_norm_cpp(par0.xx, xx, n, 1),
                       lp_norm_cpp(par0.xx, xx, n, 2)};
  }

  // Return
  return List::create(Named("res_Vit") = res_Vit,
                      Named("res_QATS") = res_QATS);
}

double lp_norm_cpp(arma::ivec& xx_0, arma::ivec& xx_1, int n, int p) {
  double res;
  if (p > 0) {
    res = pow(sum(pow(abs(xx_1 - xx_0), p)), 1.0/p) / n;
  } else {
    arma::uvec tmp = find(xx_1 != xx_0);
    res = tmp.n_elem * 1.0 / n;
  }
  return res;
}
