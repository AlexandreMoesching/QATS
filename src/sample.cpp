#include "sample.h"

int SampleFromCumsumProbVec(int m, const arma::rowvec& CumsumProbVec) {
  double u = arma::randu();
  int res = 0;
  for (int i = 0; i < m; i++) {
    if (u > CumsumProbVec[i]) {
      res = i + 1;
    } else {
      break;
    }
  }
  return res;
}

par0 sample_norm_HMM_cpp(int n, int m,
                         const arma::vec& Pi,
                         const arma::mat& pp,
                         const arma::vec& mu,
                         const arma::vec& sigma) {
  par0 par;
  par.n = n;
  par.m = m;
  par.Pi = Pi;
  par.logPi = log(Pi);
  par.pp = pp;
  par.qq = log(pp);

  // Declare row-cumulative version of pp
  arma::mat pp_cumsum(m, m);
  for (int i = 0; i < m; i++) {
    pp_cumsum.row(i) = cumsum(pp.row(i));
  }

  // Declare and compute xx and yy
  arma::ivec xx(n, arma::fill::zeros);
  arma::vec yy(n, arma::fill::randn); // Initially white noise Eps, then mu + sigma * Eps
  xx[0] = SampleFromCumsumProbVec(m, cumsum(Pi).t());
  yy[0] = mu[xx[0]] + yy[0] * sigma[xx[0]];
  for (int k = 1; k < n; k++) {
    xx[k] = SampleFromCumsumProbVec(m, pp_cumsum.row(xx[k-1]));
    yy[k] = mu[xx[k]] + sigma[xx[k]] * yy[k];
  }

  // Allocate xx and yy
  par.xx = xx;
  par.yy = yy;

  // Declare and compute f_mseq, g_mseq and GG
  arma::mat f_mseq(m, n, arma::fill::zeros);
  arma::mat g_mseq(m, n, arma::fill::zeros);
  arma::mat GG(m, n, arma::fill::zeros);
  for (int i = 0; i < m; i++) {
    f_mseq.row(i) = arma::normpdf(    yy, mu[i], sigma[i]).t();
    g_mseq.row(i) = arma::log_normpdf(yy, mu[i], sigma[i]).t();
    GG.row(i)     = cumsum(g_mseq.row(i));
  }

  // Allocate f_mseq, g_mseq and GG
  par.f_mseq = f_mseq;
  par.g_mseq = g_mseq;
  par.GG = GG;

  // Return
  return par;
}

//' @keywords internal
// [[Rcpp::export]]
List sample_norm_HMM_export_cpp(int n, int m, const arma::vec& Pi,
                                const arma::mat& pp,
                                const arma::vec& mu,
                                const arma::vec& sigma) {
  // Generate an HMM and compute necessary parameters
  par0 par = sample_norm_HMM_cpp(n, m, Pi, pp, mu, sigma);

  // Return list
  return List::create(Named("n")      = n,
                      Named("m")      = m,
                      Named("mseq")   = seq_len(m),
                      Named("Pi")     = as<std::vector<double>>(wrap(par.Pi)),
                      Named("logPi")  = as<std::vector<double>>(wrap(par.logPi)),
                      Named("pp")     = par.pp,
                      Named("qq")     = par.qq,
                      Named("xx")     = as<std::vector<double>>(wrap(par.xx + 1)),
                      Named("yy")     = as<std::vector<double>>(wrap(par.yy)),
                      Named("f_mseq") = par.f_mseq,
                      Named("g_mseq") = par.g_mseq,
                      Named("GG")     = par.GG);
}
