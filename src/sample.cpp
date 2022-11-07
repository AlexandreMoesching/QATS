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

par0 sample_norm_HMM_cpp(int n, int m, const arma::vec& Pi,
                         const arma::mat& pp, double sigma) {
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

  // Declare and compute xx
  arma::ivec xx(n, arma::fill::zeros);
  xx[0] = SampleFromCumsumProbVec(m, cumsum(Pi).t());
  for (int k = 1; k < n; k++) {
    xx[k] = SampleFromCumsumProbVec(m, pp_cumsum.row(xx[k-1]));
  }

  // At this stage, xx is still on the scale 0:(m-1)

  // Declare and compute yy
  arma::vec yy(n, arma::fill::randn);
  yy = yy * sigma + xx + 1; // + 1 to shift to the scale 1:m

  // Allocate xx and yy
  par.xx = xx;
  par.yy = yy;

  // Declare and compute f_mseq, g_mseq and GG
  arma::mat f_mseq(m, n, arma::fill::zeros);
  arma::mat g_mseq(m, n, arma::fill::zeros);
  arma::mat GG(m, n, arma::fill::zeros);
  for (int i = 0; i < m; i++) {
    f_mseq.row(i) = arma::normpdf(    yy, i + 1.0, sigma).t();
    g_mseq.row(i) = arma::log_normpdf(yy, i + 1.0, sigma).t();
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
                                const arma::mat& pp, double sigma) {
  // Generate an HMM and compute necessary parameters
  par0 par = sample_norm_HMM_cpp(n, m, Pi, pp, sigma);

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
