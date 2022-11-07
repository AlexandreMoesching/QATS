#include "log_prob.h"

double G0(const arma::vec& xx, int n, const par& par) {
  // Right-most index in each part of the partition
  arma::uvec kk_r = arma::find(arma::diff(xx) != 0.0);
  arma::uvec SS_r(kk_r.n_elem + 1);
  std::copy(kk_r.begin(), kk_r.end(), SS_r.begin());
  SS_r[kk_r.n_elem] = n - 1;
  // Left-most index in each part of the partition
  arma::uvec kk_l = kk_r + 1;
  arma::uvec SS_l(kk_l.n_elem + 1);
  std::copy(kk_l.begin(), kk_l.end(), SS_l.begin() + 1);
  SS_l[0] = 0;
  // Useful quantities
  arma::vec qqdiag = par.qq.diag();
  // Initial state
  double tmp = par.logPi.at(xx(0));
  // Transitions from one state to itself
  tmp += sum((SS_r - SS_l) % qqdiag(arma::conv_to<arma::uvec>::from(xx(SS_l))));
  // Transitions from one state to another
  for (int k = 0; k < int(kk_l.n_elem); k++) {
    tmp += par.qq.at(xx(SS_l(k)), xx(SS_r(k+1)));
  }
  // Observations at each state
  for (int k = 0; k < int(kk_l.n_elem) + 1; k++) {
    tmp += par.GG.at(xx(SS_l(k)), SS_r(k));
    if (k > 0) {
      tmp -= par.GG.at(xx(SS_l(k)), SS_r(k-1));
    }
  }
  // Return result
  return tmp;
}

double G1(int i1, const lrx0& lrx0, const par& par) {
  double tmp = (lrx0.r - lrx0.l) * par.qq.at(i1, i1) + par.GG.at(i1, lrx0.r);
  if (lrx0.l > 0) {
    tmp += par.qq.at(lrx0.x0, i1) - par.GG.at(i1, lrx0.l - 1);
  } else {
    tmp += par.logPi.at(i1);
  }
  return tmp;
}

double G2(int k1, int i1, int i2, const lrx0& lrx0, const par& par) {
  double tmp =
    (k1 - lrx0.l -1) * par.qq.at(i1, i1) + par.qq.at(i1, i2) +
    (lrx0.r - k1)    * par.qq.at(i2, i2) +
    par.GG.at(i1, k1 - 1) +
    par.GG.at(i2, lrx0.r) - par.GG.at(i2, k1 - 1);
  if (lrx0.l > 0) {
    tmp += par.qq.at(lrx0.x0, i1) - par.GG.at(i1, lrx0.l - 1);
  } else {
    tmp += par.logPi.at(i1);
  }
  return tmp;
}

double G3(int k1, int k2, int i1, int i2, int i3,
          const lrx0& lrx0, const par& par) {
  double tmp =
    (k1 - lrx0.l -1) * par.qq.at(i1, i1) + par.qq.at(i1, i2) +
    (k2 - k1 - 1)    * par.qq.at(i2, i2) + par.qq.at(i2, i3) +
    (lrx0.r - k2)    * par.qq.at(i3, i3) +
    par.GG.at(i1, k1 - 1) +
    par.GG.at(i2, k2 - 1) - par.GG.at(i2, k1 - 1) +
    par.GG.at(i3, lrx0.r) - par.GG.at(i3, k2 - 1);
  if (lrx0.l > 0) {
    tmp += par.qq.at(lrx0.x0, i1) - par.GG.at(i1, lrx0.l - 1);
  } else {
    tmp += par.logPi.at(i1);
  }
  return tmp;
}

double H1(const lrx0& lrx0, const par& par) {
  double h_star = R_NegInf, tmp = 0.0;
  for (int i1 = 0; i1 < par.m; i1++) {
    tmp = G1(i1, lrx0, par);
    if (tmp > h_star) {
      h_star = tmp;
    }
  }
  return h_star;
}

double H2(int k1, const lrx0& lrx0, const par& par) {
  double h_star = R_NegInf, tmp = 0.0;
  for (int i1 = 0; i1 < par.m; i1++) {
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 != i2) {
        tmp = G2(k1, i1, i2, lrx0, par);
        if (tmp > h_star) {
          h_star = tmp;
        }
      }
    }
  }
  return h_star;
}

double H3(int k1, int k2, const lrx0& lrx0, const par& par) {
  double h_star = R_NegInf, tmp = 0.0;
  for (int i1 = 0; i1 < par.m; i1++) {
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 != i2) {
        for (int i3 = 0; i3 < par.m; i3++) {
          if (i2 != i3) {
            tmp = G3(k1, k2, i1, i2, i3, lrx0, par);
            if (tmp > h_star) {
              h_star = tmp;
            }
          }
        }
      }
    }
  }
  return h_star;
}

kih argH1(const lrx0& lrx0, const par& par) {
  kih res;
  double tmp = 0.0;
  for (int i1 = 0; i1 < par.m; i1++) {
    tmp = G1(i1, lrx0, par);
    if (tmp > res.h_star) {
      res.i_star[0] = i1;
      res.h_star = tmp;
    }
  }
  return res;
}

kih argH2(int k1, const lrx0& lrx0, const par& par) {
  kih res; res.k_star[0] = k1;
  double tmp = 0.0;
  for (int i1 = 0; i1 < par.m; i1++) {
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 != i2) {
        tmp = G2(k1, i1, i2, lrx0, par);
        if (tmp > res.h_star) {
          res.i_star[0] = i1;
          res.i_star[1] = i2;
          res.h_star = tmp;
        }
      }
    }
  }
  return res;
}

kih argH3(int k1, int k2, const lrx0& lrx0, const par& par) {
  kih res; res.k_star = {k1, k2};
  double tmp = 0.0;
  for (int i1 = 0; i1 < par.m; i1++) {
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 != i2) {
        for (int i3 = 0; i3 < par.m; i3++) {
          if (i2 != i3) {
            tmp = G3(k1, k2, i1, i2, i3, lrx0, par);
            if (tmp > res.h_star) {
              res.i_star[0] = i1;
              res.i_star[1] = i2;
              res.i_star[2] = i3;
              res.h_star = tmp;
            }
          }
        }
      }
    }
  }
  return res;
}

void argH1_ref(kih& res, const lrx0& lrx0, const par& par) {
  res.h_star = R_NegInf;
  double tmp = 0.0;
  for (int i1 = 0; i1 < par.m; i1++) {
    tmp = G1(i1, lrx0, par);
    if (tmp > res.h_star) {
      res.i_star[0] = i1;
      res.h_star = tmp;
    }
  }
}

void argH2_ref(kih& res, const lrx0& lrx0, const par& par) {
  res.h_star = R_NegInf;
  double tmp = 0.0;
  for (int i1 = 0; i1 < par.m; i1++) {
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 != i2) {
        tmp = G2(res.k_star[0], i1, i2, lrx0, par);
        if (tmp > res.h_star) {
          res.i_star[0] = i1;
          res.i_star[1] = i2;
          res.h_star = tmp;
        }
      }
    }
  }
}

void argH3_ref(kih& res, const lrx0& lrx0, const par& par) {
  res.h_star = R_NegInf;
  double tmp = 0.0;
  for (int i1 = 0; i1 < par.m; i1++) {
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 != i2) {
        for (int i3 = 0; i3 < par.m; i3++) {
          if (i2 != i3) {
            tmp = G3(res.k_star[0], res.k_star[1], i1, i2, i3, lrx0, par);
            if (tmp > res.h_star) {
              res.i_star[0] = i1;
              res.i_star[1] = i2;
              res.i_star[2] = i3;
              res.h_star = tmp;
            }
          }
        }
      }
    }
  }
}

//' @keywords internal
// [[Rcpp::export]]
double H1_dbl_cpp(int l, int r, int x0, int m, const arma::vec& logPi,
                  const arma::mat& qq, const arma::mat& GG) {
  par par = {m, logPi, qq, GG};
  lrx0 lrx0 = {l, r, x0};
  double res = H1(lrx0, par);
  return res;
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec H2_vec_cpp(int l, int r, int x0, int m, const arma::vec& logPi,
                     const arma::mat& qq, const arma::mat& GG) {
  par par = {m, logPi, qq, GG};
  lrx0 lrx0 = {l, r, x0};
  int d = r - l + 1;
  arma::vec res(d, arma::fill::value(NA_REAL));
  for (int k1 = (l + 1); k1 <= r; k1++) {
    res(k1 - l) = H2(k1, lrx0, par);
  }
  return res;
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat H3_mat_cpp(int l, int r, int x0, int m, const arma::vec& logPi,
                     const arma::mat& qq, const arma::mat& GG) {
  par par = {m, logPi, qq, GG};
  lrx0 lrx0 = {l, r, x0};
  int d = r - l + 1;
  arma::mat res(d, d, arma::fill::value(NA_REAL));
  for (int k1 = (l + 1); k1 < r; k1++) {
    for (int k2 = (k1 + 1); k2 <= r; k2++) {
      res((k1 - l), (k2 - l)) = H3(k1, k2, lrx0, par);
    }
  }
  return res;
}
