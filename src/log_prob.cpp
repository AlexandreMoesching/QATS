#include "log_prob.h"
using namespace Rcpp;

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
  return tmp;
}

// ── Entry-cost helper: logPi[i] if l==0, else qq[x0,i] - GG[i, l-1] ──────────
// Called once per i1 in every H/argH function; keeps the branch out of inner loops.
inline double entry(int i, lrx0 win, const par& par) {
  return (win.l > 0) ? par.qq.at(win.x0, i) - par.GG.at(i, win.l - 1)
                     : par.logPi.at(i);
}

// ── G1/G2/G3 kept for reference and for H1_dbl_cpp/H2_vec_cpp/H3_mat_cpp ─────

double G1(int i1, lrx0 win, const par& par) {
  return (win.r - win.l) * par.qq.at(i1, i1)
       + par.GG.at(i1, win.r)
       + entry(i1, win, par);
}

double G2(int k1, int i1, int i2, lrx0 win, const par& par) {
  return (k1 - win.l - 1) * par.qq.at(i1, i1) + par.qq.at(i1, i2)
       + (win.r - k1)     * par.qq.at(i2, i2)
       + par.GG.at(i1, k1 - 1)
       + par.GG.at(i2, win.r) - par.GG.at(i2, k1 - 1)
       + entry(i1, win, par);
}

double G3(int k1, int k2, int i1, int i2, int i3, lrx0 win, const par& par) {
  return (k1 - win.l - 1) * par.qq.at(i1, i1) + par.qq.at(i1, i2)
       + (k2 - k1 - 1)    * par.qq.at(i2, i2) + par.qq.at(i2, i3)
       + (win.r - k2)      * par.qq.at(i3, i3)
       + par.GG.at(i1, k1 - 1)
       + par.GG.at(i2, k2 - 1) - par.GG.at(i2, k1 - 1)
       + par.GG.at(i3, win.r)   - par.GG.at(i3, k2 - 1)
       + entry(i1, win, par);
}

// ── H1: single-segment score maximised over states ───────────────────────────

double H1(lrx0 win, const par& par) {
  double h_star = R_NegInf, tmp;
  for (int i1 = 0; i1 < par.m; i1++) {
    tmp = G1(i1, win, par);
    if (tmp > h_star) h_star = tmp;
  }
  return h_star;
}

// ── H2: two-segment score with i1-level terms hoisted out of the i2 loop ─────

double H2(int k1, lrx0 win, const par& par) {
  double h_star = R_NegInf;
  const double c_r_k1 = win.r - k1;
  for (int i1 = 0; i1 < par.m; i1++) {
    const double c1 = (k1 - win.l - 1) * par.qq.at(i1, i1)
                    + par.GG.at(i1, k1 - 1)
                    + entry(i1, win, par);
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 == i2) continue;
      const double tmp = c1 + par.qq.at(i1, i2)
                            + c_r_k1 * par.qq.at(i2, i2)
                            + par.GG.at(i2, win.r) - par.GG.at(i2, k1 - 1);
      if (tmp > h_star) h_star = tmp;
    }
  }
  return h_star;
}

// ── H3: three-segment score with i1- and (i1,i2)-level terms hoisted ─────────

double H3(int k1, int k2, lrx0 win, const par& par) {
  double h_star = R_NegInf;
  const double c_r_k2 = win.r - k2;
  for (int i1 = 0; i1 < par.m; i1++) {
    const double c1 = (k1 - win.l - 1) * par.qq.at(i1, i1)
                    + par.GG.at(i1, k1 - 1)
                    + entry(i1, win, par);
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 == i2) continue;
      const double c2 = c1 + par.qq.at(i1, i2)
                           + (k2 - k1 - 1) * par.qq.at(i2, i2)
                           + par.GG.at(i2, k2 - 1) - par.GG.at(i2, k1 - 1);
      for (int i3 = 0; i3 < par.m; i3++) {
        if (i2 == i3) continue;
        const double tmp = c2 + par.qq.at(i2, i3)
                               + c_r_k2 * par.qq.at(i3, i3)
                               + par.GG.at(i3, win.r) - par.GG.at(i3, k2 - 1);
        if (tmp > h_star) h_star = tmp;
      }
    }
  }
  return h_star;
}

// ── argH* by-value variants (return kih) ──────────────────────────────────────

kih argH1(lrx0 win, const par& par) {
  kih res;
  double tmp;
  for (int i1 = 0; i1 < par.m; i1++) {
    tmp = G1(i1, win, par);
    if (tmp > res.h_star) { res.i_star[0] = i1; res.h_star = tmp; }
  }
  return res;
}

kih argH2(int k1, lrx0 win, const par& par) {
  kih res;
  res.k_star[0] = k1;
  const double c_r_k1 = win.r - k1;
  for (int i1 = 0; i1 < par.m; i1++) {
    const double c1 = (k1 - win.l - 1) * par.qq.at(i1, i1)
                    + par.GG.at(i1, k1 - 1)
                    + entry(i1, win, par);
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 == i2) continue;
      const double tmp = c1 + par.qq.at(i1, i2)
                            + c_r_k1 * par.qq.at(i2, i2)
                            + par.GG.at(i2, win.r) - par.GG.at(i2, k1 - 1);
      if (tmp > res.h_star) {
        res.i_star[0] = i1; res.i_star[1] = i2; res.h_star = tmp;
      }
    }
  }
  return res;
}

kih argH3(int k1, int k2, lrx0 win, const par& par) {
  kih res;
  res.k_star = {k1, k2};
  const double c_r_k2 = win.r - k2;
  for (int i1 = 0; i1 < par.m; i1++) {
    const double c1 = (k1 - win.l - 1) * par.qq.at(i1, i1)
                    + par.GG.at(i1, k1 - 1)
                    + entry(i1, win, par);
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 == i2) continue;
      const double c2 = c1 + par.qq.at(i1, i2)
                           + (k2 - k1 - 1) * par.qq.at(i2, i2)
                           + par.GG.at(i2, k2 - 1) - par.GG.at(i2, k1 - 1);
      for (int i3 = 0; i3 < par.m; i3++) {
        if (i2 == i3) continue;
        const double tmp = c2 + par.qq.at(i2, i3)
                               + c_r_k2 * par.qq.at(i3, i3)
                               + par.GG.at(i3, win.r) - par.GG.at(i3, k2 - 1);
        if (tmp > res.h_star) {
          res.i_star[0] = i1; res.i_star[1] = i2; res.i_star[2] = i3;
          res.h_star = tmp;
        }
      }
    }
  }
  return res;
}

// ── argH*_ref in-place variants (reuse pre-allocated kih) ────────────────────

void argH1_ref(kih& res, lrx0 win, const par& par) {
  res.h_star = R_NegInf;
  double tmp;
  for (int i1 = 0; i1 < par.m; i1++) {
    tmp = G1(i1, win, par);
    if (tmp > res.h_star) { res.i_star[0] = i1; res.h_star = tmp; }
  }
}

void argH2_ref(kih& res, lrx0 win, const par& par) {
  res.h_star = R_NegInf;
  const int k1 = res.k_star[0];
  const double c_r_k1 = win.r - k1;
  for (int i1 = 0; i1 < par.m; i1++) {
    const double c1 = (k1 - win.l - 1) * par.qq.at(i1, i1)
                    + par.GG.at(i1, k1 - 1)
                    + entry(i1, win, par);
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 == i2) continue;
      const double tmp = c1 + par.qq.at(i1, i2)
                            + c_r_k1 * par.qq.at(i2, i2)
                            + par.GG.at(i2, win.r) - par.GG.at(i2, k1 - 1);
      if (tmp > res.h_star) {
        res.i_star[0] = i1; res.i_star[1] = i2; res.h_star = tmp;
      }
    }
  }
}

void argH3_ref(kih& res, lrx0 win, const par& par) {
  res.h_star = R_NegInf;
  const int k1 = res.k_star[0], k2 = res.k_star[1];
  const double c_r_k2 = win.r - k2;
  for (int i1 = 0; i1 < par.m; i1++) {
    const double c1 = (k1 - win.l - 1) * par.qq.at(i1, i1)
                    + par.GG.at(i1, k1 - 1)
                    + entry(i1, win, par);
    for (int i2 = 0; i2 < par.m; i2++) {
      if (i1 == i2) continue;
      const double c2 = c1 + par.qq.at(i1, i2)
                           + (k2 - k1 - 1) * par.qq.at(i2, i2)
                           + par.GG.at(i2, k2 - 1) - par.GG.at(i2, k1 - 1);
      for (int i3 = 0; i3 < par.m; i3++) {
        if (i2 == i3) continue;
        const double tmp = c2 + par.qq.at(i2, i3)
                               + c_r_k2 * par.qq.at(i3, i3)
                               + par.GG.at(i3, win.r) - par.GG.at(i3, k2 - 1);
        if (tmp > res.h_star) {
          res.i_star[0] = i1; res.i_star[1] = i2; res.i_star[2] = i3;
          res.h_star = tmp;
        }
      }
    }
  }
}

// ── R-facing exports ──────────────────────────────────────────────────────────

//' @keywords internal
// [[Rcpp::export]]
double H1_dbl_cpp(int l, int r, int x0, int m, const arma::vec& logPi,
                  const arma::mat& qq, const arma::mat& GG) {
  return H1({l, r, x0}, {m, logPi, qq, GG});
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec H2_vec_cpp(int l, int r, int x0, int m, const arma::vec& logPi,
                     const arma::mat& qq, const arma::mat& GG) {
  const par   params = {m, logPi, qq, GG};
  const lrx0  win    = {l, r, x0};
  const int   d      = r - l + 1;
  arma::vec   res(d, arma::fill::value(NA_REAL));
  for (int k1 = l + 1; k1 <= r; k1++) res(k1 - l) = H2(k1, win, params);
  return res;
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat H3_mat_cpp(int l, int r, int x0, int m, const arma::vec& logPi,
                     const arma::mat& qq, const arma::mat& GG) {
  const par   params = {m, logPi, qq, GG};
  const lrx0  win    = {l, r, x0};
  const int   d      = r - l + 1;
  arma::mat   res(d, d, arma::fill::value(NA_REAL));
  for (int k1 = l + 1; k1 < r; k1++)
    for (int k2 = k1 + 1; k2 <= r; k2++)
      res(k1 - l, k2 - l) = H3(k1, k2, win, params);
  return res;
}
