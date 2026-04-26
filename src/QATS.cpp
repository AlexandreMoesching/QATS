#include "QATS.h"
using namespace Rcpp;
using namespace std::chrono;

//' Quick Adaptive Ternary Segmentation, timing in C++
//'
//' @param d0 Smallest search interval
//' @param n_seeds Number of seeds for the optimistic search
//' @param rotate Indicates whether or not the gain functions have to be rotated
//' @param n_rep Number of repetition (for timing)
//' @param n Length of the sequence
//' @param m Cardinality of the state space
//' @param logPi Log-initial distribution
//' @param qq Log-transition probabilities
//' @param GG Cumulative log-emission densities
//' @param SS Initial partition
//' @param UU Length of the initial partition
//'
//' @keywords internal
// [[Rcpp::export]]
List QATS_timer_cpp(int d0, int n_seeds, bool rotate, int n_rep,
                    int n, int m,
                    const arma::vec& logPi, const arma::mat& qq,
                    const arma::mat& GG, arma::ivec SS, int UU) {
  // Preallocation
  arma::ivec xx(n), zz(n);
  SS.resize(n);
  const par    params  = {m, logPi, qq, GG};
  const opts   options = {d0, n_seeds, rotate};
  // Get starting time point
  auto start = high_resolution_clock::now();
  // Do stuff
  for (int i = 0; i < n_rep; i++) {
    QATS_cpp(xx, zz, SS, UU, n, params, options);
  }
  // Get ending time point
  auto stop = high_resolution_clock::now();
  // Compute difference
  duration<double>              time_s  = stop - start;
  duration<double, std::micro>  time_ms = stop - start;
  // Return
  return List::create(Named("xx")      = xx + 1,
                      Named("time")    = time_s.count()  / n_rep,
                      Named("time_ms") = time_ms.count() / n_rep);
}

void QATS_cpp(arma::ivec& xx, arma::ivec& zz, arma::ivec& SS, int& UU,
              int n, const par& par, const opts& opts) {
  // Declaration
  int u = 0, du, k_star_n = 1;
  kih tmp_1, tmp_2, tmp_3, tmp_star;
  lrx0 win;
  // While-loop
  while (u < UU) {
    // Set window parameters
    win.l = SS[u];
    win.r = ((u + 1) < UU) ? SS[u + 1] - 1 : n - 1;
    du = win.r - win.l + 1;
    if (u > 0) win.x0 = zz[u - 1];
    // Look for a new change point
    if (du == 1) {
      tmp_star.h_star = R_NegInf;
    } else if (du == 2) {
      tmp_star.k_star[0] = win.r;
      argH2_ref(tmp_star, win, par);
      k_star_n = 1;
    } else { // du >= 3
      OSH2_cpp(tmp_2, win, par, opts);
      OSH3_cpp(tmp_3, win, par, opts);
      if (tmp_2.h_star > tmp_3.h_star) {
        tmp_star = tmp_2;
        k_star_n = 1;
      } else {
        tmp_star = tmp_3;
        k_star_n = 2;
      }
    }
    // Best constant path for comparison
    argH1_ref(tmp_1, win, par);
    // Update
    if (tmp_star.h_star > tmp_1.h_star + 1e-7) {
      // Shift existing entries right to make room
      if (u < UU-1) {
        for (int j = 1; j < (UU - u); j++) {
          SS[UU - j + k_star_n] = SS[UU - j];
          zz[UU - j + k_star_n] = zz[UU - j];
        }
      }
      // Insert new change points and states
      for (int j = 0; j < k_star_n; j++) SS[u + 1 + j] = tmp_star.k_star[j];
      for (int j = 0; j <= k_star_n; j++) zz[u + j]    = tmp_star.i_star[j];
      // Update UU
      UU += k_star_n;
    } else {
      zz[u] = tmp_1.i_star[0];
      u++;
    }
  }
  // Reconstruct the full sequence from partition
  rep_each_cpp(xx, zz, SS, UU, n);
}

void OSH2_cpp(kih& res, const lrx0& lrx0, const par& par, const opts& opts) {
  int L, R, M, W;
  const double nu = 0.5;
  double HM, HW, a = 0.0;
  kih tmp;
  res.h_star = R_NegInf;
  // Initialization
  L = lrx0.l + 1;
  R = lrx0.r;
  M = (int)floor((L + nu * R) / (1 + nu));
  if (opts.rotate && lrx0.l + 1 < lrx0.r) {
    a = (H2(lrx0.r, lrx0, par) - H2(lrx0.l + 1, lrx0, par))
        / (lrx0.r - lrx0.l - 1);
  }
  HM = H2(M, lrx0, par) - a * M;
  // Narrow the search interval
  while (R - L > opts.d0) {
    bool right_larger = (R - M > M - L);
    W = right_larger ? (int)ceil(R - nu * (R - M))
                     : (int)ceil(L + nu * (M - L));
    if (W == L || W == R) break;
    HW = H2(W, lrx0, par) - a * W;
    if (HW > HM) {
      if (right_larger) L = M; else R = M;
      M = W; HM = HW;
    } else {
      if (right_larger) R = W; else L = W;
    }
  }
  // Full search on the final narrow interval
  for (int k1 = L; k1 <= R; k1++) {
    tmp = argH2(k1, lrx0, par);
    if (tmp.h_star > res.h_star) res = tmp;
  }
}

void OSH3_cpp(kih& res, const lrx0& lrx0, const par& par, const opts& opts) {
  kih tmp;
  res.h_star = R_NegInf;
  if (lrx0.r - lrx0.l > opts.d0) {
    for (int i = 1; i <= opts.n_seeds; i++) {
      tmp.k_star[0] = -1;
      tmp.k_star[1] = (int)floor(lrx0.l + 2
        + (i * 1.0) / (opts.n_seeds + 1) * (lrx0.r - lrx0.l - 1.0));
      OSH3_k0_cpp(tmp, lrx0, par, opts);
      if (tmp.h_star > res.h_star) res = tmp;
    }
  } else {
    // Full grid search for small intervals
    for (int k1 = (lrx0.l + 1); k1 < lrx0.r; k1++) {
      for (int k2 = (k1 + 1); k2 <= lrx0.r; k2++) {
        tmp = argH3(k1, k2, lrx0, par);
        if (tmp.h_star > res.h_star) res = tmp;
      }
    }
  }
}

void OSH3_k0_cpp(kih& tmp, const lrx0& lrx0, const par& par, const opts& opts) {
  double h_new = R_NegInf, h_old = R_NegInf;
  int tau = 0;
  int v = 1, v0 = 20;
  // Alternating horizontal and vertical coordinate searches
  while (((h_old < h_new) && (v < v0)) || (v == 1)) {
    h_old = h_new;
    OSH3_k0_dir_cpp(tmp, tau, lrx0, par, opts);
    // If on the diagonal, also optimise k -> H3(k, k+1)
    if (tmp.k_star[0] + 1 == tmp.k_star[1]) {
      OSH3_k0_diag_cpp(tmp, lrx0, par, opts);
    }
    h_new = tmp.h_star;
    tau   = 1 - tau;
    v++;
  }
  // Compute the final argmax for the selected change-point pair
  argH3_ref(tmp, lrx0, par);
}

void OSH3_k0_dir_cpp(kih& tmp, int tau,
                     const lrx0& lrx0, const par& par, const opts& opts) {
  int L, R, M, W;
  const double nu = 0.5;
  double HM, HW, h_tmp, a = 0.0;
  tmp.h_star = R_NegInf;

  if (tau == 0) {
    // Search over k1 with k2 fixed
    L = lrx0.l + 1;
    R = tmp.k_star[1] - 1;
    M = (tmp.k_star[0] == -1) ? (int)floor((L + nu * R) / (1 + nu))
                               : tmp.k_star[0];
    if (opts.rotate && L < R) {
      a = (H3(R, tmp.k_star[1], lrx0, par) - H3(L, tmp.k_star[1], lrx0, par))
          / (R - L);
    }
    HM = H3(M, tmp.k_star[1], lrx0, par) - a * M;
    while (R - L > opts.d0) {
      bool right_larger = (R - M > M - L);
      W = right_larger ? (int)ceil(R - nu * (R - M))
                       : (int)ceil(L + nu * (M - L));
      if (W == L || W == R) break;
      HW = H3(W, tmp.k_star[1], lrx0, par) - a * W;
      if (HW > HM) {
        if (right_larger) L = M; else R = M;
        M = W; HM = HW;
      } else {
        if (right_larger) R = W; else L = W;
      }
    }
    for (int k1 = L; k1 <= R; k1++) {
      h_tmp = H3(k1, tmp.k_star[1], lrx0, par);
      if (h_tmp > tmp.h_star) { tmp.k_star[0] = k1; tmp.h_star = h_tmp; }
    }
  } else {
    // Search over k2 with k1 fixed
    L = tmp.k_star[0] + 1;
    R = lrx0.r;
    M = (tmp.k_star[1] == -1) ? (int)floor((L + nu * R) / (1 + nu))
                               : tmp.k_star[1];
    if (opts.rotate && L < R) {
      a = (H3(tmp.k_star[0], R, lrx0, par) - H3(tmp.k_star[0], L, lrx0, par))
          / (R - L);
    }
    HM = H3(tmp.k_star[0], M, lrx0, par) - a * M;
    while (R - L > opts.d0) {
      bool right_larger = (R - M > M - L);
      W = right_larger ? (int)ceil(R - nu * (R - M))
                       : (int)ceil(L + nu * (M - L));
      if (W == L || W == R) break;
      HW = H3(tmp.k_star[0], W, lrx0, par) - a * W;
      if (HW > HM) {
        if (right_larger) L = M; else R = M;
        M = W; HM = HW;
      } else {
        if (right_larger) R = W; else L = W;
      }
    }
    for (int k2 = L; k2 <= R; k2++) {
      h_tmp = H3(tmp.k_star[0], k2, lrx0, par);
      if (h_tmp > tmp.h_star) { tmp.k_star[1] = k2; tmp.h_star = h_tmp; }
    }
  }
}

void OSH3_k0_diag_cpp(kih& tmp,
                      const lrx0& lrx0, const par& par, const opts& opts) {
  int L, R, M, W;
  const double nu = 0.5;
  double HM, HW, h_tmp, a = 0.0;
  tmp.h_star = R_NegInf;
  // Search along the diagonal k -> H3(k, k+1)
  L = lrx0.l + 1;
  R = lrx0.r - 1;
  M = tmp.k_star[0];
  if (opts.rotate && L < R) {
    a = (H3(lrx0.r - 1, lrx0.r,     lrx0, par)
       - H3(lrx0.l + 1, lrx0.l + 2, lrx0, par)) / (R - L);
  }
  HM = H3(M, M + 1, lrx0, par) - a * M;
  while (R - L > opts.d0) {
    bool right_larger = (R - M > M - L);
    W = right_larger ? (int)ceil(R - nu * (R - M))
                     : (int)ceil(L + nu * (M - L));
    if (W == L || W == R) break;
    HW = H3(W, W + 1, lrx0, par) - a * W;
    if (HW > HM) {
      if (right_larger) L = M; else R = M;
      M = W; HM = HW;
    } else {
      if (right_larger) R = W; else L = W;
    }
  }
  for (int k1 = L; k1 <= R; k1++) {
    h_tmp = H3(k1, k1 + 1, lrx0, par);
    if (h_tmp > tmp.h_star) {
      tmp.k_star = {k1, k1 + 1};
      tmp.h_star = h_tmp;
    }
  }
}

void rep_each_cpp(arma::ivec& xx, const arma::ivec& zz, const arma::ivec& SS,
                  int UU, int n) {
  int jj = 0, du;
  for (int i = 0; i < (UU - 1); i++) {
    du = SS[i + 1] - SS[i];
    xx.subvec(jj, jj + du - 1).fill(zz[i]);
    jj += du;
  }
  // Last segment
  du = n - SS[UU - 1];
  xx.subvec(jj, jj + du - 1).fill(zz[UU - 1]);
}
