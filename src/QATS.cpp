#include "QATS.h"

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
List QATS_timer_cpp(int d0, int n_seeds, bool rotate, double n_rep,
                    int n, int m,
                    const arma::vec& logPi, const arma::mat& qq,
                    const arma::mat& GG, arma::ivec SS, int UU) {
  // Preallocation
  arma::ivec xx(n), zz(n);
  SS.resize(n);
  const par par = {m, logPi, qq, GG};
  const opts opts = {d0, n_seeds, rotate};
  // Get starting time point
  auto start = std::chrono::high_resolution_clock::now();
  // Do stuff
  for (int i = 0; i < n_rep; i++) {
    QATS_cpp(xx, zz, SS, UU, n, par, opts);
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

void QATS_cpp(arma::ivec& xx, arma::ivec& zz, arma::ivec& SS, int& UU,
              const int& n, const par& par, const opts& opts) {
  // Declaration
  int u = 0, du, k_star_n = 1;
  kih tmp_1, tmp_2, tmp_3, tmp_star;
  lrx0 lrx0;
  // While-loop
  while (u < UU) {
    // Set window parameters
    lrx0.l = SS[u];
    if ((u + 1) < UU) {
      lrx0.r = SS[u + 1] - 1;
    } else {
      lrx0.r = n - 1;
    }
    du = lrx0.r - lrx0.l + 1;
    if (u > 0) {
      lrx0.x0 = zz[u - 1];
    }
    // Look for a new change point
    if (du == 1) {
      tmp_star.h_star = R_NegInf;
    } else if (du == 2) {
      tmp_star.k_star[0] = lrx0.r;
      argH2_ref(tmp_star, lrx0, par);
      k_star_n = 1;
    } else { // du >= 3
      OSH2_cpp(tmp_2, lrx0, par, opts);
      OSH3_cpp(tmp_3, lrx0, par, opts);
      if (tmp_2.h_star > tmp_3.h_star) {
        tmp_star = tmp_2;
        k_star_n = 1;
      } else {
        tmp_star = tmp_3;
        k_star_n = 2;
      }
    }
    // Best constant path for comparison
    argH1_ref(tmp_1, lrx0, par);
    // Update
    if (tmp_star.h_star > tmp_1.h_star + 1e-7) {
      // Update SS and zz
      if (u < UU-1) {
        for (int j = 1; j < (UU - u); j++) {
          SS[UU - j + k_star_n] = SS[UU - j];
          zz[UU - j + k_star_n] = zz[UU - j];
        }
      }
      SS.subvec(u + 1, u + k_star_n) = tmp_star.k_star.subvec(0, k_star_n - 1);
      zz.subvec(u,     u + k_star_n) = tmp_star.i_star.subvec(0, k_star_n);
      // Update UU
      UU += k_star_n;
    } else {
      zz[u] = tmp_1.i_star[0];
      u++;
    }
  }
  // Modify xx
  rep_each_cpp(xx, zz, SS, UU, n);
}

void OSH2_cpp(kih& res, const lrx0& lrx0, const par& par, const opts& opts) {
  // Declaration
  int L, R, M, W;
  double nu = 1.0/2, HM, HW, a = 0.0;
  bool tmp_bool;
  kih tmp;
  res.h_star = R_NegInf;
  // Initialization
  L = lrx0.l + 1;
  R = lrx0.r;
  M = floor((L + nu*R)/(1+nu));
  if (opts.rotate == true && lrx0.l + 1 < lrx0.r) {
    a = (H2(lrx0.r, lrx0, par) -
      H2(lrx0.l + 1, lrx0, par)) / (lrx0.r - lrx0.l - 1);
  }
  HM = H2(M, lrx0, par) - a * M;
  // Find L:R
  while (R - L > opts.d0) {
    tmp_bool = (R - M > M - L);
    if (tmp_bool) {
      W = ceil(R - nu*(R - M));
    } else {
      W = ceil(L + nu*(M - L));
    }
    if (W == L || W == R) {
      break;
    }
    HW = H2(W, lrx0, par) - a * W;
    if (HW > HM) {
      if (tmp_bool) {
        L = M;
      } else {
        R = M;
      }
      M = W;
      HM = HW;
    } else {
      if (tmp_bool) {
        R = W;
      } else {
        L = W;
      }
    }
  }
  // Full search on L:R
  for (int k1 = L; k1 <= R; k1++) {
    tmp = argH2(k1, lrx0, par);
    if (tmp.h_star > res.h_star) {
      res = tmp;
    }
  }
}

void OSH3_cpp(kih& res, const lrx0& lrx0, const par& par, const opts& opts) {
  // Declaration
  kih tmp;
  res.h_star = R_NegInf;
  // Estimation
  if (lrx0.r - lrx0.l > opts.d0) {
    for (int i = 1; i <= opts.n_seeds; i++) {
      tmp.k_star[0] = -1;
      tmp.k_star[1] = floor(lrx0.l + 2 +
        (i * 1.0)/(opts.n_seeds + 1) * (lrx0.r - lrx0.l - 1.0));
      OSH3_k0_cpp(tmp, lrx0, par, opts);
      if (tmp.h_star > res.h_star) {
        res = tmp;
      }
    }
  } else {
    for (int k1 = (lrx0.l + 1); k1 < lrx0.r; k1++) {
      for (int k2 = (k1 + 1); k2 <= lrx0.r; k2++) {
        tmp = argH3(k1, k2, lrx0, par);
        if (tmp.h_star > res.h_star) {
          res = tmp;
        }
      }
    }
  }
}

void OSH3_k0_cpp(kih& tmp, const lrx0& lrx0, const par& par, const opts& opts) {
  // Bookkeeping
  double h_new = R_NegInf, h_old = R_NegInf;
  // Initialize
  int tau = 0;
  int v = 1, v0 = 20;
  // Start alternating horizontal and vertical searches
  while ( ((h_old < h_new) && (v < v0)) || (v == 1) ) {
    // Overwrite h_old
    h_old = h_new;
    // Horizontal/Vertical search
    OSH3_k0_dir_cpp(tmp, tau, lrx0, par, opts);
    // If on the diagonal, apply OS to k \mapsto H^3(k, k+1) with k^*_1 as the
    // initial probe
    if (tmp.k_star[0] + 1 == tmp.k_star[1]) {
      OSH3_k0_diag_cpp(tmp, lrx0, par, opts);
    }
    // Update h, tau and v
    h_new = tmp.h_star;
    tau = 1 - tau;
    v++;
  }
  // Compute final fit
  argH3_ref(tmp, lrx0, par);
}

void OSH3_k0_dir_cpp(kih& tmp, int tau,
                     const lrx0& lrx0, const par& par, const opts& opts) {
  // Declaration
  int L, R, M, W;
  double nu = 1.0/2, HM, HW, h_tmp, a = 0.0;
  bool tmp_bool;
  tmp.h_star = R_NegInf;
  // Initialization
  if (tau == 0) {
    L = lrx0.l + 1;
    R = tmp.k_star[1] - 1;
    if (tmp.k_star[0] == -1) {
      M = floor((L + nu*R)/(1+nu));
    } else {
      M = tmp.k_star[0];
    }
    if (opts.rotate == true && L < R) {
      a = (H3(R, tmp.k_star[1], lrx0, par) -
        H3(L, tmp.k_star[1], lrx0, par)) / (R - L);
    }
    HM = H3(M, tmp.k_star[1], lrx0, par) - a * M;
  } else {
    L = tmp.k_star[0] + 1;
    R = lrx0.r;
    if (tmp.k_star[1] == -1) {
      M = floor((L + nu*R)/(1+nu));
    } else {
      M = tmp.k_star[1];
    }
    if (opts.rotate == true && L < R) {
      a = (H3(tmp.k_star[0], R, lrx0, par) -
        H3(tmp.k_star[0], L, lrx0, par)) / (R - L);
    }
    HM = H3(tmp.k_star[0], M, lrx0, par) - a * M;
  }
  // Find L:R
  while (R - L > opts.d0) {
    tmp_bool = (R - M > M - L);
    if (tmp_bool) {
      W = ceil(R - nu*(R - M));
    } else {
      W = ceil(L + nu*(M - L));
    }
    if (W == L || W == R) {
      break;
    }
    if (tau == 0) {
      HW = H3(W, tmp.k_star[1], lrx0, par) - a * W;
    } else {
      HW = H3(tmp.k_star[0], W, lrx0, par) - a * W;
    }
    if (HW > HM) {
      if (tmp_bool) {
        L = M;
      } else {
        R = M;
      }
      M = W;
      HM = HW;
    } else {
      if (tmp_bool) {
        R = W;
      } else {
        L = W;
      }
    }
  }
  // Full search on L:R
  if (tau == 0) {
    for (int k1 = L; k1 <= R; k1++) {
      h_tmp = H3(k1, tmp.k_star[1], lrx0, par);
      if (h_tmp > tmp.h_star) {
        tmp.k_star[0] = k1;
        tmp.h_star = h_tmp;
      }
    }
  } else {
    for (int k2 = L; k2 <= R; k2++) {
      h_tmp = H3(tmp.k_star[0], k2, lrx0, par);
      if (h_tmp > tmp.h_star) {
        tmp.k_star[1] = k2;
        tmp.h_star = h_tmp;
      }
    }
  }
}

void OSH3_k0_diag_cpp(kih& tmp,
                      const lrx0& lrx0, const par& par, const opts& opts) {
  // Declaration
  int L, R, M, W;
  double nu = 1.0/2, HM, HW, h_tmp, a = 0.0;
  bool tmp_bool;
  tmp.h_star = R_NegInf;
  // Initialization
  L = lrx0.l + 1;
  R = lrx0.r - 1;
  M = tmp.k_star[0];
  if (opts.rotate == true && L < R) {
    a = (H3(lrx0.r - 1, lrx0.r,     lrx0, par) -
      H3(lrx0.l + 1, lrx0.l + 2, lrx0, par)) / (R - L);
  }
  HM = H3(M, M + 1, lrx0, par) - a * M;
  // Find L:R
  while (R - L > opts.d0) {
    tmp_bool = (R - M > M - L);
    if (tmp_bool) {
      W = ceil(R - nu*(R - M));
    } else {
      W = ceil(L + nu*(M - L));
    }
    if (W == L || W == R) {
      break;
    }
    HW = H3(W, W + 1, lrx0, par) - a * W;
    if (HW > HM) {
      if (tmp_bool) {
        L = M;
      } else {
        R = M;
      }
      M = W;
      HM = HW;
    } else {
      if (tmp_bool) {
        R = W;
      } else {
        L = W;
      }
    }
  }
  // Full search on L:R
  for (int k1 = L; k1 <= R; k1++) {
    h_tmp = H3(k1, k1 + 1, lrx0, par);
    if (h_tmp > tmp.h_star) {
      tmp.k_star = {k1, k1 + 1};
      tmp.h_star = h_tmp;
    }
  }
}

void rep_each_cpp(arma::ivec& xx, const arma::ivec& zz, const arma::ivec& SS,
                  int& UU, int n) {
  int jj = 0, du = 0;
  for (int i = 0; i < (UU - 1); i++) {
    du = SS[i + 1] - SS[i];
    xx.subvec(jj, jj + du - 1).fill(zz[i]);
    jj += du;
  }
  // Setting last part of xx, i.e. for u = UU - 1
  du = n - SS[UU - 1];
  xx.subvec(jj, jj + du - 1).fill(zz[UU - 1]);
}

/*
# add 'star star R' to make it run!
# Some debeuging tricks:
# Add '// Rcout << "Try... ";' before and '// Rcout << "Worked!\n ";' after any
# function call until the problem is found!

 library(rHMM)
 opts  <- list(n.seeds = 5, d0 = 3,
 rotate = TRUE, n.rep = 1)
 SS <- 0; UU <- 1
 n <- 1e3 + 1; m <- 5; sigma <- 1; K <- 1

 for (i in 1:1e4) {
 set.seed(i)
 par <- sample.norm.HMM(n, m, sigma, K)

 res.QATS <- QATS_timer_cpp(opts$d0, opts$n.seeds,
 opts$rotate, opts$n.rep,
 n, m, par$logPi, par$qq, par$GG,
 SS, UU)
 }
 */
