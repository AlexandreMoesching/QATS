# ── Internal helpers ──────────────────────────────────────────────────────────

#' Validate and fill in default QATS options
#'
#' @param opts A list of options (partially or fully specified).
#' @param par  The HMM parameter list (used only to check SS bounds).
#' @return A fully validated and completed opts list.
#' @keywords internal
.validate_opts <- function(opts, par) {
  stopifnot(is.list(opts))

  if (is.null(opts$n.seeds)) {
    opts$n.seeds <- 3L
  } else {
    stopifnot(is.numeric(opts$n.seeds), length(opts$n.seeds) == 1L,
              opts$n.seeds >= 1L)
  }

  if (is.null(opts$d0)) {
    opts$d0 <- 3L
  } else {
    stopifnot(is.numeric(opts$d0), length(opts$d0) == 1L, opts$d0 >= 1L)
  }

  if (is.null(opts$rotate)) {
    opts$rotate <- FALSE
  } else {
    stopifnot(is.logical(opts$rotate), length(opts$rotate) == 1L)
  }

  opts
}

#' Validate and normalise the SS / UU initial-partition fields
#'
#' @param opts   Opts list (already validated by .validate_opts).
#' @param par    HMM parameter list.
#' @param cpp    Logical; if TRUE, return C++-style 0-based left endpoints.
#' @return A list with elements \code{SS} and \code{UU}.
#' @keywords internal
.parse_SS <- function(opts, par, cpp = FALSE) {
  if (is.null(opts$SS)) {
    if (cpp) {
      return(list(SS = 0L, UU = 1L))
    } else {
      return(list(SS = matrix(c(1L, par$n), nrow = 1L), UU = 1L))
    }
  }
  stopifnot(
    is.matrix(opts$SS), ncol(opts$SS) == 2L,
    all(opts$SS[, 1L] >= 1L & opts$SS[, 2L] <= par$n),
    all(opts$SS[, 1L] <= opts$SS[, 2L])
  )
  if (cpp) {
    list(SS = opts$SS[, 1L] - 1L, UU = nrow(opts$SS))
  } else {
    list(SS = opts$SS, UU = nrow(opts$SS))
  }
}

#' Validate the common par object
#' @keywords internal
.validate_par <- function(par) {
  stopifnot(
    is.list(par),
    is.numeric(par$n),     length(par$n)     == 1L, par$n >= 1L,
    is.numeric(par$m),     length(par$m)     == 1L, par$m >= 1L,
    is.numeric(par$logPi), length(par$logPi) == par$m,
    is.matrix(par$qq),     all(dim(par$qq)   == c(par$m, par$m)),
    is.matrix(par$GG),     all(dim(par$GG)   == c(par$m, par$n))
  )
}

# ── Main exported functions ───────────────────────────────────────────────────

#' Quick Adaptive Ternary Segmentation (pure R)
#'
#' Given an HMM parameter object, do a greedy ternary segmentation of the
#' observation sequence by repeatedly inserting up to two change-points per
#' segment where they most improve the log-likelihood.
#'
#' @param par  A list with components
#'   \describe{
#'     \item{n}{(integer) length of the sequence}
#'     \item{m}{(integer) number of hidden states}
#'     \item{logPi}{(numeric vector length m) initial log-probabilities}
#'     \item{qq}{(matrix m×m) log-transition matrix}
#'     \item{GG}{(matrix m×n) cumulative log-emissions}
#'   }
#' @param opts  A list of options
#'   \describe{
#'     \item{n.seeds}{(integer ≥ 1) number of seeds for optimistic search (default 3)}
#'     \item{d0}{(integer ≥ 1) minimum window width (default 3)}
#'     \item{SS}{(integer matrix, n_segments × 2) optional initial partition; each row \[l, r\]}
#'     \item{rotate}{(logical) apply linear rotation adjustment (default FALSE)}
#'   }
#' @return A list with
#'   \describe{
#'     \item{xx}{(integer vector length n) estimated hidden path}
#'     \item{logp}{(numeric) log-probability of \code{xx}}
#'     \item{time}{\code{difftime} object giving CPU time}
#'   }
#' @examples
#' \dontrun{
#' set.seed(1)
#' par <- sample.HMM(
#'   n = 200, m = 3, K = 5,
#'   emi.dist = "normal",
#'   emi.param = list(mu = 1:3, sigma = rep(1, 3))
#' )
#' res <- QATS.R(par)
#' str(res)
#' }
#' @export
QATS.R <- function(par,
                   opts = list()) {
  .validate_par(par)
  opts  <- .validate_opts(opts, par)
  part  <- .parse_SS(opts, par, cpp = FALSE)
  SS    <- part$SS
  UU    <- part$UU

  x0 <- 1L
  u  <- 1L
  time <- Sys.time()

  tmp <- argH1(1L, par$n, x0, par)
  zz  <- rep(tmp$i_star, UU)

  while (u <= UU) {
    lu <- SS[u, 1L]; ru <- SS[u, 2L]; du <- ru - lu + 1L
    x0 <- if (u > 1L) zz[u - 1L] else 1L

    if (du == 1L) {
      k_star <- NULL; i_star <- 1L; h_star <- -Inf
    } else if (du == 2L) {
      tmp1   <- OSH2(lu, ru, x0, par, opts)
      k_star <- tmp1$k_star; i_star <- tmp1$i_star; h_star <- tmp1$h_star
    } else {
      tmp1 <- OSH2(lu, ru, x0, par, opts)
      tmp2 <- OSH3(lu, ru, x0, par, opts)
      if (tmp1$h_star >= tmp2$h_star) {
        k_star <- tmp1$k_star; i_star <- tmp1$i_star; h_star <- tmp1$h_star
      } else {
        k_star <- tmp2$k_star; i_star <- tmp2$i_star; h_star <- tmp2$h_star
      }
    }

    tmp          <- argH1(lu, ru, x0, par)
    h_star.const <- tmp$h_star
    i_star.const <- tmp$i_star

    if (h_star > h_star.const + 1e-7) {
      SS <- rbind(
        if (u > 1L) SS[seq_len(u - 1L), , drop = FALSE] else NULL,
        cbind(c(lu, k_star), c(k_star - 1L, ru)),
        if (u < UU) SS[seq(u + 1L, UU), , drop = FALSE] else NULL
      )
      zz  <- c(zz[seq_len(u - 1L)], i_star, zz[u + seq_len(UU - u)])
      UU  <- length(zz)
    } else {
      zz[u] <- i_star.const
      u      <- u + 1L
    }
  }

  time <- difftime(Sys.time(), time, units = "secs")
  xx   <- xx_SS.zz(SS, zz)
  logp <- G0(xx, par)
  list(xx = xx, logp = logp, time = time)
}

#' Quick Adaptive Ternary Segmentation (C++ backend)
#'
#' Exactly like \code{QATS.R()}, but delegates the inner loop to the C++
#' implementation for maximum performance, and optionally times repeated runs.
#'
#' @inheritParams QATS.R
#' @param opts  Same as \code{QATS.R()} with one extra field:
#'   \describe{
#'     \item{n.rep}{(integer ≥ 1) repetitions for timing (default 1)}
#'   }
#'
#' @return A list with
#'   \describe{
#'     \item{xx}{(integer vector length n) estimated hidden path}
#'     \item{logp}{(numeric) log-probability of \code{xx}}
#'     \item{time}{(numeric) average seconds per call}
#'     \item{time_ms}{(numeric) average microseconds per call}
#'   }
#' @export
QATS.CPP <- function(par,
                     opts = list()) {
  .validate_par(par)
  opts <- .validate_opts(opts, par)

  if (is.null(opts$n.rep)) {
    opts$n.rep <- 1L
  } else {
    stopifnot(is.numeric(opts$n.rep), length(opts$n.rep) == 1L,
              opts$n.rep >= 1L)
  }

  part <- .parse_SS(opts, par, cpp = TRUE)

  res  <- QATS_timer_cpp(
    opts$d0, opts$n.seeds, opts$rotate, opts$n.rep,
    par$n, par$m, par$logPi, par$qq, par$GG,
    part$SS, part$UU
  )
  logp <- G0(c(res$xx), par)
  list(xx = c(res$xx), logp = logp, time = res$time, time_ms = res$time_ms)
}

#' Quick Adaptive Ternary Segmentation – Step-by-Step Display
#'
#' An interactive variant of QATS that pauses after each partition update and
#' plots the current segmentation, the Viterbi path, and the true hidden
#' sequence.
#'
#' @param xx.0   Integer vector of length \code{n}: true hidden-state sequence.
#' @param xx.Vit Integer vector of length \code{n}: Viterbi-estimated sequence.
#' @inheritParams QATS.R
#'
#' @return A list with
#'   \describe{
#'     \item{xx}{(integer vector length n) estimated hidden path}
#'     \item{logp}{(numeric) log-probability of \code{xx}}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' par    <- sample.HMM(n = 100, m = 3, K = 5,
#'                      emi.dist = "normal",
#'                      emi.param = list(mu = 1:3, sigma = rep(1, 3)))
#' resVit <- Viterbi.CPP(par)
#' QATS.display(par$xx, resVit$xx, par)
#' }
#' @export
QATS.display <- function(xx.0, xx.Vit, par,
                         opts = list()) {
  stopifnot(
    is.numeric(xx.0), is.numeric(xx.Vit),
    length(xx.0) == length(xx.Vit)
  )
  .validate_par(par)
  stopifnot(par$n == length(xx.0))
  opts <- .validate_opts(opts, par)
  part <- .parse_SS(opts, par, cpp = FALSE)
  SS   <- part$SS
  UU   <- part$UU

  x0 <- 1L
  u  <- 1L
  n  <- par$n
  CP <- (2L:n)[xx.0[1L:(n - 1L)] != xx.0[2L:n]]

  tmp <- argH1(1L, par$n, x0, par)
  zz  <- rep(tmp$i_star, UU)

  xx <- xx_SS.zz(SS, zz)
  display.0(xx.0, xx.Vit, xx, par, SS)
  cat("Next search interval: ", u, "\nPress [enter] to continue\n", sep = "")
  readline()

  while (u <= UU) {
    lu <- SS[u, 1L]; ru <- SS[u, 2L]; Su <- lu:ru; du <- length(Su)
    x0 <- if (u > 1L) zz[u - 1L] else 1L

    if (du > 1L) {
      tmp.plot1 <- H2_vec(Su, du, x0, par)
      tmp.plot2 <- if (du > 2L) H3_mat(Su, du, x0, par)
                   else         list(res = NULL, t_star = NULL)
      display.1(xx.0, xx.Vit, xx, par, SS, NULL, CP, lu,
                tmp.plot1$res, tmp.plot1$t_star, NULL,
                tmp.plot2$res, tmp.plot2$t_star, NULL)
      cat("Function to maximize with exact maximum (black)\n",
          "Press [enter] to continue\n", sep = "")
      readline()
    } else {
      cat("Length of partition = 1, nothing to look for...\n",
          "Press [enter] to continue\n", sep = "")
      readline()
    }

    if (du == 1L) {
      k_star <- NULL; i_star <- 1L; h_star <- -Inf
    } else if (du == 2L) {
      tmp1   <- OSH2(lu, ru, x0, par, opts)
      k_star <- tmp1$k_star; i_star <- tmp1$i_star; h_star <- tmp1$h_star
    } else {
      tmp1 <- OSH2(lu, ru, x0, par, opts)
      tmp2 <- OSH3(lu, ru, x0, par, opts)
      if (tmp1$h_star >= tmp2$h_star) {
        k_star <- tmp1$k_star; i_star <- tmp1$i_star; h_star <- tmp1$h_star
      } else {
        k_star <- tmp2$k_star; i_star <- tmp2$i_star; h_star <- tmp2$h_star
      }
    }

    if (du > 1L) {
      if (length(k_star) == 0L) {
        cat("No maximum found\nPress [enter] to continue\n")
      } else {
        k2_arg <- if (du == 2L) NULL else t(tmp2$k_star)
        display.1(xx.0, xx.Vit, xx, par, SS, NULL, CP, lu,
                  tmp.plot1$res, tmp.plot1$t_star, tmp1$k_star,
                  tmp.plot2$res, tmp.plot2$t_star, k2_arg)
        cat(if (length(k_star) == 1L) "One change point found.\n"
            else                       "Two change points found.\n")
      }
      readline()
    }

    tmp          <- argH1(lu, ru, x0, par)
    h_star.const <- tmp$h_star
    i_star.const <- tmp$i_star

    if (h_star > h_star.const + 1e-7) {
      SS <- rbind(
        if (u > 1L) SS[seq_len(u - 1L), , drop = FALSE] else NULL,
        cbind(c(lu, k_star), c(k_star - 1L, ru)),
        if (u < UU) SS[seq(u + 1L, UU), , drop = FALSE] else NULL
      )
      zz  <- c(zz[seq_len(u - 1L)], i_star, zz[u + seq_len(UU - u)])
      UU  <- length(zz)
    } else {
      zz[u] <- i_star.const
      u      <- u + 1L
      k_star <- NULL
    }

    xx <- xx_SS.zz(SS, zz)
    display.0(xx.0, xx.Vit, xx, par, SS, k_star)
    if (!is.null(k_star)) cat("New change point(s) added (green)\n")
    else                   cat("No new change point\n")
    if (u <= UU) {
      cat("Next search interval: ", u, "\nPress [enter] to continue\n", sep = "")
      readline()
    }
  }

  xx   <- xx_SS.zz(SS, zz)
  logp <- G0(xx, par)
  cat("Final fit\n")
  graphics::par(mfrow = c(1, 1), mar = c(4.2, 4.2, 0.2, 0.2))
  display.result(xx.0, xx.Vit, xx, par)
  list(xx = xx, logp = logp)
}

# ── Optimistic search primitives (internal) ───────────────────────────────────

#' Basic one-dimensional optimistic (ternary) search
#'
#' @param L      Integer: left boundary (inclusive).
#' @param R      Integer: right boundary (inclusive), must satisfy R > L.
#' @param d0     Integer ≥ 1: minimum search window width before full scan.
#' @param fun    Function of one integer k returning a numeric score.
#' @param argfun Function of one integer k returning a list with at least
#'   \code{h_star} (numeric score) and \code{i_star} (associated states).
#' @param M      Integer: initial probe position; computed automatically if < 0.
#'
#' @return A list with \code{k_star}, \code{i_star}, and \code{h_star}.
#' @keywords internal
OS <- function(L, R, d0, fun, argfun, M = -1L) {
  nu <- 0.5
  if (M < 0L) M <- floor((L + nu * R) / (1 + nu))
  HM <- fun(M)

  while (R - L > d0) {
    right_larger <- (R - M > M - L)
    W <- if (right_larger) ceiling(R - nu * (R - M)) else ceiling(L + nu * (M - L))
    if (W == L || W == R) break
    HW <- fun(W)
    if (HW > HM) {
      if (right_larger) L <- M else R <- M
      M <- W; HM <- HW
    } else {
      if (right_larger) R <- W else L <- W
    }
  }

  k_star <- integer(0L); i_star <- integer(0L); h_star <- -Inf
  for (k in seq.int(L, R)) {
    tmp <- argfun(k)
    if (tmp$h_star > h_star) {
      k_star <- k; i_star <- tmp$i_star; h_star <- tmp$h_star
    }
  }
  list(k_star = k_star, i_star = i_star, h_star = h_star)
}

#' 1-dimensional optimistic search for H2
#' @keywords internal
OSH2 <- function(l, r, x0, par, opts) {
  if (opts$rotate && l + 1L < r) {
    H2_L <- H2(l, l + 1L, r, x0, par)
    H2_R <- H2(l, r,      r, x0, par)
    a    <- (H2_R - H2_L) / (r - l - 1L)
    OS(l + 1L, r, opts$d0,
       function(k) H2(l, k, r, x0, par) - a * k,
       function(k) argH2(l, k, r, x0, par))
  } else {
    OS(l + 1L, r, opts$d0,
       function(k) H2(l, k, r, x0, par),
       function(k) argH2(l, k, r, x0, par))
  }
}

#' 2-dimensional optimistic search for H3
#' @keywords internal
OSH3 <- function(l, r, x0, par, opts) {
  if (r - l > opts$d0) {
    kk <- cbind(
      -1L,
      unique(sort(floor(
        l + 2L + seq_len(opts$n.seeds) / (opts$n.seeds + 1L) * (r - l - 1L)
      )))
    )
    hh <- rep(-Inf, nrow(kk))
    for (d in seq_len(nrow(kk))) {
      tmp      <- OSH3_k0(kk[d, ], l, r, x0, par, opts)
      kk[d, ]  <- tmp$k_star
      hh[d]    <- tmp$h_star
    }
    k_star <- kk[which.max(hh), ]
  } else {
    h_star <- -Inf; k_star <- rep(NA_integer_, 2L)
    for (k1 in seq.int(l + 1L, r - 1L)) {
      for (k2 in seq.int(k1 + 1L, r)) {
        tmp <- H3(l, k1, k2, r, x0, par)
        if (tmp > h_star) { h_star <- tmp; k_star <- c(k1, k2) }
      }
    }
  }
  res <- argH3(l, k_star[1L], k_star[2L], r, x0, par)
  list(k_star = k_star, i_star = res$i_star, h_star = res$h_star)
}

#' 2-dimensional optimistic search for H3 with a starting point
#' @keywords internal
OSH3_k0 <- function(k0, l, r, x0, par, opts) {
  h_new <- h_old <- -Inf
  tau <- 0L; v <- 1L; v0 <- 20L

  while (((h_old < h_new) & (v < v0)) || v == 1L) {
    h_old <- h_new

    if (tau == 0L) {
      if (opts$rotate && l + 1L < k0[2L] - 1L) {
        a <- (H3(l, k0[2L] - 1L, k0[2L], r, x0, par) -
              H3(l, l + 1L,      k0[2L], r, x0, par)) / (k0[2L] - l - 2L)
        tmp <- OS(l + 1L, k0[2L] - 1L, opts$d0,
                  function(k) H3(l, k, k0[2L], r, x0, par) - a * k,
                  function(k) argH3(l, k, k0[2L], r, x0, par),
                  k0[1L])
      } else {
        tmp <- OS(l + 1L, k0[2L] - 1L, opts$d0,
                  function(k) H3(l, k, k0[2L], r, x0, par),
                  function(k) argH3(l, k, k0[2L], r, x0, par),
                  k0[1L])
      }
    } else {
      if (opts$rotate && k0[1L] + 1L < r) {
        a <- (H3(l, k0[1L], r,         r, x0, par) -
              H3(l, k0[1L], k0[1L]+1L, r, x0, par)) / (r - k0[1L] - 1L)
        tmp <- OS(k0[1L] + 1L, r, opts$d0,
                  function(k) H3(l, k0[1L], k, r, x0, par) - a * k,
                  function(k) argH3(l, k0[1L], k, r, x0, par),
                  k0[2L])
      } else {
        tmp <- OS(k0[1L] + 1L, r, opts$d0,
                  function(k) H3(l, k0[1L], k, r, x0, par),
                  function(k) argH3(l, k0[1L], k, r, x0, par),
                  k0[2L])
      }
    }
    k0[1L + tau] <- tmp$k_star

    if (k0[1L] + 1L == k0[2L]) {
      if (opts$rotate && l + 2L < r) {
        a <- (H3(l, r - 1L, r,         r, x0, par) -
              H3(l, l + 1L, l + 2L,    r, x0, par)) / (r - l - 2L)
        tmp <- OS(l + 1L, r - 1L, opts$d0,
                  function(k) H3(l, k, k + 1L, r, x0, par) - a * k,
                  function(k) argH3(l, k, k + 1L, r, x0, par),
                  k0[1L])
      } else {
        tmp <- OS(l + 1L, r - 1L, opts$d0,
                  function(k) H3(l, k, k + 1L, r, x0, par),
                  function(k) argH3(l, k, k + 1L, r, x0, par),
                  k0[1L])
      }
      k0 <- c(tmp$k_star, tmp$k_star + 1L)
    }

    h_new <- H3(l, k0[1L], k0[2L], r, x0, par)
    tau   <- 1L - tau
    v     <- v + 1L
  }

  res <- argH3(l, k0[1L], k0[2L], r, x0, par)
  list(k_star = k0, i_star = res$i_star, h_star = res$h_star)
}
