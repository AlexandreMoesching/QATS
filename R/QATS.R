#' Quick Adaptive Ternary Segmentation (pure R)
#'
#' Given an HMM parameter object, do a greedy ternary segmentation of the
#' observation sequence by repeatedly inserting up to two change‐points per
#' segment where they most improve the log-likelihood.
#'
#' @param par  A list with components
#'   \describe{
#'     \item{n}{(integer) length of the sequence}
#'     \item{m}{(integer) number of hidden states}
#'     \item{yy}{(numeric vector) observations}
#'     \item{logPi}{(numeric vector length m) initial log-probabilities}
#'     \item{qq}{(matrix m×m) log-transition matrix}
#'     \item{GG}{(matrix m×n) cumulative log-emissions}
#'   }
#' @param opts  A list of options
#'   \describe{
#'     \item{n.seeds}{(integer ≥ 1) number of seeds for optimistic search}
#'     \item{d0}{(integer ≥ 1) minimum window width}
#'     \item{SS}{(integer matrix 2×k) optional initial partition rows = \[l,r\]}
#'     \item{rotate}{(logical) whether to apply linear “rotation” adjustment}
#'     \item{n.rep}{integer ≥ 1; number of repetitions for timing of C++ version (default 1)}
#'   }
#' @return A list with
#'   \describe{
#'     \item{xx}{(integer vector length n) estimated hidden path}
#'     \item{logp}{(numeric) log-probability of \code{xx}}
#'     \item{time}{\code{difftime} object giving CPU time}
#'   }
#' @examples
#' \dontrun{
#' par <- sample.HMM(
#'   n = 100, m = 3, K = 5,
#'   emi.dist = "normal",
#'   emi.param = list(mu = 1:3, sigma = rep(1, 3))
#' )
#' res <- QATS.R(par)
#' str(res)
#' }
#' @export
QATS.R <- function(par,
                   opts = list(
                     n.seeds = 3L,
                     d0 = 3L,
                     SS = NULL,
                     rotate = FALSE
                   )) {
  # Input validation
  stopifnot(
    is.list(par),
    is.numeric(par$n),     length(par$n)     == 1,       par$n >= 1,
    is.numeric(par$m),     length(par$m)     == 1,       par$m >= 1,
    is.numeric(par$logPi), length(par$logPi) == par$m,
    is.matrix(par$qq),     all(dim(par$qq)   == c(par$m, par$m)),
    is.matrix(par$GG),     all(dim(par$GG)   == c(par$m, par$n)),
    is.list(opts)
  )

  # Set missing options
  if (is.null(opts$n.seeds)) {
    opts$n.seeds <- 3L
  } else {
    stopifnot(is.numeric(opts$n.seeds), length(opts$n.seeds) == 1L, opts$n.seeds >= 1L)
  }

  if (is.null(opts$d0)) {
    opts$d0 <- 3L
  } else {
    stopifnot(is.numeric(opts$d0), length(opts$d0) == 1L, opts$d0 >= 1L)
  }

  if (is.null(opts$SS)) {
    SS <- matrix(c(1L, par$n), nrow = 1L)
    UU <- 1L
  } else {
    stopifnot(
      is.matrix(opts$SS), ncol(opts$SS) == 2L,
      all(opts$SS[,1L] >= 1L & opts$SS[,2L] <= par$n),
      all(opts$SS[,1L] <= opts$SS[,2L])
    )
    SS <- opts$SS
    UU <- nrow(SS)
  }

  if (is.null(opts$rotate)) {
    opts$rotate <- FALSE
  } else {
    stopifnot(is.logical(opts$rotate),  length(opts$rotate) == 1L)
  }

  # Initialize
  x0 <- 1L
  u <- 1L

  # Find best constant vector
  time <- Sys.time()
  tmp <- argH1(1L, par$n, x0, par)
  zz <- rep(tmp$i_star, UU)

  # Loop
  while (u <= UU) {
    # Current part of the partition under investigation
    lu <- SS[u, 1L]
    ru <- SS[u, 2L]
    du <- ru - lu + 1L

    # Last value of the path before the current partition
    x0 <- ifelse(u > 1L, zz[u - 1L], 1L)

    # Look for 2 CPs (if du >= 3), 1 CP (if du = 2), or none (if du = 1)
    if (du == 1L) {
      k_star <- NULL
      i_star <- 1L
      h_star <- -Inf
    } else if (du == 2L) {
      tmp1 <- OSH2(lu, ru, x0, par, opts)
      k_star <- tmp1$k_star
      i_star <- tmp1$i_star
      h_star <- tmp1$h_star
    } else { # du >= 3
      tmp1 <- OSH2(lu, ru, x0, par, opts)
      tmp2 <- OSH3(lu, ru, x0, par, opts)
      if (tmp1$h_star >= tmp2$h_star) {
        k_star <- tmp1$k_star
        i_star <- tmp1$i_star
        h_star <- tmp1$h_star
      } else {
        k_star <- tmp2$k_star
        i_star <- tmp2$i_star
        h_star <- tmp2$h_star
      }
    }

    # Best constant path for comparison
    tmp <- argH1(lu, ru, x0, par)
    h_star.const <- tmp$h_star
    i_star.const <- tmp$i_star

    # Update
    if (h_star > h_star.const + 1e-7) {
      SS <- rbind(
        SS[-(u:UU), ],
        cbind(c(lu, k_star), c(k_star - 1L, ru)),
        SS[-(1L:u), ]
      )
      zz <- c(zz[-(u:UU)], i_star, zz[-(1L:u)])
      UU <- length(zz)
    } else {
      zz[u] <- i_star.const
      u <- u + 1L
      k_star <- NULL
    }
  }

  time <- difftime(Sys.time(), time, units = "secs")
  xx <- xx_SS.zz(SS, zz)
  logp <- G0(xx, par)
  return(list(
    xx = xx,
    logp = logp,
    time = time
  ))
}

#' Quick Adaptive Ternary Segmentation (C++ backend)
#'
#' Exactly like \code{QATS.R()}, but calls the underlying C++ implementation
#' for the inner loop, and optionally times it.
#'
#' @inheritParams QATS.R
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
                     opts = list(
                       n.seeds = 3L,
                       d0 = 3L,
                       SS = NULL,
                       rotate = FALSE,
                       n.rep = 1L
                     )) {
  # Input validation
  stopifnot(
    is.list(par),
    is.numeric(par$n),     length(par$n)     == 1,       par$n >= 1,
    is.numeric(par$m),     length(par$m)     == 1,       par$m >= 1,
    is.numeric(par$logPi), length(par$logPi) == par$m,
    is.matrix(par$qq),     all(dim(par$qq)   == c(par$m, par$m)),
    is.matrix(par$GG),     all(dim(par$GG)   == c(par$m, par$n)),
    is.list(opts)
  )

  # Set missing options
  if (is.null(opts$n.seeds)) {
    opts$n.seeds <- 3L
  } else {
    stopifnot(is.numeric(opts$n.seeds), length(opts$n.seeds) == 1L, opts$n.seeds >= 1L)
  }

  if (is.null(opts$d0)) {
    opts$d0 <- 3L
  } else {
    stopifnot(is.numeric(opts$d0), length(opts$d0) == 1L, opts$d0 >= 1L)
  }

  if (is.null(opts$SS)) {
    SS <- 0L
    UU <- 1L
  } else {
    stopifnot(
      is.matrix(opts$SS), ncol(opts$SS) == 2L,
      all(opts$SS[,1L] >= 1L & opts$SS[,2L] <= par$n),
      all(opts$SS[,1L] <= opts$SS[,2L])
    )
    SS <- opts$SS[, 1L] - 1L
    UU <- length(SS)
  }

  if (is.null(opts$rotate)) {
    opts$rotate <- FALSE
  } else {
    stopifnot(is.logical(opts$rotate),  length(opts$rotate) == 1L)
  }

  if (is.null(opts$n.rep)) {
    opts$n.rep <- 1L
  } else {
    stopifnot(is.numeric(opts$n.rep),  length(opts$n.rep) == 1L)
  }

  # Start estimation
  res <- QATS_timer_cpp(
    opts$d0, opts$n.seeds, opts$rotate, opts$n.rep,
    par$n, par$m, par$logPi, par$qq, par$GG, SS, UU
  )

  # Compute path probability
  logp <- G0(c(res$xx), par)
  return(list(
    xx = c(res$xx),
    logp = logp,
    time = res$time,
    time_ms = res$time_ms
  ))
}

#' Quick Adaptive Ternary Segmentation – Step‐by‐Step Display
#'
#' An interactive variant of QATS that pauses after each partition update and
#' plots the current segmentation, the Viterbi path, and the true hidden sequence.
#'
#' @param xx.0  Integer (or numeric) vector of length \code{n}: the true hidden‐state sequence.
#' @param xx.Vit  Integer (or numeric) vector of length \code{n}: the Viterbi‐estimated sequence.
#' @inheritParams QATS.R
#'
#' @details After initializing to the best constant fit on each segment, \code{QATS.display()}
#' enters a greedy loop.  At each step it:
#' \enumerate{
#'   \item Highlights the sub‐interval currently under search.
#'   \item Plots the marginal gain functions (\code{H2}, \code{H3}) and their maxima.
#'   \item Pauses for user input before inserting up to two new change‐points.
#' }
#'
#' @return A list with
#'   \describe{
#'     \item{xx}{(integer vector length n) estimated hidden path}
#'     \item{logp}{(numeric) log-probability of \code{xx}}
#'   }
#'
#' @examples
#' \dontrun{
#' par <- sample.HMM(
#'   n = 100, m = 3, K = 5,
#'   emi.dist = "normal",
#'   emi.param = list(mu = 1:3, sigma = rep(1, 3))
#' )
#' xx0 <- par$xx
#' resVit <- Viterbi.CPP(par)
#' QATS.display(xx0, resVit$xx, par)
#' }
#'
#' @export
QATS.display <- function(xx.0, xx.Vit, par,
                         opts = list(
                           n.seeds = 3L,
                           d0 = 3L,
                           SS = NULL,
                           rotate = FALSE
                         )) {
  # Input validation
  stopifnot(
    is.numeric(xx.0), is.numeric(xx.Vit), length(xx.0) == length(xx.Vit),
    is.list(par),
    is.numeric(par$n),     length(par$n)     == 1,       par$n >= 1,
    par$n == length(xx.0),
    is.numeric(par$m),     length(par$m)     == 1,       par$m >= 1,
    is.numeric(par$logPi), length(par$logPi) == par$m,
    is.matrix(par$qq),     all(dim(par$qq)   == c(par$m, par$m)),
    is.matrix(par$GG),     all(dim(par$GG)   == c(par$m, par$n)),
    is.list(opts)
  )

  # Set missing options
  if (is.null(opts$n.seeds)) {
    opts$n.seeds <- 3L
  } else {
    stopifnot(is.numeric(opts$n.seeds), length(opts$n.seeds) == 1L, opts$n.seeds >= 1L)
  }

  if (is.null(opts$d0)) {
    opts$d0 <- 3L
  } else {
    stopifnot(is.numeric(opts$d0), length(opts$d0) == 1L, opts$d0 >= 1L)
  }

  if (is.null(opts$SS)) {
    SS <- matrix(c(1L, par$n), nrow = 1L)
    UU <- 1L
  } else {
    stopifnot(
      is.matrix(opts$SS), ncol(opts$SS) == 2L,
      all(opts$SS[,1L] >= 1L & opts$SS[,2L] <= par$n),
      all(opts$SS[,1L] <= opts$SS[,2L])
    )
    SS <- opts$SS
    UU <- nrow(SS)
  }

  if (is.null(opts$rotate)) {
    opts$rotate <- FALSE
  } else {
    stopifnot(is.logical(opts$rotate),  length(opts$rotate) == 1L)
  }

  # Initialize
  x0 <- 1L
  u <- 1L

  # Number of change points
  n <- par$n
  CP <- (2L:n)[xx.0[1L:(n - 1L)] != xx.0[2L:n]]

  # Find best constant vector
  tmp <- argH1(1L, par$n, x0, par)
  zz <- rep(tmp$i_star, UU)

  # For plotting purposes
  xx <- xx_SS.zz(SS, zz)
  display.0(xx.0, xx.Vit, xx, par, SS)
  cat("Next search interval: ", u, "\n",
      "Press [enter] to continue\n",
      sep = ""
  )
  line <- readline()

  # Loop
  while (u <= UU) {
    # Current part of the partition under investigation
    lu <- SS[u, 1L]
    ru <- SS[u, 2L]
    Su <- lu:ru
    du <- length(Su)

    # Last value of the path before the current partition
    x0 <- ifelse(u > 1L, zz[u - 1L], 1L)

    # Plot the functions to maximize
    if (du > 1L) {
      tmp.plot1 <- H2_vec(Su, du, x0, par)
      if (du > 2L) {
        tmp.plot2 <- H3_mat(Su, du, x0, par)
      } else {
        tmp.plot2 <- list(res = NULL, t_star = NULL)
      }
      display.1(
        xx.0, xx.Vit, xx, par, SS, NULL, CP, lu,
        tmp.plot1$res, tmp.plot1$t_star, NULL,
        tmp.plot2$res, tmp.plot2$t_star, NULL
      )
      cat("Function to maximize with exact maximum (black)\n",
          "Press [enter] to continue\n",
          sep = ""
      )
      line <- readline()
    } else {
      cat("Length of partition = 1, nothing to look for...\n",
          "Press [enter] to continue\n",
          sep = ""
      )
      line <- readline()
    }

    # Look for 2 CPs (if du >= 3), 1 CP (if du = 2), or none (if du = 1)
    if (du == 1L) {
      k_star <- NULL
      i_star <- 1L
      h_star <- -Inf
    } else if (du == 2L) {
      tmp1 <- OSH2(lu, ru, x0, par, opts)
      k_star <- tmp1$k_star
      i_star <- tmp1$i_star
      h_star <- tmp1$h_star
    } else { # du >= 3
      tmp1 <- OSH2(lu, ru, x0, par, opts)
      tmp2 <- OSH3(lu, ru, x0, par, opts)
      if (tmp1$h_star >= tmp2$h_star) {
        k_star <- tmp1$k_star
        i_star <- tmp1$i_star
        h_star <- tmp1$h_star
      } else {
        k_star <- tmp2$k_star
        i_star <- tmp2$i_star
        h_star <- tmp2$h_star
      }
    }

    if (du > 1L) {
      if (length(k_star) == 0L) {
        cat("No maximum found (Does this ever happen?)\n",
            "Press [enter] to continue\n",
            sep = ""
        )
      } else {
        if (du == 2L) {
          tmp2$k_star <- NULL
        }
        display.1(
          xx.0, xx.Vit, xx, par, SS, NULL, CP, lu,
          tmp.plot1$res, tmp.plot1$t_star, tmp1$k_star,
          tmp.plot2$res, tmp.plot2$t_star, t(tmp2$k_star)
        )
        if (length(k_star) == 1L) {
          cat("One change point found.\n")
        } else {
          cat("Two change points found.\n")
        }
      }
      line <- readline()
    }

    # Best constant path for comparison
    tmp <- argH1(lu, ru, x0, par)
    h_star.const <- tmp$h_star
    i_star.const <- tmp$i_star

    # Update
    if (h_star > h_star.const + 1e-7) {
      SS <- rbind(
        SS[-(u:UU), ],
        cbind(c(lu, k_star), c(k_star - 1L, ru)),
        SS[-(1L:u), ]
      )
      zz <- c(zz[-(u:UU)], i_star, zz[-(1L:u)])
      UU <- length(zz)
    } else {
      zz[u] <- i_star.const
      u <- u + 1L
      k_star <- NULL
    }

    # For plotting purposes
    xx <- xx_SS.zz(SS, zz)
    display.0(xx.0, xx.Vit, xx, par, SS, k_star)
    if (!is.null(k_star)) {
      cat("New change point(s) added (green)\n")
    } else {
      cat("No new change point\n")
    }
    if (u <= UU) {
      cat("Next search interval: ", u, "\n",
          "Press [enter] to continue\n",
          sep = ""
      )
      line <- readline()
    }
  }

  xx <- xx_SS.zz(SS, zz)
  logp <- G0(xx, par)

  # For plotting purposes
  cat("Final fit\n")
  graphics::par(mfrow = c(1, 1), mar = c(4.2, 4.2, 0.2, 0.2))
  display.result(xx.0, xx.Vit, xx, par)
  return(list(
    xx = xx,
    logp = logp
  ))
}

#' Basic one-dimensional optimistic search
#'
#' Perform a ternary search over an interval \eqn{L}{L} to \eqn{R}{R} using an objective and an argument-returning function.
#'
#' @param L Integer scalar: left boundary (inclusive).
#' @param R Integer scalar: right boundary (inclusive), must satisfy R > L.
#' @param d0 Integer scalar ≥ 1: minimum search window width.
#' @param fun Function of one integer argument k; returns a numeric score.
#' @param argfun Function of one integer argument k; returns a list with at least:
#'   \describe{
#'     \item{h_star}{numeric score at k}
#'     \item{i_star}{integer or vector of associated state(s)}
#'   }
#' @param M Integer scalar: initial probe position; if < 0, computed automatically.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{k_star}{integer, best breakpoint(s)}
#'     \item{i_star}{integer or vector, associated state(s)}
#'     \item{h_star}{numeric, maximum score}
#'   }
#' @keywords internal
OS <- function(L, R, d0, fun, argfun, M = -1L) {
  # Initialize
  nu <- 1 / 2
  if (M < 0L) {
    M <- floor((L + nu * R) / (1 + nu))
  }

  # Determine HM
  HM <- fun(M)

  # Look for L and R
  while (R - L > d0) {
    test <- (R - M > M - L)
    W <- if (test) ceiling(R - nu*(R - M)) else ceiling(L + nu*(M - L))
    if (W == L || W == R) break
    HW <- fun(W)
    if (HW > HM) {
      if (test) L <- M else R <- M
      M  <- W
      HM <- HW
    } else {
      if (test) R <- W else L <- W
    }
  }

  # Full search
  k_star <- integer(0L)
  i_star <- integer(0L)
  h_star <- -Inf
  for (k in seq.int(L, R)) {
    tmp <- argfun(k)
    if (tmp$h_star > h_star) {
      k_star <- k
      i_star <- tmp$i_star
      h_star <- tmp$h_star
    }
  }
  return(list(
    k_star = k_star,
    i_star = i_star,
    h_star = h_star
  ))
}

#' 1-dimensional optimistic search
#'
#' @param l Left-most index
#' @param r Right-most index
#' @param x0 Previous state
#' @param par Model parameters
#' @param opts A list of options for the estimation
#'
#' @return Triplet
#' @keywords internal
OSH2 <- function(l, r, x0, par, opts) {
  # Apply base OS
  if (opts$rotate == TRUE && l + 1L < r) {
    H2_L <- H2(l, l + 1L, r, x0, par)
    H2_R <- H2(l, r, r, x0, par)
    a <- (H2_R - H2_L) / (r - l - 1L)
    res <- OS(
      l + 1L, r, opts$d0,
      function(k) (H2(l, k, r, x0, par) - a * k),
      function(k) argH2(l, k, r, x0, par)
    )
  } else {
    res <- OS(
      l + 1L, r, opts$d0,
      function(k) H2(l, k, r, x0, par),
      function(k) argH2(l, k, r, x0, par)
    )
  }
  return(res)
}

#' 2-dimensional optimistic search
#'
#' @param l Left-most index
#' @param r Right-most index
#' @param x0 Previous state
#' @param par Model parameters
#' @param opts A list of options for the estimation
#'
#' @return Triplet
#' @keywords internal
OSH3 <- function(l, r, x0, par, opts) {
  if (r - l > opts$d0) {
    # Prepare possible starting point(s) k0 for OSH3
    kk <-
      cbind(
        -1L,
        unique(
          sort(
            floor(l + 2L + (1L:opts$n.seeds) / (opts$n.seeds + 1L) * (r - l - 1L))
          )
        )
      )
    kk.n <- nrow(kk)

    # Resulting score(s)
    hh <- rep(-Inf, kk.n)

    # Loop over all starting point(s)
    for (d in 1L:kk.n) {
      tmp <- OSH3_k0(kk[d, ], l, r, x0, par, opts)
      kk[d, ] <- tmp$k_star
      hh[d] <- tmp$h_star
    }

    # Take the final fit with the largest score
    k_star <- kk[which.max(hh), ]
  } else {
    h_star <- -Inf
    k_star <- rep(NA, 2L)
    for (k1 in (l + 1L):(r - 1L)) {
      for (k2 in (k1 + 1L):r) {
        tmp <- H3(l, k1, k2, r, x0, par)
        if (tmp > h_star) {
          h_star <- tmp
          k_star <- c(k1, k2)
        }
      }
    }
  }

  # Compute final fit and return
  res <- argH3(l, k_star[1L], k_star[2L], r, x0, par)
  return(list(
    k_star = k_star,
    i_star = res$i_star,
    h_star = res$h_star
  ))
}

#' 2-dimensional optimistic search, with starting point k0
#'
#' @param k0 Starting point
#' @param l Left-most index
#' @param r Right-most index
#' @param x0 Previous state
#' @param par Model parameters
#' @param opts A list of options for the estimation
#'
#' @return Triplet
#' @keywords internal
OSH3_k0 <- function(k0, l, r, x0, par, opts) {
  # Bookkeeping
  h_new <- h_old <- -Inf

  # Initialize
  tau <- 0L # Whether a horizontal (0) or vertical (1) search is performed
  v <- 1L # Number of iterations
  v0 <- 20L # Maximum number of iterations

  # Start alternating horizontal and vertical searches
  while ( ((h_old < h_new) & (v < v0)) || v == 1L ) {
    # (0) Overwrite h_old
    h_old <- h_new
    # (1) Horizontal/Vertical search
    if (tau == 0L) {
      if (opts$rotate == TRUE && l + 1L < k0[2L] - 1L) {
        H3_L <- H3(l, l + 1L, k0[2L], r, x0, par)
        H3_R <- H3(l, k0[2L] - 1L, k0[2L], r, x0, par)
        a <- (H3_R - H3_L) / (k0[2L] - l - 2L)
        tmp <- OS(
          l + 1L, k0[2L] - 1L, opts$d0,
          function(k) (H3(l, k, k0[2L], r, x0, par) - a * k),
          function(k) argH3(l, k, k0[2L], r, x0, par),
          k0[1L]
        )
      } else {
        tmp <- OS(
          l + 1L, k0[2L] - 1L, opts$d0,
          function(k) H3(l, k, k0[2L], r, x0, par),
          function(k) argH3(l, k, k0[2L], r, x0, par),
          k0[1L]
        )
      }
    } else { # tau == 1
      if (opts$rotate == TRUE && k0[1L] + 1L < r) {
        H3_L <- H3(l, k0[1L], k0[1L] + 1L, r, x0, par)
        H3_R <- H3(l, k0[1L], r, r, x0, par)
        a <- (H3_R - H3_L) / (r - k0[1L] - 1L)
        tmp <- OS(
          k0[1L] + 1L, r, opts$d0,
          function(k) (H3(l, k0[1L], k, r, x0, par) - a * k),
          function(k) argH3(l, k0[1L], k, r, x0, par),
          k0[2L]
        )
      } else {
        tmp <- OS(
          k0[1L] + 1L, r, opts$d0,
          function(k) H3(l, k0[1L], k, r, x0, par),
          function(k) argH3(l, k0[1L], k, r, x0, par),
          k0[2L]
        )
      }
    }
    k0[1L + tau] <- tmp$k_star

    # (2) If on diagonal, apply OS to k \mapsto H^3(k, k+1) with k^*_1 as the
    #     initial probe
    if (k0[1L] + 1L == k0[2L]) {
      if (opts$rotate == TRUE && l + 2L < r) {
        H3_L <- H3(l, l + 1L, l + 2L, r, x0, par)
        H3_R <- H3(l, r - 1L, r, r, x0, par)
        a <- (H3_R - H3_L) / (r - l - 2L)
        tmp <- OS(
          l + 1L, r - 1L, opts$d0,
          function(k) (H3(l, k, k + 1L, r, x0, par) - a * k),
          function(k) argH3(l, k, k + 1L, r, x0, par),
          k0[1L]
        )
      } else {
        tmp <- OS(
          l + 1L, r - 1L, opts$d0,
          function(k) H3(l, k, k + 1L, r, x0, par),
          function(k) argH3(l, k, k + 1L, r, x0, par),
          k0[1L]
        )
      }
      k0 <- c(tmp$k_star, tmp$k_star + 1L)
    }

    # (3) Update h, tau and v
    h_new <- H3(l, k0[1L], k0[2L], r, x0, par)
    tau <- 1L - tau
    v <- v + 1L
  }

  # Compute final fit and return
  res <- argH3(l, k0[1L], k0[2L], r, x0, par)
  return(list(
    k_star = k0,
    i_star = res$i_star,
    h_star = res$h_star
  ))
}
