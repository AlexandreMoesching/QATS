#' Quick Adaptive Ternary Segmentation, R version
#'
#' @param par Model parameters
#' @param opts A list of options for the estimation. Default is \code{opts =
#' list(n.seeds = 3, d0 = 3, SS = NULL, rotate = FALSE)} where \code{n.seeds} is
#' the number of seeds for the optimistic search, \code{d0} is the size of the
#' smallest search interval, \code{SS} is the initial partition and
#' \code{rotate} indicates whether or not the gain functions have to be rotated.
#'
#' @return Estimated sequence, its probability and estimation time
#' @export
QATS.R <- function(par,
                   opts = list(
                     n.seeds = 3,
                     d0 = 3,
                     SS = NULL,
                     rotate = FALSE
                   )) {
  # Set missing options
  if (is.null(opts$n.seeds)) {
    opts$n.seeds <- 3
  }
  if (is.null(opts$d0)) {
    opts$d0 <- 3
  }
  if (is.null(opts$SS)) {
    SS <- matrix(c(1, par$n), nrow = 1)
  } else {
    SS <- opts$SS
  }
  if (is.null(opts$rotate)) {
    opts$rotate <- FALSE
  }

  # Initialize
  UU <- nrow(SS)
  x0 <- 1
  u <- 1

  # Find best constant vector
  time <- Sys.time()
  tmp <- argH1(1, par$n, x0, par)
  zz <- rep(tmp$i_star, UU)

  # Loop
  while (u <= UU) {
    # Current part of the partition under investigation
    lu <- SS[u, 1]
    ru <- SS[u, 2]
    du <- ru - lu + 1

    # Last value of the path before the current partition
    x0 <- ifelse(u > 1, zz[u - 1], 1)

    # Look for 2 CPs (if du >= 3), 1 CP (if du = 2), or none (if du = 1)
    if (du == 1) {
      k_star <- NULL
      i_star <- 1
      h_star <- -Inf
    } else if (du == 2) {
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
        cbind(c(lu, k_star), c(k_star - 1, ru)),
        SS[-(1:u), ]
      )
      zz <- c(zz[-(u:UU)], i_star, zz[-(1:u)])
      UU <- length(zz)
    } else {
      zz[u] <- i_star.const
      u <- u + 1
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

#' Quick Adaptive Ternary Segmentation, C++ version
#'
#' @param par Model parameters
#' @param opts A list of options for the estimation. Default is \code{opts =
#' list(n.seeds = 3, d0 = 3, rotate = FALSE, n.rep = 1)} where
#' \code{n.seeds} is the number of seeds for the optimistic search,
#' \code{d0} is the size of the smallest search interval,
#' \code{SS} is the initial partition,
#' \code{rotate} indicates whether or not the gain functions have to be rotated,
#' \code{n.rep} is the umber of repetitions (for timing).
#'
#' @return Estimated sequence, its probability and estimation time
#' @export
QATS.CPP <- function(par,
                     opts = list(
                       n.seeds = 3,
                       d0 = 3,
                       SS = NULL,
                       rotate = FALSE,
                       n.rep = 1
                     )) {
  # Set missing options
  if (is.null(opts$n.seeds)) {
    opts$n.seeds <- 3
  }
  if (is.null(opts$d0)) {
    opts$d0 <- 3
  }
  if (is.null(opts$SS)) {
    SS <- 0
    UU <- 1
  } else {
    SS <- opts$SS[, 1] - 1
    UU <- length(SS)
  }
  if (is.null(opts$rotate)) {
    opts$rotate <- FALSE
  }
  if (is.null(opts$n.rep)) {
    opts$n.rep <- 1
  }

  # Load some variables
  n <- par$n
  m <- par$m
  logPi <- par$logPi
  qq <- par$qq
  GG <- par$GG

  # Start estimation
  res <- QATS_timer_cpp(
    opts$d0, opts$n.seeds, opts$rotate, opts$n.rep,
    n, m, logPi, qq, GG, SS, UU
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

#' Quick Adaptive Ternary Segmentation, with display
#'
#' @param xx.0 True hidden sequence
#' @param xx.Vit Hidden sequence estimated by Viterbi
#' @param par Model parameters
#' @param opts A list of options for the estimation. Default is \code{opts =
#' list(n.seeds = 3, d0 = 3, SS = NULL, rotate = FALSE)} where \code{n.seeds} is
#' the number of seeds for the optimistic search, \code{d0} is the size of the
#' smallest search interval, \code{SS} is the initial partition and
#' \code{rotate} indicates whether or not the gain functions have to be rotated.
#'
#' @return Estimated sequence and its probability
#' @export
QATS.display <- function(xx.0, xx.Vit, par,
                         opts = list(
                           n.seeds = 3,
                           d0 = 3,
                           SS = NULL,
                           rotate = FALSE
                         )) {
  # Set missing options
  if (is.null(opts$n.seeds)) {
    opts$n.seeds <- 3
  }
  if (is.null(opts$d0)) {
    opts$d0 <- 3
  }
  if (is.null(opts$SS)) {
    SS <- matrix(c(1, par$n), nrow = 1)
  } else {
    SS <- opts$SS
  }
  if (is.null(opts$rotate)) {
    opts$rotate <- FALSE
  }

  # Initialize
  UU <- nrow(SS)
  x0 <- 1
  u <- 1

  # Number of change points
  n <- par$n
  CP <- (2:n)[xx.0[1:(n - 1)] != xx.0[2:n]]

  # Find best constant vector
  tmp <- argH1(1, par$n, x0, par)
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
    lu <- SS[u, 1]
    ru <- SS[u, 2]
    Su <- lu:ru
    du <- length(Su)

    # Last value of the path before the current partition
    x0 <- ifelse(u > 1, zz[u - 1], 1)

    # Plot the functions to maximize
    if (du > 1) {
      tmp.plot1 <- H2_vec(Su, du, x0, par)
      if (du > 2) {
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
    if (du == 1) {
      k_star <- NULL
      i_star <- 1
      h_star <- -Inf
    } else if (du == 2) {
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

    if (du > 1) {
      if (length(k_star) == 0) {
        cat("No maximum found (Does this ever happen?)\n",
            "Press [enter] to continue\n",
            sep = ""
        )
      } else {
        if (du == 2) {
          tmp2$k_star <- NULL
        }
        display.1(
          xx.0, xx.Vit, xx, par, SS, NULL, CP, lu,
          tmp.plot1$res, tmp.plot1$t_star, tmp1$k_star,
          tmp.plot2$res, tmp.plot2$t_star, t(tmp2$k_star)
        )
        if (length(k_star) == 1) {
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
        cbind(c(lu, k_star), c(k_star - 1, ru)),
        SS[-(1:u), ]
      )
      zz <- c(zz[-(u:UU)], i_star, zz[-(1:u)])
      UU <- length(zz)
    } else {
      zz[u] <- i_star.const
      u <- u + 1
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
#' @param L Initial left index
#' @param R Initial right index
#' @param d0 Smallest search interval
#' @param fun Function to optimize
#' @param argfun Function to optimize which also returns arguments
#' @param M Initial probe (optional)
#'
#' @return Triplet
#' @keywords internal
OS <- function(L, R, d0, fun, argfun, M = -1) {
  # Initialize
  nu <- 1 / 2
  if (M == -1) {
    M <- floor((L + nu * R) / (1 + nu))
  }

  # Determine HM
  HM <- fun(M)

  # Look for L and R
  while (R - L > d0) {
    test <- (R - M > M - L)
    W <- ifelse(test, ceiling(R - nu * (R - M)), ceiling(L + nu * (M - L)))
    if (W %in% c(L, R)) break

    # Determine HW
    HW <- fun(W)

    # Update
    if (HW > HM) {
      ifelse(test, L <- M, R <- M)
      M <- W
      HM <- HW
    } else {
      ifelse(test, R <- W, L <- W)
    }
  }

  # Full search
  k_star <- NULL
  i_star <- 1
  h_star <- -Inf
  for (k in L:R) {
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
  if (opts$rotate == TRUE && l + 1 < r) {
    H2_L <- H2(l, l + 1, r, x0, par)
    H2_R <- H2(l, r, r, x0, par)
    a <- (H2_R - H2_L) / (r - l - 1)
    res <- OS(
      l + 1, r, opts$d0,
      function(k) (H2(l, k, r, x0, par) - a * k),
      function(k) argH2(l, k, r, x0, par)
    )
  } else {
    res <- OS(
      l + 1, r, opts$d0,
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
        -1,
        unique(
          sort(
            floor(l + 2 + (1:opts$n.seeds) / (opts$n.seeds + 1) * (r - l - 1))
          )
        )
      )
    kk.n <- nrow(kk)

    # Resulting score(s)
    hh <- rep(-Inf, kk.n)

    # Loop over all starting point(s)
    for (d in 1:kk.n) {
      tmp <- OSH3_k0(kk[d, ], l, r, x0, par, opts)
      kk[d, ] <- tmp$k_star
      hh[d] <- tmp$h_star
    }

    # Take the final fit with the largest score
    k_star <- kk[which.max(hh), ]
  } else {
    h_star <- -Inf
    k_star <- rep(NA, 2)
    for (k1 in (l + 1):(r - 1)) {
      for (k2 in (k1 + 1):r) {
        tmp <- H3(l, k1, k2, r, x0, par)
        if (tmp > h_star) {
          h_star <- tmp
          k_star <- c(k1, k2)
        }
      }
    }
  }

  # Compute final fit and return
  res <- argH3(l, k_star[1], k_star[2], r, x0, par)
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
  tau <- 0 # Whether a horizontal (0) or vertical (1) search is performed
  v <- 1 # Number of iterations
  v0 <- 20 # Maximum number of iterations

  # Start alternating horizontal and vertical searches
  while ( ((h_old < h_new) & (v < v0)) || v == 1 ) {
    # (0) Overwrite h_old
    h_old <- h_new
    # (1) Horizontal/Vertical search
    if (tau == 0) {
      if (opts$rotate == TRUE && l + 1 < k0[2] - 1) {
        H3_L <- H3(l, l + 1, k0[2], r, x0, par)
        H3_R <- H3(l, k0[2] - 1, k0[2], r, x0, par)
        a <- (H3_R - H3_L) / (k0[2] - l - 2)
        tmp <- OS(
          l + 1, k0[2] - 1, opts$d0,
          function(k) (H3(l, k, k0[2], r, x0, par) - a * k),
          function(k) argH3(l, k, k0[2], r, x0, par),
          k0[1]
        )
      } else {
        tmp <- OS(
          l + 1, k0[2] - 1, opts$d0,
          function(k) H3(l, k, k0[2], r, x0, par),
          function(k) argH3(l, k, k0[2], r, x0, par),
          k0[1]
        )
      }
    } else { # tau == 1
      if (opts$rotate == TRUE && k0[1] + 1 < r) {
        H3_L <- H3(l, k0[1], k0[1] + 1, r, x0, par)
        H3_R <- H3(l, k0[1], r, r, x0, par)
        a <- (H3_R - H3_L) / (r - k0[1] - 1)
        tmp <- OS(
          k0[1] + 1, r, opts$d0,
          function(k) (H3(l, k0[1], k, r, x0, par) - a * k),
          function(k) argH3(l, k0[1], k, r, x0, par),
          k0[2]
        )
      } else {
        tmp <- OS(
          k0[1] + 1, r, opts$d0,
          function(k) H3(l, k0[1], k, r, x0, par),
          function(k) argH3(l, k0[1], k, r, x0, par),
          k0[2]
        )
      }
    }
    k0[1 + tau] <- tmp$k_star

    # (2) If on diagonal, apply OS to k \mapsto H^3(k, k+1) with k^*_1 as the
    #     initial probe
    if (k0[1] + 1 == k0[2]) {
      if (opts$rotate == TRUE && l + 2 < r) {
        H3_L <- H3(l, l + 1, l + 2, r, x0, par)
        H3_R <- H3(l, r - 1, r, r, x0, par)
        a <- (H3_R - H3_L) / (r - l - 2)
        tmp <- OS(
          l + 1, r - 1, opts$d0,
          function(k) (H3(l, k, k + 1, r, x0, par) - a * k),
          function(k) argH3(l, k, k + 1, r, x0, par),
          k0[1]
        )
      } else {
        tmp <- OS(
          l + 1, r - 1, opts$d0,
          function(k) H3(l, k, k + 1, r, x0, par),
          function(k) argH3(l, k, k + 1, r, x0, par),
          k0[1]
        )
      }
      k0 <- c(tmp$k_star, tmp$k_star + 1)
    }

    # (3) Update h, tau and v
    h_new <- H3(l, k0[1], k0[2], r, x0, par)
    tau <- 1 - tau
    v <- v + 1
  }

  # Compute final fit and return
  res <- argH3(l, k0[1], k0[2], r, x0, par)
  return(list(
    k_star = k0,
    i_star = res$i_star,
    h_star = res$h_star
  ))
}
