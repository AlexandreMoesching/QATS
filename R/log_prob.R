#' Log-probability of a hidden path
#'
#' @param xx Vector of hidden states
#' @param par Model parameters
#'
#' @return Scalar log-probability
#' @keywords internal
G0 <- function(xx, par) {
  kk <- which(diff(xx) != 0)
  SS <- cbind(c(1, kk + 1), c(kk, par$n))
  SS.n <- nrow(SS)

  # Initial state
  tmp0 <- par$logPi[xx[1]]

  # Transitions from i to i
  tmp1 <- sum((SS[, 2] - SS[, 1]) * diag(par$qq)[xx[SS[, 1]]])

  # Transitions from i to j
  tmp2 <- sum(par$qq[cbind(xx[SS[-SS.n, 2]], xx[SS[-1, 1]])])

  # Observations when at state i
  tmp3 <- rep(0, SS.n)
  for (k in 1:SS.n) {
    xk <- xx[SS[k, 1]]
    tmp3[k] <- par$GG[xk, SS[k, 2]] -
      ifelse(k == 1, 0, par$GG[xk, SS[k - 1, 2]])
  }
  tmp3 <- sum(tmp3)

  return(tmp0 + tmp1 + tmp2 + tmp3)
}

#' Score of a one-segment path on an interval
#'
#' @param l Left endpoint of the interval
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param i1 State of the segment
#' @param par Model parameters
#'
#' @return Scalar score
#' @keywords internal
G1 <- function(l, r, x0, i1, par) {
  qq <- par$qq
  GG <- par$GG
  tmp <- rep(0, 3)
  tmp[1] <- (l > 1) * qq[x0, i1] + (l == 1) * par$logPi[i1]
  tmp[2] <- (r - l) * qq[i1, i1]
  tmp[3] <- GG[i1, r] - ifelse(l > 1, GG[i1, l - 1], 0)
  return(sum(tmp))
}

#' Score of a two-segment path on an interval
#'
#' @param l Left endpoint of the interval
#' @param k1 Change point
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param i1 State of the first segment
#' @param i2 State of the second segment
#' @param par Model parameters
#'
#' @return Scalar score
#' @keywords internal
G2 <- function(l, k1, r, x0, i1, i2, par) {
  qq <- par$qq
  GG <- par$GG
  tmp <- rep(0, 6)
  tmp[1] <- (l > 1) * qq[x0, i1] + (l == 1) * par$logPi[i1]
  tmp[2] <- (k1 - l - 1) * qq[i1, i1]
  tmp[3] <- qq[i1, i2]
  tmp[4] <- (r - k1) * qq[i2, i2]
  tmp[5] <- GG[i1, k1 - 1] - ifelse(l > 1, GG[i1, l - 1], 0)
  tmp[6] <- GG[i2, r] - GG[i2, k1 - 1]
  return(sum(tmp))
}

#' Score of a three-segment path on an interval
#'
#' @param l Left endpoint of the interval
#' @param k1 First change point
#' @param k2 Second change point
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param i1 State of the first segment
#' @param i2 State of the second segment
#' @param i3 State of the third segment
#' @param par Model parameters
#'
#' @return Scalar score
#' @keywords internal
G3 <- function(l, k1, k2, r, x0, i1, i2, i3, par) {
  qq <- par$qq
  GG <- par$GG
  tmp <- rep(0, 9)
  tmp[1] <- (l > 1) * qq[x0, i1] + (l == 1) * par$logPi[i1]
  tmp[2] <- (k1 - l - 1) * qq[i1, i1]
  tmp[3] <- qq[i1, i2]
  tmp[4] <- (k2 - k1 - 1) * qq[i2, i2]
  tmp[5] <- qq[i2, i3]
  tmp[6] <- (r - k2) * qq[i3, i3]
  tmp[7] <- GG[i1, k1 - 1] - ifelse(l > 1, GG[i1, l - 1], 0)
  tmp[8] <- GG[i2, k2 - 1] - GG[i2, k1 - 1]
  tmp[9] <- GG[i3, r] - GG[i3, k2 - 1]
  return(sum(tmp))
}

#' Best one-segment score on an interval, maximised over states
#'
#' @param l Left endpoint of the interval
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Scalar best score
#' @keywords internal
H1 <- function(l, r, x0, par) {
  h_star <- -Inf
  for (i1 in par$mseq) {
    h_star.new <- G1(l, r, x0, i1, par)
    if (h_star.new > h_star) {
      h_star <- h_star.new
    }
  }
  return(h_star)
}

#' Best two-segment score on an interval, maximised over states
#'
#' @param l Left endpoint of the interval
#' @param k1 Change point
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Scalar best score
#' @keywords internal
H2 <- function(l, k1, r, x0, par) {
  mseq <- par$mseq
  h_star <- -Inf
  for (i1 in mseq) {
    for (i2 in mseq[-i1]) {
      h_star.new <- G2(l, k1, r, x0, i1, i2, par)
      if (h_star.new > h_star) {
        h_star <- h_star.new
      }
    }
  }
  return(h_star)
}

#' Best three-segment score on an interval, maximised over states
#'
#' @param l Left endpoint of the interval
#' @param k1 First change point
#' @param k2 Second change point
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Scalar best score
#' @keywords internal
H3 <- function(l, k1, k2, r, x0, par) {
  mseq <- par$mseq
  h_star <- -Inf
  for (i1 in mseq) {
    for (i2 in mseq[-i1]) {
      for (i3 in mseq[-i2]) {
        h_star.new <- G3(l, k1, k2, r, x0, i1, i2, i3, par)
        if (h_star.new > h_star) {
          h_star <- h_star.new
        }
      }
    }
  }
  return(h_star)
}

#' Dispatch to H1, H2, or H3 depending on the change-point configuration
#'
#' @param l Left endpoint of the interval
#' @param k1 First candidate change point
#' @param k2 Second candidate change point
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Scalar best score
#' @keywords internal
HD_wrap <- function(l, k1, k2, r, x0, par) {
  kk <- unique(sort(c(k1, k2)))
  if (length(kk) == 2) {
    if (l < kk[1]) {
      h_star <- H3(l, kk[1], kk[2], r, x0, par)
    } else {
      h_star <- H2(l, kk[2], r, x0, par)
    }
  } else {
    if (l < kk) {
      h_star <- H2(l, kk, r, x0, par)
    } else {
      h_star <- H1(l, r, x0, par)
    }
  }
  return(h_star)
}

#' Compute best path with one segment on the interval l:r
#'
#' @param l Left endpoint of the interval
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Change point (\code{NULL}), state value and H1-value
#' @export
argH1 <- function(l, r, x0, par) {
  k_star <- NULL
  i_star <- 1
  h_star <- -Inf
  for (i1 in par$mseq) {
    h_star.new <- G1(l, r, x0, i1, par)
    if (h_star.new > h_star) {
      k_star <- NULL
      i_star <- i1
      h_star <- h_star.new
    }
  }
  return(list(
    k_star = k_star,
    i_star = i_star,
    h_star = h_star
  ))
}

#' Compute best path with two segments and jump at k1 on the interval l:r
#'
#' @param l Left endpoint of the interval
#' @param k1 Change point
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Change point (\code{k1}), state values and H2-value
#' @export
argH2 <- function(l, k1, r, x0, par) {
  mseq <- par$mseq
  k_star <- NULL
  i_star <- 1
  h_star <- -Inf
  for (i1 in mseq) {
    for (i2 in mseq[-i1]) {
      h_star.new <- G2(l, k1, r, x0, i1, i2, par)
      if (h_star.new > h_star) {
        k_star <- k1
        i_star <- c(i1, i2)
        h_star <- h_star.new
      }
    }
  }
  return(list(
    k_star = k_star,
    i_star = i_star,
    h_star = h_star
  ))
}

#' Compute best path with three segments and jumps at k1, k2 on the interval l:r
#'
#' @param l Left endpoint of the interval
#' @param k1 First change point
#' @param k2 Second change point
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return Change points (\code{k1, k2}), state values and H3-value
#' @export
argH3 <- function(l, k1, k2, r, x0, par) {
  mseq <- par$mseq
  k_star <- NULL
  i_star <- 1
  h_star <- -Inf
  for (i1 in mseq) {
    for (i2 in mseq[-i1]) {
      for (i3 in mseq[-i2]) {
        h_star.new <- G3(l, k1, k2, r, x0, i1, i2, i3, par)
        if (h_star.new > h_star) {
          k_star <- c(k1, k2)
          i_star <- c(i1, i2, i3)
          h_star <- h_star.new
        }
      }
    }
  }
  return(list(
    k_star = k_star,
    i_star = i_star,
    h_star = h_star
  ))
}

#' Dispatch to argH1, argH2, or argH3 depending on the change-point
#' configuration
#'
#' @param l Left endpoint of the interval
#' @param k1 First candidate change point
#' @param k2 Second candidate change point
#' @param r Right endpoint of the interval
#' @param x0 Previous state
#' @param par Model parameters
#'
#' @return List with k_star, i_star, h_star
#' @keywords internal
argHD_wrap <- function(l, k1, k2, r, x0, par) {
  kk <- unique(sort(c(k1, k2)))
  if (length(kk) == 2) {
    if (l < kk[1]) {
      tmp <- argH3(l, kk[1], kk[2], r, x0, par)
    } else {
      tmp <- argH2(l, kk[2], r, x0, par)
    }
  } else {
    if (l < kk) {
      tmp <- argH2(l, kk, r, x0, par)
    } else {
      tmp <- argH1(l, r, x0, par)
    }
  }
  return(tmp)
}
