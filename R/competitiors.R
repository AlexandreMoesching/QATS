#' Viterbi decoder, R version
#'
#' @param par Model parameters
#'
#' @return Estimated sequence, its probability and estimation time
#' @export
Viterbi.R <- function(par) {
  # Load some variables
  m <- par$m
  mseq <- par$mseq
  n <- par$n
  logPi <- par$logPi
  qq <- par$qq
  g_mseq <- par$g_mseq

  # Start timer
  time <- Sys.time()

  # Start estimation
  zeta <- matrix(0, nrow = m, ncol = n)
  rho <- logPi + g_mseq[, 1]
  rho.new <- rho
  if (n > 1) {
    for (k in 2:n) {
      for (i in mseq) {
        tmp <- rho + qq[, i]
        zeta[i, k - 1] <- which.max(tmp)
        rho.new[i] <- max(tmp) + g_mseq[i, k]
      }
      rho <- rho.new
    }
  }
  zeta[, n] <- which.max(rho)
  xx <- rep(0, n)
  xx[n] <- zeta[1, n]
  if (n > 1) {
    for (k in (n - 1):1) {
      xx[k] <- zeta[xx[k + 1], k]
    }
  }

  # End timer
  time <- difftime(Sys.time(), time, units = "secs")

  # Compute path probability
  logp <- G0(xx, par)
  return(list(
    xx = xx,
    logp = logp,
    time = time
  ))
}

#' Viterbi decoder, C++ version
#'
#' @param par Model parameters
#' @param n.rep Number of repetition (for timing)
#'
#' @return Estimated sequence, its probability and estimation time
#' @export
Viterbi.CPP <- function(par, n.rep = 1) {
  # Load some variables
  m <- par$m
  n <- par$n
  logPi <- par$logPi
  qq <- par$qq
  g_mseq <- par$g_mseq

  # Start estimation
  res <- Viterbi_timer_cpp(n.rep, n, m, logPi, qq, g_mseq)

  # Compute path probability
  logp <- G0(c(res$xx), par)
  return(list(
    xx = c(res$xx),
    logp = logp,
    time = res$time,
    time_ms = res$time_ms
  ))
}

#' G-classifier, R version
#'
#' @param C1 Constant 1
#' @param C2 Constant 2
#' @param C3 Constant 3
#' @param C4 Constant 4
#' @param par Model parameters
#'
#' @return Estimated sequence, its probability and estimation time
#' @export
G_classifier.R <- function(par, C1 = 0, C2 = 1, C3 = 0, C4 = 0)
  # C1 = C3 = C4 = 0  <->  Viterbi
  # C2 = C3 = C4 = 0  <->  pointwise MAP
  # C1 = C2 = C3 = 0  <->  maximum prior probability
  # C1 = C2 = C4 = 0  <->  marginal prior mode
  # C1 = C3 = 0       <->  generalized Viterbi with suppressed contribution of data
{
  # Load some variables
  m <- par$m
  mseq <- par$mseq
  n <- par$n
  Pi <- par$Pi
  logPi <- par$logPi
  pp <- par$pp
  qq <- par$qq
  f_mseq <- par$f_mseq
  g_mseq <- par$g_mseq

  # Start timer
  time <- Sys.time()

  # Pre-computing some useful quantities
  logprob_x <- matrix(0, nrow = m, ncol = n)
  logprob_x[, 1] <- logPi
  prob_x <- Pi
  if (n > 1) {
    for (k in 2:n) {
      prob_x <- t(pp) %*% prob_x
      logprob_x[, k] <- log(prob_x)
    }
  }

  if (C1 > 0) {
    # >>>>>> (1) First pass forward
    cc <- rep(0, n)
    alpha <- matrix(0, nrow = m, ncol = n)
    tmp <- f_mseq[, 1] * Pi
    cc[1] <- sum(tmp)
    alpha[, 1] <- tmp / cc[1]
    if (n > 1) {
      for (k in 2:n) {
        for (i in mseq) {
          tmp[i] <- f_mseq[i, k] * sum(pp[, i] * alpha[, k - 1])
        }
        cc[k] <- sum(tmp)
        alpha[, k] <- tmp / cc[k]
      }
    }

    # <<<<<< (1) First pass backward
    tmp <- rep(0, m)
    beta_bar <- rep(1, m)
    alpha_bar <- matrix(0, nrow = m, ncol = n)
    alpha_bar[, n] <- alpha[, n]
    if (n > 1) {
      for (k in (n - 1):1) {
        for (i in mseq) {
          tmp[i] <- sum(pp[i, ] * f_mseq[, k] * beta_bar)
        }
        beta_bar <- tmp / cc[k + 1]
        alpha_bar[, k] <- beta_bar * alpha[, k]
      }
    }
  }

  # >>>>>> (2) Last pass forward
  C24 <- C2 + C4
  C234 <- C2 + C3 + C4
  tmp <- rep(0, m)
  zeta <- matrix(0, nrow = m, ncol = n)
  if (C1 > 0) {
    hh <- C1 * log(alpha_bar) + C2 * g_mseq + C3 * logprob_x
    rho <- C1 * log(alpha_bar[, 1]) + C2 * g_mseq[, 1] + C234 * logprob_x[, 1]
  } else {
    hh <- C2 * g_mseq + C3 * logprob_x
    rho <- C2 * g_mseq[, 1] + C234 * logprob_x[, 1]
  }
  rho.new <- rho
  if (n > 1) {
    for (k in 2:n) {
      for (i in mseq) {
        tmp <- rho + C24 * qq[, i]
        zeta[i, k - 1] <- which.max(tmp)
        rho.new[i] <- max(tmp) + hh[i, k]
      }
      rho <- rho.new
    }
  }
  zeta[, n] <- which.max(rho)

  # <<<<<< (2) Last pass backward
  xx <- rep(0, n)
  xx[n] <- zeta[1, n]
  if (n > 1) {
    for (k in (n - 1):1) {
      xx[k] <- zeta[xx[k + 1], k]
    }
  }

  # Stop timer
  time <- difftime(Sys.time(), time, units = "secs")

  # Compute path probability
  logp <- G0(xx, par)
  return(list(
    xx = xx,
    logp = logp,
    time = time
  ))
}

#' G-classifier, C++ version
#'
#' @param C1 Constant 1
#' @param C2 Constant 2
#' @param C3 Constant 3
#' @param C4 Constant 4
#' @param par Model parameters
#'
#' @return Estimated sequence, its probability and estimation time
#' @export
G_classifier.CPP <- function(par, C1 = 0, C2 = 1, C3 = 0, C4 = 0)
  # C1 = C3 = C4 = 0  <->  Viterbi
  # C2 = C3 = C4 = 0  <->  pointwise MAP
  # C1 = C2 = C3 = 0  <->  maximum prior probability
  # C1 = C2 = C4 = 0  <->  marginal prior mode
  # C1 = C3 = 0       <->  generalized Viterbi with suppressed contribution of data
{
  # Load some variables
  m <- par$m
  n <- par$n
  Pi <- par$Pi
  logPi <- par$logPi
  pp <- par$pp
  qq <- par$qq
  f_mseq <- par$f_mseq
  g_mseq <- par$g_mseq

  # Start timer
  time <- Sys.time()

  # Start estimation
  xx <- c(G_classifier_cpp(
    C1, C2, C3, C4, n, m,
    Pi, logPi, pp, qq, f_mseq, g_mseq
  ))

  # Stop timer
  time <- difftime(Sys.time(), time, units = "secs")

  # Compute path probability
  logp <- G0(xx, par)
  return(list(
    xx = xx,
    logp = logp,
    time = time
  ))
}

#' PMAP decoder, C++ version
#'
#' @param par Model parameters
#' @param n.rep Number of repetition (for timing)
#'
#' @return Estimated sequence, its probability and estimation time
#' @export
PMAP.CPP <- function(par, n.rep = 1) {
  # Load some variables
  m <- par$m
  n <- par$n
  Pi <- par$Pi
  pp <- par$pp
  f_mseq <- par$f_mseq

  # Start estimation
  res <- PMAP_timer_cpp(n.rep, n, m, Pi, pp, f_mseq)

  # Compute path probability
  logp <- G0(c(res$xx), par)
  return(list(
    xx = c(res$xx),
    logp = logp,
    time = res$time,
    time_ms = res$time_ms
  ))
}

#' K-segmentation, R version
#'
#' @param K_max Maximum number of constant pieces
#' @param par Model parameters
#'
#' @return Estimated sequences, their probabilities and estimation time
#' @export
K_segmentation.R <- function(par, K_max = min(10, n)) {
  # Load some variables
  m <- par$m
  mseq <- par$mseq
  n <- par$n
  logPi <- par$logPi
  qq <- par$qq
  g_mseq <- par$g_mseq

  # Start timer
  time <- Sys.time()

  # Pre-computing some useful quantities
  Kseq <- 1:K_max

  # Forward pass
  gamma <- array(-Inf, dim = c(m, n, K_max))
  gamma[, 1, 1] <- g_mseq[, 1] + logPi
  delta_x <- array(NA, dim = c(m, n, K_max))
  delta_s <- array(NA, dim = c(m, n, K_max))
  for (k in 2:n) {
    for (i in mseq) {
      for (s in Kseq) {
        curr.max <- -Inf
        i.max <- NA
        s.max <- NA
        tmp <- -Inf
        for (j in mseq) {
          t <- ifelse(j == i, s, s - 1)
          tmp <- ifelse(t >= 1, gamma[j, k - 1, t] + qq[j, i], -Inf)
          if (tmp > curr.max) {
            curr.max <- tmp
            i.max <- j
            s.max <- t
          }
        }
        gamma[i, k, s] <- curr.max + g_mseq[i, k]
        delta_x[i, k, s] <- i.max
        delta_s[i, k, s] <- s.max
      }
    }
  }

  # Backward pass
  xx <- matrix(NA, nrow = K_max, ncol = n)
  ss <- matrix(NA, nrow = K_max, ncol = n)
  for (K in Kseq) {
    xx[K, n] <- which.max(gamma[, n, K])
    ss[K, n] <- K
  }
  for (k in (n - 1):1) {
    for (K in Kseq) {
      xx[K, k] <- delta_x[xx[K, k + 1], k + 1, ss[K, k + 1]]
      ss[K, k] <- delta_s[xx[K, k + 1], k + 1, ss[K, k + 1]]
    }
  }

  # Stop timer
  time <- difftime(Sys.time(), time, units = "secs")

  # Check that all went well
  if (any(Kseq !=
          apply(xx, 1, function(ww) sum(ww[2:n] != ww[1:(n - 1)]) + 1))) {
    cat("Problem in K-segmentation.")
  }

  # Compute path probability
  logp <- rep(0, K_max)
  for (K in Kseq) {
    logp[K] <- G0(xx[K, ], par)
  }
  return(list(
    xx = xx,
    logp = logp,
    time = time
  ))
}

#' K-segmentation, C++ version
#'
#' @param K_max Maximum number of constant pieces
#' @param par Model parameters
#'
#' @return Estimated sequences, their probabilities and estimation time
#' @export
K_segmentation.CPP <- function(par, K_max = min(10, n)) {
  # Load some variables
  m <- par$m
  n <- par$n
  logPi <- par$logPi
  qq <- par$qq
  g_mseq <- par$g_mseq
  Kseq <- 1:K_max

  # Start timer
  time <- Sys.time()

  # Start estimation
  xx <- K_segmentation_cpp(K_max, n, m, logPi, qq, g_mseq)

  # Stop timer
  time <- difftime(Sys.time(), time, units = "secs")

  # Check that all went well
  if (any(Kseq !=
          apply(xx, 1, function(ww) sum(ww[2:n] != ww[1:(n - 1)]) + 1))) {
    cat("Problem in K-segmentation.\n")
  }

  # Compute path probability
  logp <- rep(0, K_max)
  for (K in Kseq) {
    logp[K] <- G0(xx[K, ], par)
  }
  return(list(
    xx = xx,
    logp = logp,
    time = time
  ))
}
