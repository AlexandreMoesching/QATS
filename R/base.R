## usethis namespace: start
#' @useDynLib QATS, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' Precomputes parameters for the estimation
#'
#' @param yy Sequence of observations
#' @param Pi Initial distribution
#' @param pp Transition matrix
#' @param emi.dist Emission distribution
#' @param emi.param Parameters for the emission distribution
#'
#' @return A list of computed parameters: log-probabilities and cumulative sums
#' @export
set.par <- function(yy, Pi, pp,
                    emi.dist = "normal", emi.param = list(sigma = 1)) {
  n <- length(yy)
  m <- length(Pi)

  stopifnot(emi.dist %in% c("normal", "g.student"))

  f_mseq <- matrix(0, nrow = m, ncol = n)
  g_mseq <- matrix(0, nrow = m, ncol = n)
  GG <- matrix(0, nrow = m, ncol = n)

  if (emi.dist == "normal") {
    stopifnot(!is.null(emi.param$sigma))
    for (i in 1:m) {
      f_mseq[i, ] <- stats::dnorm(yy, i, emi.param$sigma)
      g_mseq[i, ] <- stats::dnorm(yy, i, emi.param$sigma, log = TRUE)
      GG[i, ] <- cumsum(g_mseq[i, ])
    }
  } else if (emi.dist == "g.student") {
    stopifnot(!is.null(emi.param$nu) && !is.null(emi.param$sigma))
    for (i in 1:m) {
      f_mseq[i, ] <-
        stats::dt((yy - i) / emi.param$sigma, emi.param$nu) / emi.param$sigma
      g_mseq[i, ] <-
        stats::dt((yy - i) / emi.param$sigma, emi.param$nu, log = TRUE) -
        log(emi.param$sigma)
      GG[i, ] <- cumsum(g_mseq[i, ])
    }
  }
  return(list(
    yy = yy,
    n = n,
    m = m,
    mseq = 1:m,
    Pi = Pi,
    logPi = log(Pi),
    pp = pp,
    qq = log(pp),
    f_mseq = f_mseq,
    g_mseq = g_mseq,
    GG = GG,
    emi.dist = emi.dist,
    emi.param = emi.param
  ))
}

#' Generates an HMM process
#'
#' @param n Length of the sequence to generate
#' @param m Cardinality of the state space
#' @param Pi Initial state distribution
#' @param K Expected number of change points
#' @param pp Transition matrix
#' @param emi.dist Emission distribution
#' @param emi.param Parameters for the emission distribution
#'
#' @return An HMM process and all parameters for the estimation
#' @export
sample.HMM <- function(n, m, Pi = NULL, K = NULL, pp = NULL,
                       emi.dist = "normal", emi.param = list(sigma = 1)) {
  # Length of the sequence to generate
  stopifnot(n > 1)
  # Cardinality of the state space
  stopifnot(m > 1)
  # Initial state distribution
  stopifnot(length(Pi) %in% c(0, m))
  if (length(Pi) == 0) {
    Pi <- rep(1, m) / m
  } else {
    stopifnot(sum(Pi) == 1)
  }
  # Transition matrix
  if (!is.null(K)) {
    stopifnot(
      K < n - 1,
      is.null(pp)
    )
    pp <- matrix(K / ((n - 1) * (m - 1)), nrow = m, ncol = m)
    diag(pp) <- 1 - K / (n - 1)
  } else {
    stopifnot(
      !is.null(pp),
      !is.null(dim(pp)),
      dim(pp) == c(m, m),
      all(rowSums(pp) == 1)
    )
    SD <- eigen(t(pp) - diag(m))
    SS_prob <- abs(SD$vectors[, which(round(SD$values, digits = 13) == 0)])
    SS_prob <- SS_prob / sum(SS_prob) # Steady-state probabilities
    Prop_Stays_and_Jumps <- SS_prob * pp
    K <- (n - 1) * (1 - sum(diag(Prop_Stays_and_Jumps)))
  }
  # Emission distribution
  stopifnot(emi.dist %in% c("normal", "g.student"))
  if (emi.dist == "normal") {
    stopifnot(!is.null(emi.param$sigma))
    return(sample_norm_HMM_export_cpp(n, m, Pi, pp, emi.param$sigma))
  } else if (emi.dist == "g.student") {
    stopifnot(!is.null(emi.param$nu) && !is.null(emi.param$sigma))
    return(sample.t.HMM(n, m, Pi, pp, emi.param$nu, emi.param$sigma))
  }
}

#' Generates an HMM process with student noise and precomputes parameters for
#' estimation
#'
#' @param n Length of the sequence to generate
#' @param m Cardinality of the state space
#' @param Pi Initial state distribution
#' @param pp Transition matrix
#' @param nu Degrees of freedom
#' @param sigma Scale parameter
#'
#' @return The hidden sequence and the sequence of observations
#' @keywords internal
sample.t.HMM <- function(n, m, Pi, pp, nu, sigma) {
  # Declare logPi
  logPi <- log(Pi)

  # Declare qq
  qq <- log(pp)

  # Declare xx
  xx <- rep(0, n)
  mseq <- 1:m
  xx[1] <- sample(mseq, 1, prob = Pi)
  for (i in 2:n) {
    xx[i] <- sample(mseq, 1, prob = pp[xx[i - 1], ])
  }

  # Declare yy
  yy <- xx + sigma * stats::rt(n, nu)

  # Declare and compute f_mseq, g_mseq and GG
  f_mseq <- matrix(0, nrow = m, ncol = n)
  g_mseq <- matrix(0, nrow = m, ncol = n)
  GG <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    f_mseq[i, ] <- stats::dt((yy - i) / sigma, nu) / sigma
    g_mseq[i, ] <- stats::dt((yy - i) / sigma, nu, log = TRUE) - log(sigma)
    GG[i, ] <- cumsum(g_mseq[i, ])
  }

  # For information: Final standard deviation is
  StdDev <- sigma * sqrt(nu / (nu - 2))

  # Return
  return(list(
    n = n,
    m = m,
    mseq = mseq,
    Pi = Pi,
    logPi = logPi,
    pp = pp,
    qq = qq,
    xx = xx,
    yy = yy,
    f_mseq = f_mseq,
    g_mseq = g_mseq,
    GG = GG,
    StdDev = StdDev
  ))
}

#' Transform a sequence into the partition it creates and its levels
#'
#' @param xx Sequence of hidden states
#' @param n Length of the sequence
#'
#' @return Partition and corresponding levels
#' @keywords internal
SS.zz_xx <- function(xx, n) {
  SS <- c((1:(n - 1))[xx[1:(n - 1)] != xx[2:n]], n)
  UU <- length(SS)
  SS <- cbind(c(1, SS[-UU] + 1), SS)
  zz <- xx[SS[, 2]]
  return(list(SS = SS, zz = zz))
}

#' Transform a partition and its levels to a sequence
#'
#' @param SS Partition in the form of a 2-by-U matrix
#' @param zz Vector with the UU different levels
#'
#' @return Corresponding vector of hidden states
#' @keywords internal
xx_SS.zz <- function(SS, zz) {
  return(rep(zz, times = SS[, 2] - SS[, 1] + 1))
}

#' V-measure
#'
#' @param xx.0 Reference vector
#' @param xx.1 Estimated vector
#'
#' @return V-measure of xx.0 and xx.1
#' @export
V_measure <- function(xx.0, xx.1) {
  n <- length(xx.0)
  # a_{ij} is the number of elements of TRUE class i (i.e. X_k = i) that
  # were classified in class j (i.e. \hat{X}_k = j).
  aa <- table(xx.0, xx.1)
  a_io <- rowSums(aa)
  a_oj <- colSums(aa)

  H_CK <- -sum(aa * log(t(t(aa) / a_oj)), na.rm = TRUE) / n
  H_C <- -sum(a_io * log(a_io / n), na.rm = TRUE) / n
  H <- ifelse(abs(H_C) > 1e-16, 1 - H_CK / H_C, 1)

  H_KC <- -sum(aa * log(aa / a_io), na.rm = TRUE) / n
  H_K <- -sum(a_oj * log(a_oj / n), na.rm = TRUE) / n
  C <- ifelse(abs(H_K) > 1e-16, 1 - H_KC / H_K, 1)

  return(2 * H * C / (H + C))
}

#' lp norm
#'
#' @param xx.0 Reference vector
#' @param xx.1 Estimated vector
#' @param p Number indicating which l-norm should be computed
#'
#' @return lp-norm of xx.0 - xx.1
#' @export
lp_norm <- function(xx.0, xx.1, p) {
  if (p > 0) {
    return((sum(abs(xx.0 - xx.1)^p))^(1 / p) / length(xx.0))
  } else {
    return(sum(xx.1 != xx.0) / length(xx.0))
  }
}
