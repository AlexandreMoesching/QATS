
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QATS (Quick Adaptive Ternary Segmentation)

<!-- badges: start -->
<!-- badges: end -->

In hidden Markov models (HMM), one observes a noisy version of an
unobservable finite-state Markov chain. One of the main components of
HMM is to decode that signal, that is to estimate the sequence of hidden
states at the origin of the observation sequence. Existing decoding
algorithms such as [Viterbi
(1967)](https://doi.org/10.1109/TIT.1967.1054010) algorithm have
computational complexity at best linear in the size of the observed
sequence, and sub-quadratic in the size of the state space.

Assuming that the observation sequence is stored as specific cumulative
sums, we present Quick Adaptive Ternary Segmentation (QATS), a procedure
which decodes the hidden sequence in polylogarithmic computational
complexity in the size of the sequence, and cubic in the size of the
state space. In essence, the estimated sequence of states sequentially
maximizes local likelihood scores among all local paths with at most
three segments. The latter search is performed only approximately using
an adaptive search procedure. The resulting sequence is admissible in
the sense that all transitions occur with positive probability.

QATS and Viterbi algorithm, as well as algorithms to compute the
generalized risk-based estimator of [Lember and Koloydenko
(2014)](https://dl.acm.org/doi/10.5555/2627435.2627436) and the
K-segment approach of [Titsias, Holmes and Yau
(2016)](https://doi.org/10.1080/01621459.2014.998762) are implemented in
both **R** and **C++**.

Our implementation of QATS targets the homogeneous hidden‐Markov‐model
decoding problem, meaning that the transition matrix
$\boldsymbol{p}^{(k)} = \boldsymbol{p}$, for all $k=1,\ldots,n-1$,
remains fixed over time. However, our implementation imposes no
homogeneity requirement on the emission side: it only needs the
cumulative sums of the log‐emission densities at each time point (which
may vary arbitrarily). If you ever need a fully non‐homogeneous
HMM–where $\boldsymbol{p}^{(k)}$ changes with $k$–you can swap out each
use of the static $\boldsymbol{p}$ in the code for the time-indexed
$\boldsymbol{p}^{(k)}$ when computing transition log‐probabilities. This
should come at no resulting increase in computational time.

## Installation

You can install the development version of **QATS** from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("AlexandreMoesching/QATS")
```

## Usage

This is a basic example which shows how to solve a common problem:

``` r
# Load package
library(QATS)

# Set any seed, just for the example...
set.seed(123)

# Cardinality of the space of hidden states
m <- 5

# Means and standard deviations of the normal emission distributions
mu <- 1:m
sigma <- rep(0.5, m)

# Size of the observation sequence
n <- 1e3 + 1

# Expected number of change points
K <- 7

# Generate a sequence
par <- sample.HMM(
  n = n, m = m, K = K,
  emi.dist = "normal",
  emi.param = list(mu = mu, sigma = sigma)
)

# Extract the hidden (xx.0) and observed (yy) sequences
xx.0 <- par$xx
yy <- par$yy

# Actual number of change points of the generated sequence:
sum(diff(xx.0) != 0)
#> [1] 6
```

The function `sample.HMM()` already pre-computes useful quantities for
later such as the log-probability transition matrix, log-initial
distribution, densities and log-densities evaluated at the observations,
cumulative sums of log-densities, etc. The list of parameters will then
be passed to the different functions.

We can now plot the generated hidden sequence (black line) and
observations (grey points):

``` r
display.result(xx.0, par = par, yy = yy)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="90%" />

To use Viterbi algorithm in order to estimate the hidden sequence,
simply do:

``` r
res <- Viterbi.CPP(par)
xx.1 <- res$xx
(tt <- as.vector(res$time)) # Time in seconds
#> [1] 6.2959e-05
```

To use the PMAP algorithm in order to estimate the hidden sequence,
simply do:

``` r
res <- PMAP.CPP(par)
xx.1.PMAP <- res$xx
(tt <- as.vector(res$time)) # Time in seconds
#> [1] 7.0291e-05
```

To use QATS, run:

``` r
res <- QATS.CPP(par)
xx.2 <- res$xx
(tt <- as.vector(res$time)) # Time in seconds
#> [1] 0.000512708
```

Here is a plot of the true path (black), Viterbi-path (blue) and
QATS-path (red):

``` r
display.result(xx.0, xx.1, xx.2, par)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="90%" />

The following line shows that Viterbi and QATS produce the same path,
and that both did not estimate correctly 5 states out of the 1001.

``` r
c(sum(xx.0 != xx.1), sum(xx.1 != xx.1.PMAP), sum(xx.1 != xx.2))
#> [1] 6 2 0
```

## Comparison of QATS, Viterbi, risk-based estimators and K-segmentation

``` r
# Set a larger observation sequence, the other parameters remain unchanged
n <- 1e6 + 1

# Generate a sequence
par <- sample.HMM(
  n = n, m = m, K = K,
  emi.dist = "normal",
  emi.param = list(mu = mu, sigma = sigma)
)
xx.0 <- par$xx
sum(diff(xx.0) != 0)
#> [1] 12

# 1. Viterbi
res <- Viterbi.CPP(par)
xx <- matrix(res$xx, nrow = 1)
tt <- as.vector(res$time)
rownames(xx)[1] <- names(tt)[1] <- "Viterbi"

# 2. Pointwise MAP - manual
res <- PMAP.CPP(par)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "pMAP-man"

# 3. Pointwise MAP
res <- G_classifier.CPP(par, 1, 0, 0, 0)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "pMAP"

# 4. Maximum prior probability
res <- G_classifier.CPP(par, 0, 0, 0, 1)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "MPP"

# 5. Marginal prior mode
res <- G_classifier.CPP(par, 0, 0, 1, 0)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "MPM"

# 6. Generalized Viterbi
res <- G_classifier.CPP(par, 0, 1, 0, 1)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "gViterbi"

# 7. K-segmentation
K_max <- K + 4
res <- K_segmentation.CPP(par, K_max)
xx <- rbind(xx, res$xx)
tt <- c(tt, rep(as.vector(res$time), K_max))
rownames(xx)[length(tt) - ((K_max - 1):0)] <-
  names(tt)[length(tt) - ((K_max - 1):0)] <-
  paste(rep("K-seg", K_max), 1:K_max)

# 8. QATS
res <- QATS.CPP(par)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "QATS"

# Compute errors
n.estim <- nrow(xx)
fit_eval <- matrix(0, nrow = n.estim, ncol = 4)
colnames(fit_eval) <- c("l0", "l1", "l2", "V-Measure")

for (w in 1:n.estim) {
  fit_eval[w, 1] <- lp_norm(xx.0, xx[w, ], 0)
  fit_eval[w, 2] <- lp_norm(xx.0, xx[w, ], 1)
  fit_eval[w, 3] <- lp_norm(xx.0, xx[w, ], 2)
  fit_eval[w, 4] <- V_measure(xx.0, xx[w, ])
}
rownames(fit_eval) <- names(tt)

cbind(fit_eval, tt)
#>                    l0           l1           l2 V-Measure          tt
#> Viterbi  5.999994e-06 5.999994e-06 2.449487e-06 0.9999501 0.055746166
#> pMAP-man 5.999994e-06 5.999994e-06 2.449487e-06 0.9999501 0.061127459
#> pMAP     6.999993e-06 6.999993e-06 2.645749e-06 0.9999410 0.253034115
#> MPP      8.613451e-01 2.139248e+00 2.507768e-03 0.0000000 0.157251120
#> MPM      8.613451e-01 2.139248e+00 2.507768e-03 0.0000000 0.165350199
#> gViterbi 5.999994e-06 5.999994e-06 2.449487e-06 0.9999501 0.158349991
#> K-seg 1  5.979274e-01 9.759230e-01 1.316021e-03 0.0000000 1.263525963
#> K-seg 2  4.597595e-01 6.995873e-01 1.085929e-03 0.2491195 1.263525963
#> K-seg 3  6.750173e-01 7.692662e-01 9.786537e-04 0.3080688 1.263525963
#> K-seg 4  3.211047e-01 4.222776e-01 7.903308e-04 0.5280141 1.263525963
#> K-seg 5  2.476548e-01 2.753777e-01 5.751724e-04 0.6969048 1.263525963
#> K-seg 6  2.199318e-01 2.199318e-01 4.689686e-04 0.7686577 1.263525963
#> K-seg 7  1.715908e-01 1.715908e-01 4.142350e-04 0.7796504 1.263525963
#> K-seg 8  1.364359e-01 1.364359e-01 3.693721e-04 0.8115278 1.263525963
#> K-seg 9  9.011891e-02 9.011891e-02 3.001980e-04 0.8600035 1.263525963
#> K-seg 10 5.496395e-02 5.496395e-02 2.344438e-04 0.8910202 1.263525963
#> K-seg 11 4.427496e-02 4.427496e-02 2.104160e-04 0.9038261 1.263525963
#> QATS     5.999994e-06 5.999994e-06 2.449487e-06 0.9999501 0.002587375
```

## Step-by-step computation of QATS

For that, make sure that the sample size is not too large ($n < 10^3$).

``` r
set.seed(1234)
n <- 3e2 + 1
par <- sample.HMM(
  n = n, m = m, K = K,
  emi.dist = "normal",
  emi.param = list(mu = mu, sigma = sigma)
)
res.Vit <- Viterbi.CPP(par)
res.QATS <- QATS.display(par$xx, res.Vit$xx, par)
```
