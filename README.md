
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QATS

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

# Standard deviation of the normal emission distributions
sigma <- 0.5

# Size of the observation sequence
n <- 1e3 + 1

# Expected number of change points
K <- 7

# Generate a sequence
par <- sample.HMM(
  n = n, m = m, K = K,
  emi.dist = "normal",
  emi.param = list(sigma = sigma)
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
#> [1] 5.9625e-05
```

To use the PMAP algorithm in order to estimate the hidden sequence,
simply do:

``` r
res <- PMAP.CPP(par)
xx.1.PMAP <- res$xx
(tt <- as.vector(res$time)) # Time in seconds
#> [1] 5.5667e-05
```

To use QATS, run:

``` r
res <- QATS.CPP(par)
xx.2 <- res$xx
(tt <- as.vector(res$time)) # Time in seconds
#> [1] 0.000569166
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
#> [1] 5 1 0
```

## Comparison of QATS, Viterbi, risk-based estimators and K-segmentation

``` r
# Set a larger observation sequence, the other parameters remain unchanged
n <- 1e6 + 1

# Generate a sequence
par <- sample.HMM(
  n = n, m = m, K = K,
  emi.dist = "normal",
  emi.param = list(sigma = sigma)
)
xx.0 <- par$xx
sum(diff(xx.0) != 0)
#> [1] 7

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
#>                    l0           l1           l2    V-Measure          tt
#> Viterbi  0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.057550875
#> pMAP-man 1.999998e-06 1.999998e-06 1.414212e-06 9.999776e-01 0.056712542
#> pMAP     9.999990e-07 9.999990e-07 9.999990e-07 9.999886e-01 0.239171982
#> MPP      1.000000e+00 1.563793e+00 1.675358e-03 0.000000e+00 0.177110910
#> MPM      6.677863e-01 1.899365e+00 2.367604e-03 1.859827e-06 0.158532858
#> gViterbi 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.181471825
#> K-seg 1  8.382852e-01 1.228223e+00 1.417073e-03 0.000000e+00 1.263044834
#> K-seg 2  6.499274e-01 8.515071e-01 1.120119e-03 3.074565e-01 1.263044834
#> K-seg 3  2.194378e-01 2.194378e-01 4.684416e-04 6.983234e-01 1.263044834
#> K-seg 4  2.007808e-01 2.007808e-01 4.480855e-04 7.183576e-01 1.263044834
#> K-seg 5  5.772294e-02 5.772294e-02 2.402559e-04 9.179990e-01 1.263044834
#> K-seg 6  3.906596e-02 3.906596e-02 1.976510e-04 9.240728e-01 1.263044834
#> K-seg 7  1.865698e-02 1.865698e-02 1.365905e-04 9.507386e-01 1.263044834
#> K-seg 8  0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 1.263044834
#> K-seg 9  9.999990e-07 9.999990e-07 9.999990e-07 9.999886e-01 1.263044834
#> K-seg 10 1.999998e-06 3.999996e-06 2.828424e-06 9.999801e-01 1.263044834
#> K-seg 11 7.999992e-06 7.999992e-06 2.828424e-06 9.999199e-01 1.263044834
#> QATS     0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.001497458
```

## Step-by-step computation of QATS

For that, make sure that the sample size is not too large ($n < 10^3$).

``` r
set.seed(1234)
n <- 3e2 + 1
par <- sample.HMM(
  n = n, m = m, K = K,
  emi.dist = "normal",
  emi.param = list(sigma = sigma)
)
res.Vit <- Viterbi.CPP(par)
res.QATS <- QATS.display(par$xx, res.Vit$xx, par)
```
