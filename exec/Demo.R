# _______________________________________________________________________________
####    Parameters                                                          ####
rm(list = ls())
library(rHMM)

# User selected parameters
m <- 5
sigma <- 0.5
n <- 1e4 + 1
K <- 7

# Parameters of the HMM
pp <- matrix(K / ((n - 1) * (m - 1)), nrow = m, ncol = m)
diag(pp) <- 1 - K / (n - 1)
Pi <- rep(1 / m, m)

# Generate a sequence
par <- sample.HMM(
  n = n, m = m, K = K,
  emi.dist = "normal",
  emi.param = list(sigma = sigma)
)
xx.0 <- par$xx
yy <- par$yy
length((1:(n - 1))[xx.0[1:(n - 1)] != xx.0[2:n]])

display.result(xx.0, par = par, yy = yy)

# _______________________________________________________________________________
####    Fits WITHOUT step-by-step display                                   ####
par(mfrow = c(1, 1), mar = c(4.2, 4.2, 0.2, 0.2))

# 1. Viterbi
res <- Viterbi.CPP(par)
xx <- matrix(res$xx, nrow = 1)
tt <- as.vector(res$time)
rownames(xx)[1] <- names(tt)[1] <- "Viterbi"

# 2. G-classifier: pointwise MAP
res <- G_classifier.CPP(par, 1, 0, 0, 0)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "G-pMAP"

# 3. G-classifier: maximum prior probability
res <- G_classifier.CPP(par, 0, 0, 0, 1)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "G-MPP"

# 4. G-classifier: marginal prior mode
res <- G_classifier.CPP(par, 0, 0, 1, 0)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "G-MPM"

# 5. G-classifier: generalized Viterbi
res <- G_classifier.CPP(par, 0, 1, 0, 1)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "G-gV"

# 6. K-segmentation
K_max <- K + 4
res <- K_segmentation.CPP(par, K_max)
xx <- rbind(xx, res$xx)
tt <- c(tt, rep(as.vector(res$time), K_max))
rownames(xx)[length(tt) - ((K_max - 1):0)] <-
  names(tt)[length(tt) - ((K_max - 1):0)] <-
  paste(rep("K-seg", K_max), 1:K_max)

# 7. QATS
opts <- list(n.seeds = 1, n.rep = 1e2)
res <- QATS.CPP(par, opts)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "QATS_1"

# 8. QATS
opts$n.seeds <- 5
res <- QATS.CPP(par, opts)
xx <- rbind(xx, res$xx)
tt <- c(tt, as.vector(res$time))
rownames(xx)[length(tt)] <- names(tt)[length(tt)] <- "QATS_5"

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

# Display some estimators
display.result(xx.0, xx["Viterbi", ], xx["QATS_5", ], par)
cbind(fit_eval, tt)

# _______________________________________________________________________________
####    Fits WITH step-by-step display                                      ####
n <- 3e2 + 1
m <- 2
K <- 10
sigma <- 0.5
tmp <- sample.HMM(
  n = n, m = m, K = K,
  emi.dist = "normal",
  emi.param = list(sigma = sigma)
)
xx.0 <- tmp$xx
yy <- tmp$yy
par <- set.par(xx.0, tmp$Pi, tmp$pp,
  emi.dist = "normal",
  emi.param = list(sigma = sigma)
)
# par <- set.par(yy, tmp$Pi, tmp$pp,
#                      emi.dist = "normal",
#                      emi.param = list(sigma = sigma))

res.Vit <- Viterbi.CPP(par)
xx <- res.Vit$xx
res.QATS <- QATS.display(xx.0, xx, par)
res.QATS <- QATS.CPP(par)

which(xx != res.QATS$xx)
