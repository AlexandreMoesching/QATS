rm(list = ls())
library("QATS")

n <- 1e5 + 1; m <- 5; K <- 10; mu <- 1:m; sigma <- rep(0.5, m)

Pi <- rep(1, m) / m
pp <- matrix(K / ((n - 1) * (m - 1)), nrow = m, ncol = m)
diag(pp) <- 1 - K / (n - 1)

d0 <- 3; n_seeds <- 5; rotate <- FALSE; n_rep <- 1e1; n_sim <- 1e2

res <- QATS_vs_Viterbi_norm(n, m, Pi, pp, mu, sigma,
                            d0, n_seeds, rotate, n_rep, n_sim)

summary(res$res_Vit)
summary(res$res_QATS)

################################
rm(list = ls())
library("QATS")

n <- 1e5 + 1; m <- 5; K <- 10; mu <- 1:m; sigma <- rep(0.5, m)

Pi <- rep(1, m) / m
pp <- matrix(K / ((n - 1) * (m - 1)), nrow = m, ncol = m)
diag(pp) <- 1 - K / (n - 1)

d0 <- 3; n_seeds <- 5; rotate <- FALSE; n_rep <- 1e1; n_sim <- 1e2

res <- compare_norm(n, m, Pi, pp, mu, sigma,
                    d0, n_seeds, rotate, n_rep, n_sim)

summary(res$res_Vit)
summary(res$res_PMAP)
summary(res$res_QATS)

################################

rm(list = ls())
library("QATS")

n <- 1e5 + 1; m <- 5; K <- 10; mu <- 1:m; sigma <- rep(0.5, m)

Pi <- rep(1, m) / m
pp <- matrix(K / ((n - 1) * (m - 1)), nrow = m, ncol = m)
diag(pp) <- 1 - K / (n - 1)

d0 <- 3; n_seeds <- c(1,2,3,5); rotate <- FALSE; n_rep <- 1e1; n_sim <- 1e2

res <- QATS_nseeds_norm(n, m, Pi, pp, mu, sigma,
                        d0, n_seeds, rotate, n_rep, n_sim)

summary(res[, 1 + seq(0, 12, by = 4)])
summary(res[, 2 + seq(0, 12, by = 4)])
summary(res[, 3 + seq(0, 12, by = 4)])
summary(res[, 4 + seq(0, 12, by = 4)])
