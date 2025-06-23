# Required packages
rm(list = ls())
library(tidyverse)
library(DNAcopy)
library(RColorBrewer)
library(QATS)
theme_set(theme_bw())

# Data preparation
data(coriell)
cna <- CNA(cbind(coriell$Coriell.05296),
           coriell$Chromosome, coriell$Position,
           data.type = "logratio", sampleid = "c05296")

smo <- smooth.CNA(cna)
yy <- smo$c05296 %>% discard(is.na)
n <- length(yy)
m <- 3

# HMM initialization
pi <- rep(1/m, m)
A <- matrix(.01, m, m); diag(A) <- .99
A[2, 2] <- .98; A[1, 3] <- A[3, 1] <- 0

mu <- c(-0.4, 0, 0.4)
sigma <- rep(0.1, m)

B <- function(j, y) dnorm(y, mean = mu[j], sd = sigma[j])

# Forward and Backward
alpha <- matrix(0, m, n)
alpha[, 1] <- pi * B(1:m, yy[1])
for (t in 2:n) {
  alpha[, t] <- (alpha[, t-1] %*% A) * B(1:m, yy[t])
  alpha[, t] <- alpha[, t] / sum(alpha[, t])
}

beta <- matrix(0, m, n)
beta[, n] <- 1
for (t in (n-1):1) {
  beta[, t] <- (A %*% (B(1:m, yy[t+1]) * beta[, t+1]))
  beta[, t] <- beta[, t] / sum(beta[, t])
}

gamma <- alpha * beta
gamma <- gamma / rep(colSums(gamma), each = m)

# Visualization: Posterior state probabilities
colnames(gamma) <- paste0("t", seq_len(ncol(gamma)))
gamma %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "tt", values_to = "prob") %>%
  mutate(tt = rep(1:n, times = m), xx = rep(1:m, each = n)) %>%
  mutate(xx = paste("k =", xx)) %>%
  ggplot(aes(x = tt, y = prob)) +
  geom_line() +
  facet_grid(rows = vars(xx)) +
  labs(title = "Posterior probability", x = "Time", y = "P(Xt = k | Y)")

# Viterbi decoding
V <- matrix(0, m, n)
Ptr <- matrix(0, m, n)
V[, 1] <- pi * B(1:m, yy[1])
for (t in 2:n) {
  for (i in 1:m) {
    probs <- V[, t-1] * A[, i]
    Ptr[i, t] <- which.max(probs)
    V[i, t] <- max(probs) * B(i, yy[t])
  }
  V[, t] <- V[, t] / sum(V[, t])
}

xx.Vit <- numeric(n)
xx.Vit[n] <- which.max(V[, n])
for (t in (n-1):1) {
  xx.Vit[t] <- Ptr[xx.Vit[t+1], t+1]
}

# Visualization: Observed sequence + state
df_states <-
  tibble(tt = 1:n, yy = yy, xx.Vit = xx.Vit) %>%
  mutate(mu.Vit = mu[xx.Vit])

ggplot() +
  geom_point(data = df_states, aes(x = tt, y = yy), cex = 0.5) +
  geom_line(data = df_states, aes(x = tt, y = mu.Vit), col = "red", lwd = 2, alpha = 0.7) +
  labs(x = "Probe number", y = "Log intensity ratio")

# Re-estimation (EM Step)
xi <- array(0, dim = c(m, m, n-1))
for (t in 1:(n-1)) {
  denom <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      xi[i, j, t] <- alpha[i, t] * A[i, j] * B(j, yy[t+1]) * beta[j, t+1]
      denom <- denom + xi[i, j, t]
    }
  }
  xi[,,t] <- xi[,,t] / denom
}

pi.hat <- gamma[,1]
A.hat <- apply(xi, c(1,2), sum) / rowSums(gamma[, -n])
mu.hat <- rowSums(t(t(gamma) * yy)) / rowSums(gamma)
sigma.hat <- sapply(1:m, function(i) {
  sqrt(sum(gamma[i,] * (yy - mu.hat[i])^2) / sum(gamma[i,]))
})
B.hat <- function(j, y) dnorm(y, mean = mu.hat[j], sd = sigma.hat[j])

# Viterbi decoding
V <- matrix(0, m, n)
Ptr <- matrix(0, m, n)
V[, 1] <- pi.hat * B.hat(1:m, yy[1])
for (t in 2:n) {
  for (i in 1:m) {
    probs <- V[, t-1] * A.hat[, i]
    Ptr[i, t] <- which.max(probs)
    V[i, t] <- max(probs) * B.hat(i, yy[t])
  }
  V[, t] <- V[, t] / sum(V[, t])
}

xx.Vit.hat <- numeric(n)
xx.Vit.hat[n] <- which.max(V[, n])
for (t in (n-1):1) {
  xx.Vit.hat[t] <- Ptr[xx.Vit.hat[t+1], t+1]
}

# QATS decoding
QATS.par <- set.par(yy = yy, Pi = pi.hat, pp = A.hat,
                    emi.dist = "normal", emi.param = list(mu = mu.hat, sigma = sigma.hat))
QATS.res <- QATS.CPP(par = QATS.par, opts = list(n.seeds = 7))
xx.QATS.hat <- QATS.res$xx

# Combine all estimated paths
df_states_all_methods <-
  df_states %>%
  mutate(
    xx.Vit.hat = xx.Vit.hat,
    mu.Vit.hat = mu.hat[xx.Vit.hat],
    xx.QATS.hat = xx.QATS.hat,
    mu.QATS.hat = mu.hat[xx.QATS.hat]
  ) %>%
  pivot_longer(cols = xx.Vit:mu.QATS.hat) %>%
  mutate(
    method = case_when(
      grepl(".Vit", name) ~ "Viterbi",
      grepl(".QATS", name) ~ "QATS"
    ),
    method = factor(method, levels = c("QATS", "Viterbi")),
    parameters = case_when(
      !grepl(".hat", name) ~ "Initial",
      grepl(".hat", name) ~ "Estimated"
    ),
    variable = case_when(
      grepl("xx", name) ~ "xx",
      grepl("mu", name) ~ "mu"
    )
  ) %>%
  select(-name)

# Final Plot
cairo_pdf(filename = "Plot_aCGH.pdf", width = 7.5, height = 3)
ggplot() +
  geom_point(
    data = df_states_all_methods %>% select(tt, yy) %>% distinct(),
    mapping = aes(x = tt, y = yy),
    cex = 0.5
  ) +
  geom_line(
    data = df_states_all_methods %>% filter(variable == "mu", parameters == "Estimated"),
    mapping = aes(x = tt, y = value, col = method, lwd = method),
    alpha = 0.7
  ) +
  scale_linewidth_manual(values = c("Viterbi" = 0.5, "QATS" = 4)) +
  scale_color_manual(values = c("Viterbi" = "red", "QATS" = "orange")) +
  labs(x = "Probe number", y = "Log intensity ratio", col = "Method", lwd = "Method")
dev.off()
