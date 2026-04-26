# aCGH.R — Real-data illustration of QATS on array CGH data.
#
# Fits a 3-state Gaussian HMM to the smoothed Coriell sample c05296
# (from the DNAcopy package), runs one EM step to refine the parameters
# (Baum-Welch), then compares Viterbi and QATS decodings. Saves the
# final comparison to Plot_aCGH.pdf.
#
# Usage:
#   cd exec/real_data_example
#   Rscript aCGH.R

library(tidyverse)
library(DNAcopy)
library(QATS)
theme_set(theme_bw())

# --- 1. Data ------------------------------------------------------------------
# Sample c05296 from the Coriell cell line panel. CNA smoothing removes
# local noise before HMM fitting.

data(coriell)
cna <- CNA(cbind(coriell$Coriell.05296),
           coriell$Chromosome, coriell$Position,
           data.type = "logratio", sampleid = "c05296")
yy <- smooth.CNA(cna)$c05296 |> discard(is.na)
n  <- length(yy)

# --- 2. Initial parameters ----------------------------------------------------
# Three states: deletion / normal / gain. The transition matrix is nearly
# diagonal; direct deletion <-> gain transitions are structurally forbidden.

m     <- 3L
Pi    <- rep(1 / m, m)
A     <- matrix(0.01, m, m)
diag(A)  <- 0.99
A[2, 2]  <- 0.98
A[1, 3]  <- A[3, 1] <- 0

mu    <- c(-0.4, 0.0, 0.4)
sigma <- rep(0.1, m)

par_init <- set.par(yy, Pi, A,
                    emi.dist  = "normal",
                    emi.param = list(mu = mu, sigma = sigma))

# --- 3. Diagnostic: Viterbi decoding with initial parameters ------------------

vit_init <- Viterbi.CPP(par_init)

ggplot(tibble(t = seq_len(n), y = yy, mu_vit = mu[vit_init$xx]), aes(x = t)) +
  geom_point(aes(y = y), size = 0.5) +
  geom_line(aes(y = mu_vit), colour = "red", linewidth = 2, alpha = 0.7) +
  labs(x = "Probe number", y = "Log intensity ratio",
        title = "Viterbi decoding (initial parameters)")

# --- 4. One EM step to estimate parameters ------------------------------------
# Scaled forward-backward (Baum-Welch) to avoid numerical underflow.
# Emission probabilities and the transition matrix are taken from par_init
# to avoid recomputing them.

# Forward pass
alpha       <- matrix(0, m, n)
alpha[, 1]  <- par_init$Pi * par_init$f_mseq[, 1]
alpha[, 1]  <- alpha[, 1] / sum(alpha[, 1])
for (t in 2:n) {
  alpha[, t] <- (alpha[, t - 1] %*% par_init$pp) * par_init$f_mseq[, t]
  alpha[, t] <- alpha[, t] / sum(alpha[, t])
}

# Backward pass
beta       <- matrix(1, m, n)
for (t in (n - 1):1) {
  beta[, t] <- par_init$pp %*% (par_init$f_mseq[, t + 1] * beta[, t + 1])
  beta[, t] <- beta[, t] / sum(beta[, t])
}

# gamma[i, t] = P(X_t = i | Y_{1:n})
gamma <- alpha * beta
gamma <- sweep(gamma, 2, colSums(gamma), "/")

# Diagnostic: posterior state probabilities
as.data.frame(t(gamma)) |>
  setNames(paste("k =", seq_len(m))) |>
  mutate(t = seq_len(n)) |>
  pivot_longer(-t, names_to = "state", values_to = "prob") |>
  ggplot(aes(x = t, y = prob)) +
  geom_line() +
  facet_grid(rows = vars(state)) +
  labs(title = "Posterior state probabilities (initial parameters)",
        x = "Probe", y = expression(P(X[t] == k ~ "|" ~ Y)))

# xi[i, j, t] = P(X_t = i, X_{t+1} = j | Y_{1:n}), used for the A update
xi <- array(0, c(m, m, n - 1))
for (t in 1:(n - 1)) {
  mat        <- outer(alpha[, t], par_init$f_mseq[, t + 1] * beta[, t + 1]) * par_init$pp
  xi[, , t]  <- mat / sum(mat)
}

# M-step: closed-form parameter updates
Pi_hat    <- gamma[, 1]
A_hat     <- apply(xi, c(1, 2), sum) / rowSums(gamma[, -n])
mu_hat    <- drop(gamma %*% yy) / rowSums(gamma)
sigma_hat <- sapply(seq_len(m), function(i) {
  sqrt(sum(gamma[i, ] * (yy - mu_hat[i])^2) / sum(gamma[i, ]))
})

# --- 5. Decoding with estimated parameters ------------------------------------

par_hat  <- set.par(yy, Pi_hat, A_hat,
                    emi.dist  = "normal",
                    emi.param = list(mu = mu_hat, sigma = sigma_hat))

vit_hat  <- Viterbi.CPP(par_hat)
qats_hat <- QATS.CPP(par = par_hat, opts = list(n.seeds = 7))

# --- 6. Final comparison plot -------------------------------------------------

df_compare <- tibble(
  t       = seq_len(n),
  y       = yy,
  Viterbi = mu_hat[vit_hat$xx],
  QATS    = mu_hat[qats_hat$xx]
) |>
  pivot_longer(c(Viterbi, QATS), names_to = "method", values_to = "mu") |>
  mutate(method = factor(method, levels = c("QATS", "Viterbi")))

p_compare <- ggplot(df_compare, aes(x = t)) +
  geom_point(
    data    = distinct(df_compare, t, y),
    mapping = aes(y = y),
    size    = 0.3
  ) +
  geom_line(
    mapping = aes(y = mu, colour = method, linewidth = method),
    alpha   = 0.8
  ) +
  scale_linewidth_manual(values = c("Viterbi" = 0.5, "QATS" = 3)) +
  scale_colour_manual(values   = c("Viterbi" = "black", "QATS" = "grey")) +
  labs(x = "Probe number", y = "Log intensity ratio",
       colour = "Method", linewidth = "Method")

ggsave("Plot_aCGH.pdf", p_compare, width = 7.5, height = 3)
