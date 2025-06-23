rm(list = ls())
library(tidyverse)
library(QATS)

# ── 1. Parameters ───────────────────────────────────────────────────────────────

m      <- 5
mu     <- 1:m
sigma  <- rep(0.5, m)
n      <- 1e4 + 1
K      <- 7

# Build “true” HMM parameters
pp <- matrix(K / ((n - 1) * (m - 1)), nrow = m, ncol = m) %>%
  { diag(.) <- 1 - K/(n-1); . }

Pi <- rep(1/m, m)

# ── 2. Simulate base sequence ───────────────────────────────────────────────────

par_true <- sample.HMM(
  n        = n,
  m        = m,
  K        = K,
  emi.dist  = "normal",
  emi.param = list(mu = mu, sigma = sigma)
)

xx0 <- par_true$xx
yy  <- par_true$yy

# Quick check: number of state‐changes
nc <- sum(diff(xx0) != 0)
message("Number of changes in xx0: ", nc)

# Show the “true” fit
display.result(xx0, par = par_true, yy = yy)

# ── 3. Define & run all methods ────────────────────────────────────────────────

# A little helper to run QATS/Viterbi/G and K‐Seg
run_method <- function(label, par, ...) {
  res <- switch(
    label,
    Viterbi    = Viterbi.CPP(par),
    G_pMAP     = G_classifier.CPP(par, 1, 0, 0, 0),
    G_MPP      = G_classifier.CPP(par, 0, 0, 0, 1),
    G_MPM      = G_classifier.CPP(par, 0, 0, 1, 0),
    G_gV       = G_classifier.CPP(par, 0, 1, 0, 1),
    QATS_1     = QATS.CPP(par, list(n.seeds = 1, n.rep = 1e2)),
    QATS_5     = QATS.CPP(par, list(n.seeds = 5, n.rep = 1e2)),
    # for K-seg, '...' will supply Kmax
    K_seg      = K_segmentation.CPP(par, ...)
  )

  tibble(
    label = label,
    xx    = list(res$xx),
    time  = as.numeric(res$time)
  )
}

# List of fixed‐K methods
base_methods <- tribble(
  ~label,   ~par,       ~arg,
  "Viterbi", par_true,  NULL,
  "G_pMAP",  par_true,  NULL,
  "G_MPP",   par_true,  NULL,
  "G_MPM",   par_true,  NULL,
  "G_gV",    par_true,  NULL,
  "QATS_1",  par_true,  NULL,
  "QATS_5",  par_true,  NULL
)

# K‐segmentation for k = 1:(K+4)
K_max <- K + 4
kseg_tbl <- tibble(
  label = paste0("K_seg_", 1:K_max),
  par   = list(par_true),
  arg   = list(1:K_max)    # single vector, switch will use (...)
)

# combine
all_methods <- bind_rows(base_methods, kseg_tbl)

# run them all
results_tbl <- all_methods %>%
  mutate(
    out = pmap(list(label, par, arg), run_method)
  ) %>%
  select(out) %>%
  unnest(out)

# ── 4. Compute error metrics ───────────────────────────────────────────────────

metrics_tbl <- results_tbl %>%
  mutate(
    l0        = map_dbl(xx, ~ lp_norm(xx0, ., 0)),
    l1        = map_dbl(xx, ~ lp_norm(xx0, ., 1)),
    l2        = map_dbl(xx, ~ lp_norm(xx0, ., 2)),
    v_measure = map_dbl(xx, ~ V_measure(xx0, .))
  ) %>%
  select(label, time, l0, l1, l2, v_measure)

print(metrics_tbl)

# ── 5. Display two illustrative fits ──────────────────────────────────────────

# pull out the Viterbi & QATS_5 sequences
xx_vit  <- results_tbl %>% filter(label == "Viterbi") %>% pull(xx) %>% .[[1]]
xx_q5   <- results_tbl %>% filter(label == "QATS_5")  %>% pull(xx) %>% .[[1]]

display.result(xx0, xx_vit, xx_q5, par_true)

# ── 6. (Optional) step‐by‐step display ─────────────────────────────────────────

# For a smaller toy example...
toy <- sample.HMM(n = 3e2+1, m = 2, K = 10,
                  emi.dist  = "normal",
                  emi.param = list(mu = 1:2, sigma = rep(0.5, 2)))
par_toy <- set.par(toy$yy, toy$Pi, toy$pp,
                   emi.dist  = "normal",
                   emi.param = list(mu = 1:2, sigma = rep(0.5, 2)))

resV  <- Viterbi.CPP(par_toy)
resQ <- QATS.CPP(par_toy)

which(resV$xx != resQ$xx)

resQ  <- QATS.display(toy$xx, resV$xx, par_toy)
