# demo_misspecified.R — Interactive demo of QATS under parameter misspecification.
#
# Compares Viterbi and QATS (1 and 5 seeds) on a single long sequence
# using both the true and misspecified parameters. Prints error metrics
# and displays the decoded paths.

library(QATS)
library(tidyverse)

# ── 1. Parameters ───────────────────────��────────────────────────���───────────

m         <- 5
mu        <- 1:m
sigma     <- rep(0.5, m)
n         <- 1e6 + 1
K_true    <- 10
K_mis     <- 1e5
mu_mis    <- 1:m
sigma_mis <- rep(2.0, m)

# ── 2. Generate true and misspecified HMM parameters ─────────────────��───────

par_true <- sample.HMM(
  n = n, m = m, K = K_true,
  emi.dist  = "normal",
  emi.param = list(mu = mu, sigma = sigma)
)

par_mis0 <- sample.HMM(
  n = n, m = m, K = K_mis,
  emi.dist  = "normal",
  emi.param = list(mu = mu_mis, sigma = sigma_mis)
)

nu <- 0.9
pp_perturbed <- par_mis0$pp *
  matrix(runif(m * m, min = nu, max = 1 / nu), nrow = m)
pp_perturbed <- sweep(pp_perturbed, 1, rowSums(pp_perturbed), "/")

par_mis <- set.par(
  yy        = par_true$yy,
  Pi        = par_true$Pi,
  pp        = pp_perturbed,
  emi.dist  = "normal",
  emi.param = list(mu = mu_mis, sigma = sigma_mis)
)

# ── 3. Run decoders ─────────────────────────────────────────��────────────────

methods <- tribble(
  ~label,          ~fun,        ~opts,                           ~par,
  "Viterbi",       Viterbi.CPP, NULL,                            par_true,
  "QATS_1",        QATS.CPP,    list(n.seeds = 1, n.rep = 100),  par_true,
  "QATS_5",        QATS.CPP,    list(n.seeds = 5, n.rep = 100),  par_true,
  "Viterbi_mis",   Viterbi.CPP, NULL,                            par_mis,
  "QATS_1_mis",    QATS.CPP,    list(n.seeds = 1, n.rep = 100),  par_mis,
  "QATS_5_mis",    QATS.CPP,    list(n.seeds = 5, n.rep = 100),  par_mis
)

results <- methods |>
  mutate(
    result = pmap(list(f = fun, p = par, o = opts), function(f, p, o) {
      if (is.null(o)) f(p) else f(p, o)
    }),
    xx   = map(result, "xx"),
    time = map_dbl(result, ~ as.numeric(.x$time))
  )

# ── 4. Compute error metrics ────────────────────────────────────────────────

metrics <- results |>
  mutate(
    l0        = map_dbl(xx, ~ lp_norm(par_true$xx, ., 0)),
    l1        = map_dbl(xx, ~ lp_norm(par_true$xx, ., 1)),
    l2        = map_dbl(xx, ~ lp_norm(par_true$xx, ., 2)),
    v_measure = map_dbl(xx, ~ V_measure(par_true$xx, .))
  ) |>
  select(label, time, l0, l1, l2, v_measure)

print(metrics)

# ── 5. Display ───────────────────��───────────────────────────────────────────

vit_xx   <- results |> filter(label == "Viterbi") |> pull(xx) |> pluck(1)
qats5_xx <- results |> filter(label == "QATS_5")  |> pull(xx) |> pluck(1)

display.result(par_true$xx, vit_xx, qats5_xx, par_true) +
  ggtitle("Truth (black), Viterbi (blue), QATS (red)")
