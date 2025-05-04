rm(list = ls())
library(tidyverse)
library(QATS)

# ── 1. Parameters ───────────────────────────────────────────────────────────────

m         <- 5
sigma     <- 0.5
n         <- 1e6 + 1
K_true    <- 10
K_mis     <- 1e5
sigma_mis <- 2.0

# ── 2. Generate "true" and "misspecified" HMM parameter objects ───────────────

par_true <- sample.HMM(
  n       = n,
  m       = m,
  K       = K_true,
  emi.dist  = "normal",
  emi.param = list(sigma = sigma)
)

# start with a draw at the mis‐specified K/sigma
par_mis0 <- sample.HMM(
  n         = n,
  m         = m,
  K         = K_mis,
  emi.dist  = "normal",
  emi.param = list(sigma = sigma_mis)
)

# then randomly perturb its transition matrix and re‐set it to match yy & Pi
par_mis <- par_mis0$pp %>%
  {
    tmp <- matrix(runif(m * m), nrow = m) %>%
      sweep(1, rowSums(.), "/")
    (.) * tmp
  } %>%
  sweep(1, rowSums(.), "/") %>%
  set.par(
    yy        = par_true$yy,
    Pi        = par_true$Pi,
    pp        = .,                    # <–– placeholder forces this to be the pp argument
    emi.dist  = "normal",
    emi.param = list(sigma = sigma_mis)
  )

# ── 3. Define the grid of methods to run ───────────────────────────────────────

method_tbl <- tribble(
  ~label,            ~fun,        ~opts,                         ~par,
  "Viterbi",         Viterbi.CPP, NULL,                          par_true,
  "QATS_1",          QATS.CPP,    list(n.seeds = 1, n.rep = 1e2), par_true,
  "QATS_5",          QATS.CPP,    list(n.seeds = 5, n.rep = 1e2), par_true,
  "Viterbi_mis",     Viterbi.CPP, NULL,                          par_mis,
  "QATS_1_mis",      QATS.CPP,    list(n.seeds = 1, n.rep = 1e2), par_mis,
  "QATS_5_mis",      QATS.CPP,    list(n.seeds = 5, n.rep = 1e2), par_mis
)

# ── 4. Run each method, collect xx & timing ─────────────────────────────────

results_tbl <- method_tbl %>%
  mutate(
    result = pmap(
      list(f = fun, p = par, o = opts),
      function(f, p, o) {
        if (is.null(o)) {
          f(p)
        } else {
          f(p, o)
        }
      }
    ),
    xx   = map(result, "xx"),
    time = map_dbl(result, ~ as.numeric(.x$time))
  )

# ── 5. Compute error metrics in one go ────────────────────────────────────────

metrics_tbl <- results_tbl %>%
  mutate(
    l0        = map_dbl(xx, ~ lp_norm(par_true$xx, ., 0)),
    l1        = map_dbl(xx, ~ lp_norm(par_true$xx, ., 1)),
    l2        = map_dbl(xx, ~ lp_norm(par_true$xx, ., 2)),
    v_measure = map_dbl(xx, ~ V_measure(par_true$xx, .))
  ) %>%
  select(label, time, l0, l1, l2, v_measure)

print(metrics_tbl)

# ── 6. Display chosen sequences ────────────────────────────────────────────────

viterbi_seq <- results_tbl %>% filter(label == "Viterbi")   %>% pull(xx) %>% .[[1]]
qats5_seq   <- results_tbl %>% filter(label == "QATS_5")    %>% pull(xx) %>% .[[1]]

display.result(par_true$xx, viterbi_seq, qats5_seq, par_true)
