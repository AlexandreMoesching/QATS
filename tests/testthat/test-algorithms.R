make_par <- function(seed = 1, n = 100, m = 3, stay = 0.9, K = NULL, pp = NULL) {
  set.seed(seed)
  if (is.null(pp)) {
    pp      <- matrix((1 - stay) / (m - 1), nrow = m, ncol = m)
    diag(pp) <- stay
  }
  Pi <- rep(1/m, m)
  args <- list(n = n, m = m, Pi = Pi,
               emi.dist = "normal",
               emi.param = list(mu = seq_len(m) * 3, sigma = rep(1, m)))
  if (!is.null(K)) args$K  <- K
  else             args$pp <- pp
  do.call(sample.HMM, args)
}

# ── QATS ─────────────────────────────────────────────────────────────────────

test_that("QATS.R returns correct structure", {
  par <- make_par()
  res <- QATS.R(par)
  expect_named(res, c("xx", "logp", "time"))
  expect_length(res$xx, par$n)
  expect_true(all(res$xx %in% seq_len(par$m)))
  expect_true(is.finite(res$logp))
})

test_that("QATS.CPP returns correct structure", {
  par <- make_par()
  res <- QATS.CPP(par)
  expect_named(res, c("xx", "logp", "time", "time_ms"))
  expect_length(res$xx, par$n)
  expect_true(all(res$xx %in% seq_len(par$m)))
})

test_that("QATS.R and QATS.CPP produce the same path", {
  par     <- make_par(seed = 7)
  res_R   <- QATS.R(par)
  res_CPP <- QATS.CPP(par)
  expect_equal(res_R$xx, res_CPP$xx)
})

test_that("QATS path log-prob >= constant-path log-prob", {
  par    <- make_par(seed = 11)
  res_Q  <- QATS.R(par)
  # constant best path
  const  <- argH1(1L, par$n, 1L, par)
  xx_c   <- rep(const$i_star, par$n)
  lp_c   <- QATS:::G0(xx_c, par)
  expect_true(res_Q$logp >= lp_c - 1e-10)
})

test_that("QATS transitions are all admissible (positive probability)", {
  par <- make_par(seed = 99, n = 200, m = 4, K = 15)
  res <- QATS.CPP(par)
  xx  <- res$xx
  n   <- par$n
  transitions <- cbind(xx[-n], xx[-1L])
  # For every observed transition i->j, pp[i,j] must be > 0
  probs <- par$pp[transitions]
  expect_true(all(probs > 0))
})

# ── Viterbi ───────────────────────────────────────────────────────────────────

test_that("Viterbi.R and Viterbi.CPP produce the same path", {
  par   <- make_par(seed = 3)
  res_R <- Viterbi.R(par)
  res_C <- Viterbi.CPP(par)
  expect_equal(res_R$xx, res_C$xx)
})

test_that("Viterbi path is at least as good as QATS (Viterbi is globally optimal)", {
  # For small n, Viterbi should give a log-prob >= QATS
  par   <- make_par(seed = 21, n = 80)
  vit   <- Viterbi.CPP(par)
  qats  <- QATS.CPP(par)
  expect_true(vit$logp >= qats$logp - 1e-8)
})

# ── PMAP ─────────────────────────────────────────────────────────────────────

test_that("PMAP.CPP returns valid states", {
  par <- make_par(seed = 5)
  res <- PMAP.CPP(par)
  expect_length(res$xx, par$n)
  expect_true(all(res$xx %in% seq_len(par$m)))
})

# ── K-segmentation ────────────────────────────────────────────────────────────

test_that("K_segmentation.CPP returns K rows each with correct segment count", {
  par   <- make_par(seed = 8, n = 60)
  K_max <- 5L
  res   <- K_segmentation.CPP(par, K_max = K_max)
  expect_equal(nrow(res$xx), K_max)
  # Row K should have exactly K constant segments
  for (K in seq_len(K_max)) {
    n_segs <- sum(res$xx[K, -1] != res$xx[K, -par$n]) + 1L
    expect_equal(n_segs, K)
  }
})

test_that("K_segmentation.R and K_segmentation.CPP agree", {
  par   <- make_par(seed = 13, n = 40)
  K_max <- 4L
  res_R <- K_segmentation.R(par,   K_max = K_max)
  res_C <- K_segmentation.CPP(par, K_max = K_max)
  expect_equal(res_R$xx, res_C$xx)
})

# ── Accuracy smoke-test ───────────────────────────────────────────────────────

test_that("QATS achieves reasonable accuracy on well-separated states", {
  # Means 3 apart, sigma=1 => SNR high => error should be low
  par <- make_par(seed = 42, n = 500, m = 3, stay = 0.95)
  res <- QATS.CPP(par)
  err <- lp_norm(par$xx, res$xx, 0)
  expect_lt(err, 0.15)
})
