make_par <- function(seed = 1, n = 50, m = 3, stay = 0.9) {
  set.seed(seed)
  pp      <- matrix((1 - stay) / (m - 1), nrow = m, ncol = m)
  diag(pp) <- stay
  Pi <- rep(1/m, m)
  sample.HMM(n, m, Pi = Pi, pp = pp,
             emi.dist = "normal",
             emi.param = list(mu = seq_len(m), sigma = rep(1, m)))
}

test_that("G0 returns a finite scalar", {
  par <- make_par()
  lp  <- G0(par$xx, par)
  expect_length(lp, 1L)
  expect_true(is.finite(lp))
})

test_that("G0 agrees between R and C++ backends", {
  par     <- make_par(seed = 5)
  res_R   <- QATS.R(par)
  res_CPP <- QATS.CPP(par)
  # Both should give the same log-probability (same path)
  expect_equal(res_R$logp, res_CPP$logp, tolerance = 1e-10)
  expect_equal(res_R$xx,   res_CPP$xx)
})

test_that("argH1 / argH2 / argH3 return finite scores", {
  par <- make_par()
  l <- 1L; r <- 20L; x0 <- 1L

  h1 <- argH1(l, r, x0, par)
  expect_true(is.finite(h1$h_star))
  expect_true(h1$i_star %in% seq_len(par$m))

  h2 <- argH2(l, 10L, r, x0, par)
  expect_true(is.finite(h2$h_star))
  expect_equal(h2$k_star, 10L)

  h3 <- argH3(l, 7L, 15L, r, x0, par)
  expect_true(is.finite(h3$h_star))
  expect_equal(h3$k_star, c(7L, 15L))
})

test_that("argH2 h_star equals max G2 over states (argmax consistency)", {
  par <- make_par(seed = 3)
  l <- 1L; k1 <- 10L; r <- 20L; x0 <- 1L

  res <- argH2(l, k1, r, x0, par)
  # Verify h_star matches G2 at the returned states
  expect_equal(res$h_star, G2(l, k1, r, x0, res$i_star[1], res$i_star[2], par),
               tolerance = 1e-12)
  # Verify no other (i1, i2) pair has a higher G2
  best_brute <- max(vapply(par$mseq, function(i1)
    max(vapply(par$mseq[par$mseq != i1], function(i2)
      G2(l, k1, r, x0, i1, i2, par), numeric(1L))), numeric(1L)))
  expect_equal(res$h_star, best_brute, tolerance = 1e-12)
})
