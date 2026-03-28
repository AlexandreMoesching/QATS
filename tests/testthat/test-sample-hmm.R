make_pp <- function(m, stay = 0.9) {
  pp      <- matrix((1 - stay) / (m - 1), nrow = m, ncol = m)
  diag(pp) <- stay
  pp
}

test_that("sample.HMM returns correct structure for normal emissions", {
  set.seed(42)
  m  <- 3; n <- 100
  pp <- make_pp(m)
  Pi <- rep(1/m, m)
  par <- sample.HMM(n, m, Pi = Pi, pp = pp,
                    emi.dist = "normal",
                    emi.param = list(mu = 1:m, sigma = rep(1, m)))

  expect_equal(par$n, n)
  expect_equal(par$m, m)
  expect_length(par$xx, n)
  expect_length(par$yy, n)
  expect_true(all(par$xx %in% 1:m))
  expect_equal(dim(par$GG), c(m, n))
})

test_that("sample.HMM with K creates a valid transition matrix", {
  set.seed(7)
  m <- 3; n <- 200; K <- 10
  par <- sample.HMM(n, m, K = K,
                    emi.dist = "normal",
                    emi.param = list(mu = 1:m, sigma = rep(1, m)))
  # Each row of pp should sum to 1
  expect_equal(rowSums(par$pp), rep(1, m), tolerance = 1e-10)
})

test_that("sample.HMM rejects n = 1", {
  expect_error(sample.HMM(1, 2))
})

test_that("sample.HMM accepts Pi with 1/3 entries (FP tolerance)", {
  set.seed(3)
  m  <- 3; n <- 50
  Pi <- c(1/3, 1/3, 1/3)
  pp <- make_pp(m)
  expect_no_error(
    sample.HMM(n, m, Pi = Pi, pp = pp,
               emi.dist = "normal",
               emi.param = list(mu = 1:m, sigma = rep(1, m)))
  )
})
