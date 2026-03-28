make_pp <- function(m, stay = 0.9) {
  pp      <- matrix((1 - stay) / (m - 1), nrow = m, ncol = m)
  diag(pp) <- stay
  pp
}

test_that("set.par returns correct structure for normal emissions", {
  set.seed(1)
  m  <- 3; n <- 20
  pp <- make_pp(m)
  Pi <- rep(1/m, m)
  yy <- rnorm(n)

  par <- set.par(yy, Pi, pp, emi.dist = "normal",
                 emi.param = list(mu = 1:m, sigma = rep(1, m)))

  expect_equal(par$n, n)
  expect_equal(par$m, m)
  expect_equal(dim(par$GG),     c(m, n))
  expect_equal(dim(par$g_mseq), c(m, n))
  expect_equal(dim(par$f_mseq), c(m, n))
  # GG[i, k] must equal cumsum of g_mseq row i up to k
  expect_equal(par$GG[1, ], cumsum(par$g_mseq[1, ]))
})

test_that("set.par errors when mu is NULL (catches the old sigma×2 bug)", {
  m  <- 2; n <- 10
  pp <- make_pp(m)
  Pi <- c(0.5, 0.5)
  yy <- rnorm(n)
  expect_error(
    set.par(yy, Pi, pp, emi.dist = "normal",
            emi.param = list(sigma = rep(1, m)))  # mu missing — must error
  )
})

test_that("set.par accepts Pi summing to 1 within tolerance", {
  m  <- 3; n <- 10
  pp <- make_pp(m)
  Pi <- c(1/3, 1/3, 1/3)  # sum may differ from 1 by FP rounding
  yy <- rnorm(n)
  # Should not throw
  expect_no_error(
    set.par(yy, Pi, pp, emi.dist = "normal",
            emi.param = list(mu = 1:m, sigma = rep(1, m)))
  )
})
