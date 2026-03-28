test_that("lp_norm: p=0 counts mismatches", {
  x0 <- c(1, 1, 2, 2, 3)
  x1 <- c(1, 2, 2, 2, 3)
  expect_equal(lp_norm(x0, x1, 0), 1 / 5)
})

test_that("lp_norm: p=0 identical sequences gives 0", {
  x <- c(1, 2, 1, 3, 2)
  expect_equal(lp_norm(x, x, 0), 0)
})

test_that("lp_norm: p=1 agrees with manual calculation", {
  x0 <- c(1, 2, 3)
  x1 <- c(2, 2, 1)
  expect_equal(lp_norm(x0, x1, 1), (1 + 0 + 2) / 3)
})

test_that("V_measure: identical sequences gives 1", {
  x <- c(1, 1, 2, 2, 3, 3)
  expect_equal(V_measure(x, x), 1)
})

test_that("V_measure: output is in [0, 1]", {
  set.seed(42)
  x0 <- sample(1:3, 50, replace = TRUE)
  x1 <- sample(1:3, 50, replace = TRUE)
  v  <- V_measure(x0, x1)
  expect_true(v >= 0 && v <= 1)
})

test_that("SS.zz_xx and xx_SS.zz are inverses", {
  xx <- c(1, 1, 2, 2, 2, 3, 1, 1)
  n  <- length(xx)
  res <- SS.zz_xx(xx, n)
  xx2 <- xx_SS.zz(res$SS, res$zz)
  expect_equal(xx2, xx)
})
