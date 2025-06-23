rm(list = ls())
library("QATS")

n <- 2e2+1; m <- 5; mu <- 1:m; sigma <- rep(0.5, m); K <- 10
opts <- list(n.rep = 1, rotate = FALSE)

opts$SS <- NULL
# opts$SS <- seq(from = 1, to = n-1, by = floor(n/10))
# opts$SS <- cbind(opts$SS, c(opts$SS[-1]-1, n))

n.sim <- 5e1
time_res <- matrix(0, n.sim, 3)
colnames(time_res) <- c("C++", "R", "R/C++")

for (i.sim in 1:n.sim) {
  cat("Simulation number", i.sim, "\n")
  # Generate data
  par <- sample.HMM(n = n, m = m, K = K,
                    emi.dist = "normal",
                    emi.param = list(mu = mu, sigma = sigma))

  # Fit C++
  res_cpp <- QATS.CPP(par, opts)
  time_res[i.sim, 1] <- as.vector(res_cpp$time)

  # Fit R
  res_R <- QATS.R(par, opts)
  time_res[i.sim, 2] <- as.vector(res_R$time)

  # Test if equal
  test <- all(res_cpp$xx == res_R$xx)
  if (!test) {
    cat("Not the same result!")
    break
  }
}
time_res <- time_res[time_res[,1] != 0,]
time_res[,3] <- time_res[,2]/time_res[,1]

summary(time_res)
boxplot(time_res, log = "y")

# which(res_cpp$xx != res_R$xx)
# res_disp <-
#   QATS.display(par$xx, res_cpp$xx, par, opts)

# plot(1:n, par$xx, type = "l")
# lines(1:n, res_cpp$xx, col = 2)
# lines(1:n, res_R$xx, col = 3)
# lines(1:n, par$xx, col = 1, lty = 2)
# points(1:n, par$yy)
