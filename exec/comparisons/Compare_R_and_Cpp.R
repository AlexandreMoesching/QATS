rm(list = ls())
library(QATS)
library(ggplot2)

n <- 2e2 + 1; m <- 5; mu <- 1:m; sigma <- rep(0.5, m); K <- 10
opts <- list(n.rep = 1, rotate = FALSE)

n.sim <- 5e1
time_res <- matrix(0, n.sim, 3)
colnames(time_res) <- c("C++", "R", "R/C++")

for (i.sim in 1:n.sim) {
  cat("Simulation number", i.sim, "\n")
  par <- sample.HMM(n = n, m = m, K = K,
                    emi.dist = "normal",
                    emi.param = list(mu = mu, sigma = sigma))

  res_cpp <- QATS.CPP(par, opts)
  time_res[i.sim, 1] <- as.vector(res_cpp$time)

  res_R <- QATS.R(par, opts)
  time_res[i.sim, 2] <- as.vector(res_R$time)

  test <- all(res_cpp$xx == res_R$xx)
  if (!test) {
    cat("Not the same result!")
    break
  }
}
time_res <- time_res[time_res[, 1] != 0, ]
time_res[, 3] <- time_res[, 2] / time_res[, 1]

summary(time_res)

df <- data.frame(
  time   = c(time_res[, 1], time_res[, 2]),
  version = factor(rep(c("C++", "R"), each = nrow(time_res)),
                   levels = c("C++", "R"))
)

ggplot(df, aes(version, time)) +
  geom_boxplot(fill = "#3182bd", alpha = 0.6, width = 0.5) +
  scale_y_log10() +
  labs(x = NULL, y = "Time (seconds, log scale)",
       title = "QATS runtime: R vs C++") +
  theme_minimal(base_size = 13)
