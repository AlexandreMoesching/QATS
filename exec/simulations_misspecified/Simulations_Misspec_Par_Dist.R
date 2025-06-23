rm(list = ls())
library(QATS)
library(doSNOW)
library(tidyverse)

# cores <- parallel::detectCores()
cores <- 7                  # Number of cores used on the server
n.sim <- 1e2                # Number of simulations for a given setup
nu_vals <- c(1, 2, 5, 10, 15, 20)

# Prepare parallel computing
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)

# All results
all_res <- NULL

for (m in 2) { # c(2, 3, 5)
  for (s in 1) { # c(0.1, 1.0)
    for (nu in nu_vals) {
      # Message
      cat("--------------------------------------------------------\n")
      cat("Setting m = ", m, ", sigma = ", s, ", nu = ", nu, "\n", sep = "")
      start <- Sys.time()

      # HMM parameters
      Pi <- rep(1, m) / m

      # Prepare randomized simulations
      nKpp <- matrix(NA, nrow = 0, ncol = 4,
                     dimnames = list(NULL, c("n", "K", "pp_offdiag", "pp_diag")))
      for (n.pow in 5) { # 3:6
        n <- 10^n.pow + 1
        K.seq <- c(outer(c(1,2,5), 10^(0:7))); K.seq <- K.seq[K.seq < n/5e1]
        nKpp <- rbind(nKpp, cbind(n, K.seq, K.seq / ((n - 1) * (m - 1)), 1 - K.seq / (n - 1)))
      }
      nKpp <- do.call(rbind, replicate(n.sim, nKpp, simplify = FALSE))
      nKpp.n <- nrow(nKpp)
      nKpp <- nKpp[sample(1:nKpp.n, nKpp.n, replace = FALSE), ]

      # Prepare progress bar
      pb <- txtProgressBar(min = 1, max = nKpp.n, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      cl_opts <- list(progress = progress)

      # Start simulations
      result <- foreach(i = 1:nKpp.n,
                        .options.snow = cl_opts,
                        .combine = "rbind",
                        .packages = c("tidyverse", "QATS")) %dopar% {
                          # Setup
                          n <- nKpp[i, 1]
                          K <- nKpp[i, 2]

                          # HMM parameters
                          pp <- matrix(nKpp[i, 3], nrow = m, ncol = m)
                          diag(pp) <- nKpp[i, 4]

                          # Misspecify parameters
                          # Pi_mis <- Pi * runif(m, 1 / nu, nu)
                          # Pi_mis <- Pi_mis / sum(Pi_mis)
                          Pi_mis <- Pi

                          perturb <- matrix(runif(m * m, 1 / nu, nu), nrow = m)
                          diag(perturb) <- 1
                          pp_mis <- pp * perturb
                          pp_mis <- sweep(pp_mis, 1, rowSums(pp_mis), '/')
                          # pp_mis <- pp

                          # s_mis <- s * runif(1, 4 / (3 + nu), (3 + nu) / 4)
                          s_mis <- s

                          # Combine results
                          bind_rows(
                            tibble(
                              Parameter = "Pi",
                              `Row index` = 1:m,
                              `Col index` = 1,
                              `Well specified value` = Pi,
                              `Misspecified value` = Pi_mis
                            ),
                            tibble(
                              Parameter = "pp",
                              `Row index` = rep(1:m, times = m),
                              `Col index` = rep(1:m, each = m),
                              `Well specified value` = c(pp),
                              `Misspecified value` = c(pp_mis)
                            ),
                            tibble(
                              Parameter = "s",
                              `Row index` = 1,
                              `Col index` = 1,
                              `Well specified value` = s,
                              `Misspecified value` = s_mis
                            )
                          ) %>%
                            mutate(
                              run_ID = i,
                              n = n,
                              m = m,
                              K = K,
                              s = s,
                              nu = nu
                            )
                        }

      # Combine results
      all_res <- bind_rows(all_res, result)

      # Message
      cat("\nTotal elapsed time:", difftime(Sys.time(), start, units = "secs"),
          "seconds\n")
      cat("--------------------------------------------------------\n\n")
    }
  }
}

# Close clusters
close(pb)
stopCluster(cl)

# Make sure all is closed
showConnections(all = TRUE) %>%
  as.data.frame %>%
  rownames_to_column('con_id') %>%
  filter(grepl("<-localhost:", x = description)) %>%
  pull(con_id) %>%
  as.integer %>%
  map(getConnection) %>%
  walk(close)

# Save results
write_csv(x = all_res, file = "data/Simulations_Misspec_Par_Dist.csv")
