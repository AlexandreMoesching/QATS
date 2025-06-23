rm(list = ls())
library(QATS)
library(doSNOW)
library(tidyverse)

run_simulation <- function(m, s, n.funcall, n.sim0, cores) {
  # n.funcall = number of function calls for a given setup
  # n.sim0    = number of simulations within a function call
  # n.sim     = n.funcall * n.sim0
  #           = number of simulations for a given setup

  # Fit parameters
  d0 <- 3; n.seeds <- 3; rotate <- FALSE

  # HMM parameters
  Pi <- rep(1, m) / m
  pp <- matrix(NA, nrow = m, ncol = m)
  mu_vec <- 1:m
  s_vec <- rep(s, m)

  # Prepare randomized simulations
  nKpp <- matrix(NA, nrow = 0, ncol = 4,
                 dimnames = list(NULL, c("n", "K", "pp_offdiag", "pp_diag")))
  for (n.pow in 3:6) {
    n <- 10^n.pow + 1
    K.seq <- c(outer(c(1,2,5), 10^(0:7))); K.seq <- K.seq[K.seq < n/5e1]
    nKpp <- rbind(nKpp, cbind(n, K.seq, K.seq / ((n - 1) * (m - 1)), 1 - K.seq / (n - 1)))
  }
  nKpp <- do.call(rbind, replicate(n.funcall, nKpp, simplify = FALSE))
  nKpp.n <- nrow(nKpp)
  nKpp <- nKpp[sample(1:nKpp.n, nKpp.n, replace = FALSE), ]

  # Prepare parallel computing
  cl <- makeSOCKcluster(cores) #, outfile = ""
  registerDoSNOW(cl)
  pb <- txtProgressBar(min = 1, max = nKpp.n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  cl_opts <- list(progress = progress)

  # File name
  tot.time <- Sys.time()
  file.name <-
    paste("m", m,
          "s", floor(s/10^floor(log10(s))), "e", floor(log10(s)), "_",
          format(tot.time, "%Y%m%d_%H%M%S"), ".RData", sep = "")

  # Start simulations
  result <- foreach(i = 1:nKpp.n,
                    .options.snow = cl_opts,
                    .combine = "rbind",
                    .packages = "QATS") %dopar% {
                      # Setup
                      n <- nKpp[i, 1]
                      K <- nKpp[i, 2]

                      # HMM parameters
                      pp <- matrix(nKpp[i, 3], nrow = m, ncol = m)
                      diag(pp) <- nKpp[i, 4]

                      # Number of repetitions for precision of execution time
                      n.rep <- ceiling(1e3/(log10(n))^4)

                      # Run simulations
                      res <- compare_norm(n, m, Pi, pp, mu_vec, s_vec,
                                          d0, n.seeds, rotate,
                                          n.rep, n.sim0)

                      # Return results
                      data.frame(time_Vit  = res$res_Vit[ ,1],
                                 l0_Vit    = res$res_Vit[ ,2],
                                 l1_Vit    = res$res_Vit[ ,3],
                                 l2_Vit    = res$res_Vit[ ,4],
                                 time_PMAP = res$res_PMAP[,1],
                                 l0_PMAP   = res$res_PMAP[,2],
                                 l1_PMAP   = res$res_PMAP[,3],
                                 l2_PMAP   = res$res_PMAP[,4],
                                 time_QATS = res$res_QATS[,1],
                                 l0_QATS   = res$res_QATS[,2],
                                 l1_QATS   = res$res_QATS[,3],
                                 l2_QATS   = res$res_QATS[,4],
                                 n = n,
                                 m = m,
                                 K = K,
                                 s = s,
                                 d0 = d0,
                                 n.seeds = n.seeds,
                                 rotate = rotate)
                    }
  # Save data
  save(file = paste("data/", file.name, sep = ""), list = ls())

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
}

# cores <- parallel::detectCores()
cores <- 7                  # Number of cores used on the server
n.funcall <- 1e1            # Number of function calls for a given setup
n.sim0 <- 1e1               # Number of simulations within a function call
n.sim <- n.funcall * n.sim0 # Number of simulations for a given setup

for (m in c(2, 3, 5)) {
  for (s in c(0.1, 1.0)) {
    # Message
    cat("--------------------------------------------------------\n")
    cat("Setting m = ", m, ", sigma = ", s, "\n", sep = "")
    start <- Sys.time()

    # Run simulations
    run_simulation(m, s, n.funcall, n.sim0, cores)

    # Message
    cat("Total elapsed time:", difftime(Sys.time(), start, units = "secs"),
        "seconds\n")
    cat("--------------------------------------------------------\n\n")
  }
}

# Make sure all is closed
showConnections(all = TRUE) %>%
  as.data.frame %>%
  rownames_to_column('con_id') %>%
  filter(grepl("<-localhost:", x = description)) %>%
  pull(con_id) %>%
  as.integer %>%
  map(getConnection) %>%
  walk(close)
