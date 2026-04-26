# run_misspecified.R — Misspecified HMM simulation study.
#
# Tests robustness of Viterbi and QATS when the transition matrix is perturbed
# by a factor nu. Also generates the parameter perturbation distribution data
# (no decoding, runs in seconds).
#
# Usage:
#   cd exec/simulations
#   Rscript run_misspecified.R
#
# Output:
#   data/misspecified_<timestamp>.parquet     (main results)
#   data/perturbation_<timestamp>.parquet     (parameter perturbation data)

source("helpers.R")

run_misspecified <- function(cores   = 9,
                             n_sim   = 1e3,
                             nu_vals = c(1, 2, 5, 10, 15, 20),
                             d0      = 3,
                             n_seeds = 3,
                             rotate  = FALSE,
                             n_rep   = 5) {
  design <- build_grid(
    m_vals     = 2,
    sigma_vals = 1,
    n_powers   = 5,
    n_sim      = n_sim
  )

  message(nrow(design), " tasks x ", length(nu_vals), " nu values")

  plan(multisession, workers = cores)
  on.exit({
    plan(sequential)
    close_stale_connections()
  })

  handlers(handler_progress(
    format = "[:bar] :percent | :current/:total | ETA :eta",
    clear = FALSE
  ))

  results <- with_progress({
    p <- progressor(steps = nrow(design))

    future_pmap_dfr(
      design,
      function(m, sigma, n, K, rep) {
        par_true <- sample.HMM(
          n = n, m = m, K = K,
          emi.dist  = "normal",
          emi.param = list(mu = 1:m, sigma = rep(sigma, m))
        )

        out <- map_dfr(nu_vals, function(nu) {
          par_mis <- misspecify_par(par_true, m, sigma, nu)

          vit  <- Viterbi.CPP(par_mis, n_rep)
          qats <- QATS.CPP(par_mis, list(
            d0 = d0, n.seeds = n_seeds, rotate = rotate, n.rep = n_rep
          ))

          tibble(
            n = n, m = m, K = K, sigma = sigma, nu = nu, rep = rep,
            method = c("Viterbi", "QATS"),
            time   = c(as.numeric(vit$time), as.numeric(qats$time)),
            l0     = c(lp_norm(par_true$xx, vit$xx, 0),
                       lp_norm(par_true$xx, qats$xx, 0)),
            l2     = c(lp_norm(par_true$xx, vit$xx, 2),
                       lp_norm(par_true$xx, qats$xx, 2))
          )
        })

        p()
        out
      },
      .options = furrr_options(seed = TRUE)
    )
  })

  dir.create("data", showWarnings = FALSE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

  write_parquet(results, sprintf("data/misspecified_%s.parquet", ts))
  message("Saved ", nrow(results), " simulation rows")

  # Parameter perturbation study (fast, no decoding)
  message("Generating perturbation distribution data...")
  perturb <- generate_perturbation_data(nu_vals = nu_vals, n_sim = n_sim)
  write_parquet(perturb, sprintf("data/perturbation_%s.parquet", ts))
  message("Saved ", nrow(perturb), " perturbation rows")

  invisible(results)
}

run_misspecified()
