# run_well_specified.R — Well-specified HMM simulation study.
#
# Compares Viterbi, PMAP, and QATS across a grid of (m, sigma, n, K) settings.
# Each parallel task runs ONE simulation via the C++ compare_norm() function,
# giving optimal load balancing across workers.
#
# Usage:
#   cd exec/simulations
#   Rscript run_well_specified.R
#
# Output:
#   data/well_specified_<timestamp>.parquet

source("helpers.R")

run_well_specified <- function(cores   = 9,
                               n_sim   = 100,
                               d0      = 3,
                               n_seeds = 3,
                               rotate  = FALSE) {
  design <- build_grid(
    m_vals     = c(2, 3, 5),
    sigma_vals = c(0.1, 1.0),
    n_powers   = 3:6,
    n_sim      = n_sim
  ) |>
    mutate(
      pp_offdiag = K / ((n - 1) * (m - 1)),
      pp_diag    = 1 - K / (n - 1),
      n_rep      = timing_reps(n)
    )

  message(nrow(design), " tasks | ",
          n_distinct(design$m), " m values x ",
          n_distinct(design$sigma), " sigma values x ",
          n_distinct(design$n), " sequence lengths")

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
      function(m, sigma, n, K, rep, pp_offdiag, pp_diag, n_rep) {
        Pi <- rep(1 / m, m)
        pp <- matrix(pp_offdiag, m, m)
        diag(pp) <- pp_diag

        res <- compare_norm(
          n, m, Pi, pp, 1:m, rep(sigma, m),
          d0, n_seeds, rotate, n_rep, n_sim = 1L
        )

        p()

        tibble(
          n = n, m = m, K = K, sigma = sigma, rep = rep,
          time_Viterbi = res$res_Vit[1, 1],
          l0_Viterbi   = res$res_Vit[1, 2],
          l1_Viterbi   = res$res_Vit[1, 3],
          l2_Viterbi   = res$res_Vit[1, 4],
          time_PMAP    = res$res_PMAP[1, 1],
          l0_PMAP      = res$res_PMAP[1, 2],
          l1_PMAP      = res$res_PMAP[1, 3],
          l2_PMAP      = res$res_PMAP[1, 4],
          time_QATS    = res$res_QATS[1, 1],
          l0_QATS      = res$res_QATS[1, 2],
          l1_QATS      = res$res_QATS[1, 3],
          l2_QATS      = res$res_QATS[1, 4]
        )
      },
      .options = furrr_options(seed = TRUE)
    )
  })

  dir.create("data", showWarnings = FALSE)
  file_out <- sprintf("data/well_specified_%s.parquet",
                       format(Sys.time(), "%Y%m%d_%H%M%S"))
  write_parquet(results, file_out)
  message("Saved ", nrow(results), " rows -> ", file_out)

  invisible(results)
}

run_well_specified()
