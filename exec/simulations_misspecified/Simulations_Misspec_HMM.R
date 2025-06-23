# ------------------------------------------------------------
#  QATS misspecification simulations
# ------------------------------------------------------------

# -------------------------- libraries -----------------------
library(tidyverse)  # dplyr, purrr, tibble…
library(QATS)
library(furrr)      # parallel purrr
library(progressr)  # nice progress bars
library(arrow)      # fast columnar I/O

set.seed(123)       # reproducible across cores

# Enable ETA
progressr::handlers("cli")

# ----------------------- design grid ------------------------
build_design <- function(
    m_vals      = 2,
    sigma_vals  = 1,
    n_powers    = 5,
    n_sim       = 1e2
) {
  # K grid depends on n, so build that first
  base_K <- c(outer(c(1, 2, 5), 10^(0:7)))          # 1,2,5 × 10^0…10^7
  crossing(n_pow = n_powers, m = m_vals) %>%
    mutate(
      n      = 10^n_pow + 1,
      K_grid = map(n, \(nn) base_K[base_K < nn / 50])
    ) %>%
    unnest(K_grid) %>%
    rename(K = K_grid) %>%
    select(-n_pow) %>%
    crossing(rep = seq_len(n_sim), sigma = sigma_vals) %>%
    arrange(sample(row_number()))
}

# ------------------ miss-specify parameters -----------------
misspecify_par <- function(par_true, m, sigma, nu) {
  # Pi_mis <- par_true$Pi * runif(m, 1 / nu, nu)
  # Pi_mis <- Pi_mis / sum(Pi_mis)
  Pi_mis <- par_true$Pi

  perturb <- matrix(runif(m * m, 1 / nu, nu), nrow = m)
  diag(perturb) <- 1
  pp_mis <- par_true$pp * perturb
  pp_mis <- sweep(pp_mis, 1, rowSums(pp_mis), '/')
  # pp_mis <- par_true$pp

  mu_mis <- 1:m
  # sigma_mis <- sigma * runif(1, 4 / (3 + nu), (3 + nu) / 4)
  sigma_mis <- sigma
  sigma_mis <- rep(sigma_mis, m)

  set.par(
    yy        = par_true$yy,
    Pi        = Pi_mis,
    pp        = pp_mis,
    emi.dist  = 'normal',
    emi.param = list(mu = mu_mis, sigma = sigma_mis)
  )
}

# -------------------- run ONE simulation --------------------
simulate_once <- function(m,
                          sigma,
                          nu_vals,
                          n,
                          K,
                          rep,
                          d0 = 3,
                          n.seeds = 3,
                          rotate = FALSE) {
  # n.rep <- ceiling(1e3 / (log10(n))^4)
  n.rep <- 5

  # -------- generate the sequence ONCE ------------------
  par_true <- sample.HMM(
    n        = n,
    m        = m,
    K        = K,
    emi.dist = "normal",
    emi.param = list(mu = 1:m, sigma = rep(sigma, m))
  )

  # -------- run methods for EACH nu ---------------------
  map_dfr(nu_vals, function(nu) {

    par_mis <- misspecify_par(par_true, m, sigma, nu)
    opts_q  <- list(d0 = d0, n.seeds = n.seeds, rotate = rotate, n.rep = n.rep)

    tibble(method = c("Viterbi", "QATS"),
           fun    = list(Viterbi.CPP, QATS.CPP),
           par    = list(par_mis, par_mis),
           opts   = list(n.rep, opts_q)) %>%
      pmap_dfr(function(method, fun, par, opts) {
        seg <- fun(par, opts)
        tibble(
          method = method,
          time   = as.numeric(seg$time),
          l0     = lp_norm(par_true$xx, seg$xx, 0),
          l2     = lp_norm(par_true$xx, seg$xx, 2)
        )
      }) %>%
      mutate(
        n = n,
        m = m,
        K = K,
        s = sigma,
        nu = nu,
        d0 = d0,
        n.seeds = n.seeds,
        rotate = rotate,
        rep = rep,
        .before = 1
      )
  })
}

# ------------------- run ALL simulations --------------------
#' Run the full grid of miss-specification simulations in parallel
#'
#' @description
#' Generates a factorial design of hidden-Markov-model settings (see
#' [`build_design()`]) and runs `simulate_once()` for every row in parallel
#' using **future**/**furrr**.
#' Each call produces segmentation quality metrics (`l0`, `l2`) and run-time
#' for four decoding methods (Viterbi/QATS × true/miss-specified parameters).
#' The combined result is saved as an Apache **Parquet** file (fast, compact,
#' column-oriented) and returned invisibly.
#'
#' @details
#' The work distribution is handled by `future_pmap_dfr()` (**furrr**) on the
#' backend chosen by `future::plan()`.  A progress bar is displayed through
#' **progressr**; the bar advances once per design row, _after_ the
#' corresponding simulation finishes, so it reflects true completion (not just
#' task submission).  Reproducibility across machines is guaranteed by
#' `furrr_options(seed = TRUE)`, which seeds each worker deterministically.
#'
#' When the computations end, the function switches back to sequential
#' execution (`plan(sequential)`) so that later code in the session does not
#' keep the cluster alive by accident.  Results are persisted with
#' `arrow::write_parquet()`.  The file name encodes the current timestamp to
#' avoid accidental overwrites.
#'
#' @param cores  `integer(1)` – number of parallel workers.  Defaults to
#'   `future::availableCores()`.  **Tip:** on a shared server choose a value
#'   below the physical core count so other users are not starved.
#' @param out_dir  `character(1)` – directory where the Parquet file is written.
#'   Created on-the-fly if it does not exist.
#' @param chunk  `integer(1)` – row chunk size used by Arrow when writing the
#'   file.  Tuning this can improve write speed for extremely large outputs;
#'   the default (`5e5`) is a good compromise for 10^6–10^7 rows.
#'
#' @return
#' (Invisibly) a `tibble` with one row **per simulation × decoding method** and
#' the columns
#' \describe{
#'   \item{n, m, K, s, nu, d0, n.seedsm, rotate, rep}{design variables}
#'   \item{method}{decoding method label}
#'   \item{time}{execution time in seconds}
#'   \item{l0, l2}{segmentation error metrics}
#' }
#' Side-effect: a Parquet file `misspec_<timestamp>.parquet` in `out_dir`.
#'
#' @examples
#' \dontrun{
#' ## use 7 workers, save to default directory
#' future::plan(multisession, workers = 7)
#' res <- run_all_simulations()
#'
#' ## read the saved file later:
#' library(arrow)
#' big_tbl <- arrow::read_parquet("Data_Misspec_Refactor/misspec_20250527_153045.parquet")
#' }
#' @export
run_all_simulations <- function(
    cores   = future::availableCores(), # how many parallel R processes to use
    chunk   = 5e5,                      # row-chunk size for Arrow (I/O tuning)
    nu_vals = c(1, 2, 5, 10, 15, 20)
) {

  # ------------------------------------------------------------------
  # 1. EXPERIMENT DESIGN: one row = one simulation -------------------
  # ------------------------------------------------------------------
  design <- build_design()   # tibble with columns m, sigma, nu, n, K, rep …
  total  <- nrow(design)     # we’ll need this for the progress bar

  # ------------------------------------------------------------------
  # 2. PARALLEL BACKEND (future) -------------------------------------
  # ------------------------------------------------------------------
  # 'multisession' = launch <cores> *separate* R sessions on the same machine.
  # Every time we later call a 'future_*' function (provided by {furrr}),
  # the work will be shipped to one of these sessions.
  future::plan(multisession, workers = cores)

  # ------------------------------------------------------------------
  # 3. ENABLE PROGRESS BARS (progressr) ------------------------------
  # ------------------------------------------------------------------
  # Without this call, progress bars stay silent.  Setting 'global = TRUE'
  # means every 'with_progress({...})' in this session will display a bar.
  progressr::handlers(global = TRUE)

  # ------------------------------------------------------------------
  # 4. RUN THE SIMULATIONS -------------------------------------------
  # ------------------------------------------------------------------
  # 'with_progress()' ACTIVATES the bar for everything inside its braces.
  results <- progressr::with_progress({

    # a) construct a "progressor" with 'total' steps
    p <- progressr::progressor(steps = total)

    # b) 'future_pmap_dfr()' = parallel version of purrr::pmap_dfr()
    #    - maps an anonymous function over the *rows* of 'design'
    #    - runs each call in a parallel worker
    #    - row-binds ('_dfr') all returned tibbles into one big tibble
    furrr::future_pmap_dfr(
      design,
      ~ {
        #   Inside each worker:
        #   -------------------
        #   • simulate_once(...) does the heavy computation
        #   • p() tells the progress bar we've finished *one* design row
        #   • the resulting tibble is returned to future_pmap_dfr()
        out <- simulate_once(..., nu_vals = nu_vals)
        p()            # advance bar by 1 step
        out            # returned to furrr for row-binding
      },
      .options = furrr::furrr_options(seed = TRUE)  # reproducible RNG
    )
  })

  # ------------------------------------------------------------------
  # 5. BACK TO SINGLE PROCESS ----------------------------------------
  # ------------------------------------------------------------------
  # Prevents accidental parallelism in code that runs *after* this function.
  future::plan(sequential)

  # ------------------------------------------------------------------
  # 6. SAVE RESULTS AS PARQUET (arrow) -------------------------------
  # ------------------------------------------------------------------
  # Why Parquet?
  #   • column-oriented, highly compressed  →  file is small
  #   • random access  →  can read just a few columns later
  #   • language-agnostic  →  Python, R, Julia all read it natively
  #
  # The file name includes a timestamp to avoid overwriting previous runs.
  file_out <- file.path(
    sprintf("data/Simulations_Misspec_HMM_%s.parquet",
            format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  arrow::write_parquet(results, file_out, chunk_size = chunk)

  # Friendly confirmation for the console:
  message("✓ Saved ", nrow(results), " rows → ", file_out)

  # ------------------------------------------------------------------
  # 8. RETURN VALUE ---------------------------------------------------
  # ------------------------------------------------------------------
  # 'invisible()' prevents the large tibble from printing to the console
  # unless the caller explicitly captures it, e.g. res <- run_all_simulations()
  invisible(results)
}

# --------------------------- usage ---------------------------
plan(multisession, workers = 7)
res <- run_all_simulations()
showConnections(all = TRUE) %>%
  as.data.frame %>%
  rownames_to_column('con_id') %>%
  filter(grepl("<-localhost:", x = description)) %>%
  pull(con_id) %>%
  as.integer %>%
  map(getConnection) %>%
  walk(close)
