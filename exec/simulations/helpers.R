# helpers.R — Shared utilities for QATS simulation studies.
#
# Source this file before running any simulation or analysis script:
#   source("helpers.R")

library(QATS)
library(tidyverse)
library(furrr)
library(progressr)
library(arrow)

# --- Design grid construction ------------------------------------------------

#' Build a factorial design grid for simulation studies.
#'
#' For each (m, sigma, n) combination, K values are chosen as
#' {1,2,5} x 10^{0..7} subject to K < n/50.
#' The grid is crossed with rep = 1:n_sim and shuffled for load balancing.
build_grid <- function(m_vals, sigma_vals, n_powers, n_sim) {
  base_K <- c(outer(c(1, 2, 5), 10^(0:7)))

  expand_grid(m = m_vals, sigma = sigma_vals, n_pow = n_powers) |>
    mutate(
      n = 10^n_pow + 1,
      K_list = map(n, \(nn) base_K[base_K < nn / 50])
    ) |>
    unnest(K_list) |>
    rename(K = K_list) |>
    select(-n_pow) |>
    crossing(rep = seq_len(n_sim)) |>
    slice_sample(prop = 1)
}

# --- Timing repetitions ------------------------------------------------------

#' Number of decoder repetitions for timing precision.
#' Longer sequences need fewer reps since each decode is slower.
timing_reps <- function(n) {
  ceiling(1e3 / log10(n)^4)
}

# --- Parameter perturbation ---------------------------------------------------

#' Misspecify HMM parameters by perturbing the transition matrix.
#'
#' Each off-diagonal entry of pp is multiplied by Uniform(1/nu, nu);
#' diagonal entries are left unchanged. Rows are then renormalised.
#' nu = 1 returns the original parameters (no perturbation).
misspecify_par <- function(par_true, m, sigma, nu) {
  perturb <- matrix(runif(m * m, 1 / nu, nu), nrow = m)
  diag(perturb) <- 1
  pp_mis <- par_true$pp * perturb
  pp_mis <- sweep(pp_mis, 1, rowSums(pp_mis), "/")

  set.par(
    yy        = par_true$yy,
    Pi        = par_true$Pi,
    pp        = pp_mis,
    emi.dist  = "normal",
    emi.param = list(mu = 1:m, sigma = rep(sigma, m))
  )
}

# --- Perturbation distribution study ------------------------------------------

#' Generate parameter perturbation data (no decoding, runs in seconds).
#'
#' For each grid point and nu value, constructs the true pp from K/(n-1),
#' applies a random perturbation, and records element-wise values.
#' Used to visualise how much the transition matrix changes with nu.
generate_perturbation_data <- function(m_vals = 2, sigma_vals = 1,
                                       n_powers = 5, n_sim = 100,
                                       nu_vals = c(1, 2, 5, 10, 15, 20)) {
  design <- build_grid(m_vals, sigma_vals, n_powers, n_sim)

  pmap_dfr(design, function(m, sigma, n, K, rep) {
    pp_offdiag <- K / ((n - 1) * (m - 1))
    pp_diag    <- 1 - K / (n - 1)
    pp <- matrix(pp_offdiag, m, m)
    diag(pp) <- pp_diag

    map_dfr(nu_vals, function(nu) {
      perturb <- matrix(runif(m * m, 1 / nu, nu), nrow = m)
      diag(perturb) <- 1
      pp_mis <- pp * perturb
      pp_mis <- sweep(pp_mis, 1, rowSums(pp_mis), "/")

      tibble(
        n = n, m = m, K = K, sigma = sigma, nu = nu, rep = rep,
        row = rep(1:m, times = m),
        col = rep(1:m, each = m),
        pp_true = c(pp),
        pp_mis  = c(pp_mis)
      )
    })
  })
}

# --- Connection cleanup -------------------------------------------------------

#' Close stale localhost connections left by parallel backends.
close_stale_connections <- function() {
  cons <- showConnections(all = TRUE)
  if (nrow(cons) == 0) return(invisible())
  stale <- which(grepl("<-localhost:", cons[, "description"]))
  for (id in as.integer(rownames(cons)[stale])) {
    try(close(getConnection(id)), silent = TRUE)
  }
}
