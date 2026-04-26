# analyze_well_specified.R — Analysis and plots for the well-specified study.
#
# Reads the Parquet output from run_well_specified.R and produces:
#   1. Plot_Time_Ratios.pdf         — Median time ratios across all settings
#   2. Plot_Time_Ratios_Bands.pdf   — Time ratios with 10-90% bands (n = 1e6+1)
#   3. Plot_Errors.pdf              — l0 and l2 errors (n = 1e6+1)
#   4. Plot_Errors_Differences.pdf  — QATS - PMAP error differences
#   5. Plot_Time_Errors.pdf         — Time vs error 2D density
#
# Usage:
#   cd exec/simulations
#   Rscript analyze_well_specified.R

library(tidyverse)
library(patchwork)
library(arrow)
library(ggh4x)

# ── A. Load data ─────────────────────────────────────────────────────────────

parquet_file <- list.files("data", pattern = "^well_specified_.*\\.parquet$",
                           full.names = TRUE) |>
  sort() |>
  tail(1)

message("Loading ", parquet_file)
df_raw <- read_parquet(parquet_file)

# ── B. Pivot and summarise ───────────────────────────────────────────────────

df <- df_raw |>
  mutate(run_ID = row_number()) |>
  pivot_longer(
    cols = time_Viterbi:l2_QATS,
    names_to = c("measurement", "method"),
    names_sep = "_"
  ) |>
  mutate(p = K / (n - 1))

# Time ratios (method vs method)
df_times <- df |>
  filter(measurement == "time") |>
  pivot_wider(names_from = method) |>
  mutate(
    `Viterbi / QATS` = Viterbi / QATS,
    `PMAP / QATS`    = PMAP / QATS
  ) |>
  pivot_longer(c(Viterbi, PMAP, QATS, `Viterbi / QATS`, `PMAP / QATS`),
               names_to = "name") |>
  group_by(n, m, K, sigma, p, name) |>
  summarise(
    q10 = quantile(value, 0.1),
    q50 = quantile(value, 0.5),
    q90 = quantile(value, 0.9),
    .groups = "drop"
  )

# Error metrics
df_errors <- df |>
  filter(measurement != "time") |>
  pivot_wider(names_from = method) |>
  mutate(
    `QATS - Viterbi` = QATS - Viterbi,
    `QATS - PMAP`    = QATS - PMAP,
    `PMAP - Viterbi`  = PMAP - Viterbi
  ) |>
  pivot_longer(c(Viterbi, PMAP, QATS,
                 `QATS - Viterbi`, `QATS - PMAP`, `PMAP - Viterbi`),
               names_to = "name") |>
  mutate(value = value * if_else(measurement == "l2", sqrt(n), 1)) |>
  group_by(n, m, K, sigma, p, measurement, name) |>
  summarise(
    q10 = quantile(value, 0.1),
    q50 = quantile(value, 0.5),
    q90 = quantile(value, 0.9),
    .groups = "drop"
  )

# Display labels
label_n     <- function(n) paste0("n = 1e", log10(n - 1), "+1")
label_sigma <- function(s) if_else(s < 1, paste0("σ = ", s),
                                   paste0("σ = ", s, ".0"))

df_times  <- df_times  |> mutate(n_lab = label_n(n), sigma_lab = label_sigma(sigma))
df_errors <- df_errors |> mutate(n_lab = label_n(n), sigma_lab = label_sigma(sigma))

error_labels <- c(l0 = "Misclassification rate",
                  l1 = "Mean absolute error",
                  l2 = "Root mean squared error")

# ── C. Plot 1: Time ratios ──────────────────────────────────────────────────

tmp <- df_times |>
  filter(name %in% c("Viterbi / QATS", "PMAP / QATS"))

plt <- ggplot(tmp, aes(x = p, y = q50,
                colour = factor(m),
                group = interaction(n, m, name, sigma))) +
  geom_hline(yintercept = 1) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  geom_line(aes(linetype = n_lab), colour = "black", linewidth = 0.4) +
  scale_x_log10() +
  scale_y_log10(
    breaks = c(outer(c(1, 2, 5), 10^(-2:3)))
  ) +
  coord_cartesian(ylim = c(0.5, NA)) +
  facet_grid2(name ~ sigma_lab) +
  labs(x = "p = K / (n - 1)", y = "Time ratio (median)",
       colour = "m", linetype = "n") +
  theme_bw()
ggsave("Plot_Time_Ratios.pdf", plt, width = 6, height = 5, device = cairo_pdf)
message("Wrote Plot_Time_Ratios.pdf")

# ── D. Plot 2: Time ratios with bands (n = 1e6+1) ───────────────────────────

tmp_ratio <- df_times |>
  filter(n == 1e6 + 1, name == "Viterbi / QATS", sigma == 1.0)

p1 <- ggplot(tmp_ratio, aes(x = p, colour = factor(m), fill = factor(m))) +
  geom_hline(yintercept = 1) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.15, colour = NA) +
  geom_line(aes(y = q50), linewidth = 0.4, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10(breaks = c(outer(c(1, 2, 5), 10^(-2:3)))) +
  coord_cartesian(ylim = c(0.5, NA)) +
  labs(x = "p", y = "t(Viterbi) / t(QATS)", colour = "m", fill = "m") +
  theme_bw() +
  theme(legend.position = "none")

tmp_abs <- df_times |>
  filter(n == 1e6 + 1, name %in% c("Viterbi", "QATS", "PMAP"), sigma == 1.0)

p2 <- ggplot(tmp_abs, aes(x = p, colour = factor(m), fill = factor(m),
                           group = interaction(m, name))) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.15, colour = NA) +
  geom_line(aes(y = q50, linetype = name), linewidth = 0.4, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "p", y = "Time [s]", colour = "m", fill = "m", linetype = "Method") +
  theme_bw()

ggsave("Plot_Time_Ratios_Bands.pdf", p1 + p2, width = 7, height = 3, device = cairo_pdf)
message("Wrote Plot_Time_Ratios_Bands.pdf")

# ── E. Plot 3: Errors (l0 and l2, n = 1e6+1) ────────────────────────────────

tmp <- df_errors |>
  filter(
    name %in% c("Viterbi", "QATS", "PMAP"),
    measurement %in% c("l0", "l2"),
    n == 1e6 + 1
  ) |>
  mutate(error_label = error_labels[measurement])

plt <- ggplot(tmp, aes(x = p, group = interaction(m, sigma, measurement, name))) +
  geom_ribbon(aes(ymin = q10, ymax = q90, fill = factor(m)), alpha = 0.15) +
  geom_line(aes(y = q50, colour = factor(m)), linewidth = 0.6, alpha = 0.7) +
  geom_line(aes(y = q50, linetype = name), colour = "black", linewidth = 0.4) +
  scale_x_log10() +
  scale_y_sqrt() +
  facet_grid2(error_label ~ sigma_lab, scales = "free_y", independent = "y") +
  labs(x = "p = K / (n - 1)", y = "Error (median)",
       colour = "m", fill = "m", linetype = "Method") +
  theme_bw()
ggsave("Plot_Errors.pdf", plt, width = 7, height = 5, device = cairo_pdf)
message("Wrote Plot_Errors.pdf")

# ── F. Plot 4: Error differences (QATS - PMAP) ─────────────────────────────

tmp <- df_errors |>
  filter(
    name == "QATS - PMAP",
    measurement %in% c("l0", "l2")
  ) |>
  mutate(error_label = error_labels[measurement])

plt <- ggplot(tmp, aes(x = p, group = interaction(n, m, sigma, error_label))) +
  geom_line(aes(y = q50, colour = factor(m)), linewidth = 0.6, alpha = 0.7) +
  geom_line(aes(y = q50, linetype = n_lab), colour = "black", linewidth = 0.4) +
  scale_x_log10() +
  scale_y_sqrt() +
  facet_grid2(error_label ~ sigma_lab, scales = "free_y", independent = "y") +
  labs(x = "p = K / (n - 1)", y = "Error difference = QATS - PMAP (median)",
       colour = "m", linetype = "n") +
  theme_bw()
ggsave("Plot_Errors_Differences.pdf", plt, width = 7, height = 5, device = cairo_pdf)
message("Wrote Plot_Errors_Differences.pdf")

# ── G. Plot 5: Time vs error 2D density ─────────────────────────────────────

tmp <- df |>
  filter(n == 1e6 + 1, m == 3, sigma == 1.0,
         K / (n - 1) == 1e-4, method == "QATS",
         measurement %in% c("time", "l2")) |>
  pivot_wider(names_from = measurement)

plt <- ggplot(tmp, aes(x = l2, y = time)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette = 4, direction = 1) +
  scale_x_log10(expand = c(0, 0)) +
  scale_y_log10(expand = c(0, 0)) +
  labs(x = "Root mean squared error", y = "Time [s]", fill = "Density") +
  theme_bw()
ggsave("Plot_Time_Errors.pdf", plt, width = 4, height = 3, device = cairo_pdf)
message("Wrote Plot_Time_Errors.pdf")
