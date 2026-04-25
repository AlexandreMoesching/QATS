# analyze_misspecified.R — Analysis and plots for the misspecified study.
#
# Reads the Parquet outputs from run_misspecified.R and produces:
#   Plot_Error_Time_Misspec.pdf
#
# Usage:
#   cd exec/simulations
#   Rscript analyze_misspecified.R

library(tidyverse)
library(patchwork)
library(arrow)
library(ggh4x)

# ── A. Load simulation results ───────────────────────────────────────────────

sim_file <- list.files("data", pattern = "^misspecified_.*\\.parquet$",
                       full.names = TRUE) |> sort() |> tail(1)
pert_file <- list.files("data", pattern = "^perturbation_.*\\.parquet$",
                        full.names = TRUE) |> sort() |> tail(1)

message("Loading ", sim_file)
message("Loading ", pert_file)

df_sim  <- read_parquet(sim_file)
df_pert <- read_parquet(pert_file)

# ── B. Prepare simulation data ──────────────────────────────────────────────

df_sim <- df_sim |>
  pivot_longer(c(time, l0, l2), names_to = "meas") |>
  mutate(value = value * if_else(meas == "l2", sqrt(n), 1))

# Compute % difference from well-specified baseline (nu = 1)
df_sim <- df_sim |>
  group_by(n, m, K, sigma, method, meas, rep) |>
  mutate(
    baseline = value[nu == 1],
    pct_diff = 100 * if_else(baseline > 0, (value - baseline) / baseline, 0)
  ) |>
  ungroup()

# Reshape: raw value and % difference as separate rows
df_sim_long <- df_sim |>
  pivot_longer(
    cols      = c(value, pct_diff),
    names_to  = "value_type",
    values_to = "val"
  ) |>
  mutate(
    value_type = recode(value_type,
      value    = "Raw value",
      pct_diff = "% difference from baseline"
    )
  )

# Summarise across replications
df_sim_sum <- df_sim_long |>
  group_by(n, m, K, sigma, nu, method, meas, value_type) |>
  summarise(
    q25 = quantile(val, 0.25),
    q50 = quantile(val, 0.50),
    q75 = quantile(val, 0.75),
    .groups = "drop"
  ) |>
  mutate(
    p = K / (n - 1),
    name = recode(meas,
      time = "'Time'",
      l0   = "'Misclassification rate'",
      l2   = "'Root mean squared error'"
    ),
    name = factor(name, levels = c("'Misclassification rate'",
                                   "'Root mean squared error'",
                                   "'Time'"))
  )

# ── C. Prepare perturbation data ────────────────────────────────────────────

df_pert_sum <- df_pert |>
  filter(row == 1) |>
  mutate(
    pct_diff = 100 * if_else(pp_true > 0, (pp_mis - pp_true) / pp_true, 0)
  ) |>
  pivot_longer(
    cols      = c(pp_mis, pct_diff),
    names_to  = "value_type",
    values_to = "val"
  ) |>
  mutate(
    value_type = recode(value_type,
      pp_mis   = "Raw value",
      pct_diff = "% difference from baseline"
    )
  ) |>
  group_by(n, m, K, sigma, nu, row, col, value_type) |>
  summarise(
    q25 = quantile(val, 0.25),
    q50 = quantile(val, 0.50),
    q75 = quantile(val, 0.75),
    .groups = "drop"
  ) |>
  mutate(
    p = K / (n - 1),
    method = "Both",
    name = paste0("tilde(p)[", row, col, "]"),
    name = factor(name, levels = c("tilde(p)[11]", "tilde(p)[12]"))
  )

# ── D. Combine ──────────────────────────────────────────────────────────────

# Align column sets
df_plot <- bind_rows(
  df_sim_sum |> select(n, m, K, sigma, nu, p, method, name, value_type,
                       q25, q50, q75),
  df_pert_sum |> select(n, m, K, sigma, nu, p, method, name, value_type,
                        q25, q50, q75)
) |>
  mutate(
    nu_lab = factor(paste("\u03BD =", nu),
                    levels = paste("\u03BD =", sort(unique(nu))))
  )

# ── E. Plot ─────────────────────────────────────────────────────────────────

cairo_pdf("Plot_Error_Time_Misspec.pdf", width = 8.5, height = 8.5)
df_plot |>
  filter(
    nu > 1,
    sigma == 1,
    m == 2,
    n == 1e5 + 1,
    value_type == "% difference from baseline",
    name != "'Root mean squared error'"
  ) |>
  ggplot(aes(x = p)) +
  geom_ribbon(aes(ymin = q25, ymax = q75,
                  fill = nu_lab, colour = nu_lab, linetype = method),
              alpha = 0.2, linewidth = 0.2) +
  geom_line(aes(y = q50, colour = nu_lab, linetype = method),
            linewidth = 0.5) +
  scale_x_log10() +
  scale_linetype_manual(values = c(Both = 1, QATS = 1, Viterbi = 2)) +
  facet_nested(
    name ~ nu_lab,
    scales = "free",
    labeller = labeller(name = label_parsed, .default = label_value)
  ) +
  labs(
    x = "p = K / (n - 1)",
    y = "Median relative difference from baseline (\u03BD = 1) [%]",
    linetype = "Method",
    colour = "Perturbation\nparameter",
    fill = "Perturbation\nparameter"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
message("Wrote Plot_Error_Time_Misspec.pdf")

# ── F. Supplementary: raw values ────────────────────────────────────────────

cairo_pdf("Plot_Raw_Values_Misspec.pdf", width = 8.5, height = 8.5)
df_plot |>
  filter(
    sigma == 1,
    m == 2,
    n == 1e5 + 1,
    value_type == "Raw value",
    name != "'Root mean squared error'"
  ) |>
  ggplot(aes(x = p)) +
  geom_ribbon(aes(ymin = q25, ymax = q75,
                  fill = nu_lab, colour = nu_lab, linetype = method),
              alpha = 0.2, linewidth = 0.2) +
  geom_line(aes(y = q50, colour = nu_lab, linetype = method),
            linewidth = 0.5) +
  scale_x_log10() +
  scale_linetype_manual(values = c(Both = 1, QATS = 1, Viterbi = 2)) +
  facet_nested(
    name ~ nu_lab,
    scales = "free",
    labeller = labeller(name = label_parsed, .default = label_value)
  ) +
  labs(
    x = "p = K / (n - 1)",
    y = "Median values",
    linetype = "Method",
    colour = "Perturbation\nparameter",
    fill = "Perturbation\nparameter"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
message("Wrote Plot_Raw_Values_Misspec.pdf")
