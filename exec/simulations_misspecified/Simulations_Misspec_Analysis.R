rm(list = ls()); try(dev.off(), silent = TRUE)
library(tidyverse)
library(magrittr)
library(patchwork)

#__________________________________________________________________________ ----
# A. Load and prepare QATS and Viterbi simulations ----
df <- arrow::read_parquet("data/Simulations_Misspec_HMM_20250623_132547.parquet")

# Get important information
(df %>% pull(     d0) %>% sort() %>% unique() -> d0     )
(df %>% pull(n.seeds) %>% sort() %>% unique() -> n.seeds)
(df %>% pull( rotate) %>% sort() %>% unique() -> rotate )

# Remove unnecessary columns and pivot
df %<>%
  select(-(d0:rotate)) %>%
  pivot_longer(time:l2, names_to = "meas")

# Rescale
df %<>%
  mutate(value = value * ifelse(meas == "l2", sqrt(n), 1))

# Compute difference with baseline (nu = 1)
df %<>%
  rename(`Raw value` = value) %>%
  group_by(n, m, K, s, method, meas, rep) %>%
  mutate(
    `Baseline value` = `Raw value`[nu == 1],
    `Difference with baseline` = `Raw value` - `Baseline value`,
    `% difference from well-specified baseline` = 100 * ifelse(`Baseline value` > 0, `Difference with baseline` / `Baseline value`, 0)
  ) %>%
  select(-`Difference with baseline`, -`Baseline value`) %>%
  ungroup()

# Pivot longer
df %<>%
  pivot_longer(
    cols = c(`Raw value`, `% difference from well-specified baseline`),
    names_to = "value_type"
  )

# Summarise
df_sum_1 <-
  df %>%
  group_by(across(-c(rep, value))) %>%
  summarise(
    q_25 = quantile(value, p = 0.25),
    q_50 = quantile(value, p = 0.50),
    q_75 = quantile(value, p = 0.75),
    .groups = "drop"
  )

#__________________________________________________________________________ ----
# B. Load and prepare parameter simulations ----
df <- read_csv("data/Simulations_Misspec_Par_Dist.csv", show_col_types = FALSE)

# Filter
df %<>% filter(Parameter == "pp", `Row index` == 1)

# Rename, compute differences and pivot
df %<>%
  rename(
    rep = run_ID,
    `Baseline value` = `Well specified value`,
    `Raw value` = `Misspecified value`
  ) %>%
  mutate(
    `Difference with baseline` = `Raw value` - `Baseline value`,
    `% difference from well-specified baseline` = 100 * ifelse(`Baseline value` > 0, `Difference with baseline` / `Baseline value`, 0)
  ) %>%
  select(-`Difference with baseline`, -`Baseline value`) %>%
  pivot_longer(
    cols = c(`Raw value`, `% difference from well-specified baseline`),
    names_to = "value_type"
  )

# Compute quantiles
df_sum_2 <-
  df %>%
  group_by(across(-c(rep, value))) %>%
  summarise(
    q_25 = quantile(value, p = 0.25),
    q_50 = quantile(value, p = 0.50),
    q_75 = quantile(value, p = 0.75),
    .groups = "drop"
  )

#__________________________________________________________________________ ----
# C. Combine tables ----
df_sum_1 %<>%
  mutate(
    name =
      case_when(
        meas == "time" ~ "Time",
        meas == "l0" ~ "Misclassification rate",
        meas == "l2" ~ "Root mean squared error",
      ),
    .after = nu
  ) %>%
  select(-meas)

df_sum_2 %<>%
  mutate(
    name = paste0("p_", `Row index`, `Col index`),
    method = "Both",
    .after = nu
  ) %>%
  select(-c(Parameter:`Col index`))

# Combine tables
df_sum <-
  bind_rows(df_sum_1, df_sum_2)

# Get sequences
n_seq  <- df_sum %>% pull(n ) %>% unique() %>% sort() %>% {log10(. - 1)}
m_seq  <- df_sum %>% pull(m ) %>% unique() %>% sort()
s_seq  <- df_sum %>% pull(s ) %>% unique() %>% sort()
nu_seq <- df_sum %>% pull(nu) %>% unique() %>% sort()
df_sum %>% pull(name) %>% unique() %>% sort()

# Name parameters nicely
df_sum %<>%
  mutate(
    p = K / (n - 1),
    .after = nu
  ) %>%
  mutate(
    n_disp = paste0("n = 1e", log10(n - 1)),
    n_disp = factor(n_disp, paste0("n = 1e", n_seq)),
    .after = n
  ) %>%
  mutate(
    m_disp = paste("m =", m),
    m_disp = factor(m_disp, paste("m =", m_seq)),
    .after = m
  ) %>%
  mutate(
    s_disp = paste0("\u03C3 = ", s),
    s_disp = factor(s_disp, paste0("\u03C3 = ", s_seq)),
    .after = s
  ) %>%
  mutate(
    nu_disp = paste("\u03BD =", nu),
    nu_disp = factor(nu_disp, levels = paste("\u03BD =", nu_seq)),
    .after = nu
  ) %>%
  mutate(
    name = recode(name,
                  "p_11" = "tilde(p)[11]",          # math
                  "p_12" = "tilde(p)[12]",          # math
                  "Misclassification rate"   = "'Misclassification rate'",      # quoted
                  "Root mean squared error" = "'Root mean squared error'",      # quoted
                  "Time"                    = "'Time'"                          # quoted
    ),
    name = factor(
      name,
      levels = c("tilde(p)[11]", "tilde(p)[12]",
                 "'Misclassification rate'",
                 "'Root mean squared error'",
                 "'Time'")
    )
  )

#__________________________________________________________________________ ----
# D. Plot ----

math_labeller <- function(lbl) {
  recode(lbl,
         "p_11" = "tilde(p)[11]",
         "p_12" = "tilde(p)[12]",
         .default = lbl           # everything else as-is, no spaces issue
  )
}

cairo_pdf(filename = "Plot_Error_Time_Misspec.pdf", width = 8.5, height = 8.5)
# Plot all but for sigma = 1, since this is the only interesting setting
df_sum %>%
  filter(
    nu > 1,
    s == 1,
    m == 2,
    n == 1e5 + 1,
    value_type == "% difference from well-specified baseline",
    nu <= 20,
    name != "'Root mean squared error'"
    ) %>%
  ggplot(aes(x = p)) +
  geom_ribbon(aes(ymin = q_25, ymax = q_75, fill = nu_disp, col = nu_disp, lty = method), alpha = 0.2, lwd = 0.2) +
  geom_line(aes(y = q_50, col = nu_disp, lty = method), lwd = 0.5) +
  scale_x_log10() +
  scale_linetype_manual(values = c(
    "Both" = 1,
    "QATS" = 1,
    "Viterbi" = 2
  )) +
  # scale_y_sqrt() +
  ggh4x::facet_nested(
    name ~ nu_disp,
    scales = "free",
    # independent = "y",
    labeller = ggplot2::labeller(            # custom labeller
      name = ggplot2::label_parsed,          # parse rows
      .default = ggplot2::label_value        # columns untouched
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(
    y = "Median relative difference from well-specified baseline (\u03BD = 1) [%]", lty = "Method",
    col = "Perturbation\nparameter", fill = "Perturbation\nparameter"
  ) # + ggtitle(label = "Results for \u03C3 = 1, m = 2, and n = 1e5 + 1")
dev.off()

df_sum %>%
  filter(
    s == 1,
    m == 2,
    n == 1e5 + 1,
    value_type == "Raw value",
    nu <= 20,
    name != "'Root mean squared error'"
  ) %>%
  ggplot(aes(x = p)) +
  geom_ribbon(aes(ymin = q_25, ymax = q_75, fill = nu_disp, col = nu_disp, lty = method), alpha = 0.2, lwd = 0.2) +
  geom_line(aes(y = q_50, col = nu_disp, lty = method), lwd = 0.5) +
  #scale_x_log10() +
  scale_linetype_manual(values = c(
    "Both" = 1,
    "QATS" = 1,
    "Viterbi" = 2
  )) +
  # scale_y_sqrt() +
  ggh4x::facet_nested(
    name ~ nu_disp,
    scales = "free",
    # independent = "y",
    labeller = ggplot2::labeller(            # custom labeller
      name = ggplot2::label_parsed,          # parse rows
      .default = ggplot2::label_value        # columns untouched
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(
    y = "Median values", lty = "Method",
    col = "Perturbation\nparameter", fill = "Perturbation\nparameter"
  ) # + ggtitle(label = "Results for \u03C3 = 1, m = 2, and n = 1e5 + 1")
