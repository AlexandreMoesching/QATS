rm(list = ls()); try(dev.off(), silent = TRUE)
library(tidyverse)
library(magrittr)
library(patchwork)

#__________________________________________________________________________ ----
# A. Combine data ----
all.file.names <- list.files(path = "data")
all.file.names.n <- length(all.file.names)
df <- NULL
for (i in 1:all.file.names.n) {
  cat(
    "Loading file ",
    ifelse(i < 10, "0", ""), i,
    "/", all.file.names.n, "\n", sep = "")
  load(paste0("data/", all.file.names[i]))
  df %<>% bind_rows(as_tibble(result))
}

#__________________________________________________________________________ ----
# B. Create summary tables ----
# Get important information
(df %>% pull(     d0) %>% sort() %>% unique() -> d0     )
(df %>% pull(n.seeds) %>% sort() %>% unique() -> n.seeds)
(df %>% pull( rotate) %>% sort() %>% unique() -> rotate )

# Remove unnecessary columns
df %<>% select(-(d0:rotate))

# Add row number
df %<>% mutate(run_ID = row_number(), .before = everything())

# Pivot longer
df %<>%
  pivot_longer(
    cols = time_Vit:l2_QATS,
    names_to = c("measurement", "method"),
    names_sep = "_"
  )

# Rename Vit to Viterbi
df %<>% mutate(method = ifelse(method == "Vit", "Viterbi", method))

# Time ratios
df_times <-
  df %>%
  filter(measurement == "time") %>%
  pivot_wider(names_from = method) %>%
  mutate(
    Viterbi_QATS = Viterbi / QATS,
    Viterbi_PMAP = Viterbi / PMAP,
    PMAP_QATS = PMAP / QATS
  ) %>%
  pivot_longer(Viterbi:PMAP_QATS) %>%
  group_by(n, m, K, s, name) %>%
  summarise(
    q_1 = quantile(value, 0.1),
    q_5 = quantile(value, 0.5),
    q_9 = quantile(value, 0.9)
  ) %>%
  ungroup() %>%
  pivot_longer(cols = q_1:q_9, names_to = "quantile")

# Error differences
df_errors <-
  df %>%
  filter(measurement != "time") %>%
  pivot_wider(names_from = method) %>%
  mutate(
    QATS_Viterbi = QATS - Viterbi,
    QATS_PMAP = QATS - PMAP,
    PMAP_Viterbi = PMAP - Viterbi
  ) %>%
  pivot_longer(Viterbi:PMAP_Viterbi) %>%
  group_by(n, m, K, s, measurement, name) %>%
  summarise(
    q_1 = quantile(value, 0.1),
    q_5 = quantile(value, 0.5),
    q_9 = quantile(value, 0.9)
  ) %>%
  ungroup() %>%
  pivot_longer(cols = q_1:q_9, names_to = "quantile")

# Rescale errors
df_errors %<>%
  mutate(
    value = value * ifelse(measurement == "l2", sqrt(n), 1)
  )

# Compute p
df_times  %<>% mutate(p = K / (n - 1), .after = K)
df_errors %<>% mutate(p = K / (n - 1), .after = K)

# Rename errors
df_errors %<>%
  mutate(
    error_type_display =
      case_when(
        measurement == "l0" ~ "Misclassification rate",
        measurement == "l1" ~ "Mean absolute error",
        measurement == "l2" ~ "Root mean squared error",
      ),
    error_type_display =
      paste0(
        error_type_display,
        ifelse(grepl("_", name), " difference", "")
      ),
    .after = measurement
  )

#__________________________________________________________________________ ----
# C. Plot times ----

# Make sure n is a factor ----
df_times %<>% mutate(n = paste0("1e", log10(n - 1), "+1"))

# Make sure m is a factor ----
df_times %<>% mutate(m = factor(m, levels = c(2, 3, 5, 10)))

# Modify names ----
df_times %<>%
  mutate(
    name_displayed = gsub("_", " / ", name),
    name_displayed = gsub("PMAP", "t(PMAP)", name_displayed),
    name_displayed = gsub("QATS", "t(QATS)", name_displayed),
    name_displayed = gsub("Viterbi", "t(Viterbi)", name_displayed),
    .after = name
  )

# Modify standard deviations ----
df_times %<>%
  mutate(
    s_displayed =
      case_when(
        s < 1 ~ paste0("\u03C3 = ", s),
        TRUE ~ paste0("\u03C3 = ", s, ".0"),
      )
  )

# Plot time ratios ----
tmp <-
  df_times %>%
  filter(
    name_displayed %in% c("t(Viterbi) / t(QATS)", "t(PMAP) / t(QATS)"),
    s %in% c(0.1, 1.0),
    quantile == "q_5"
  ) %>%
  rename(`Time ratio` = value)

cairo_pdf("Plot_Time_Ratios.pdf", width = 6, height = 5)
ggplot() +
  theme_bw() +
  geom_hline(yintercept = 1) +
  geom_line(
    mapping = aes(x = p, y = `Time ratio`, col = m,
                  group = interaction(n, m, name_displayed, s_displayed)),
    data = tmp,
    lwd = 2,
    alpha = 0.7
  ) +
  geom_line(
    mapping = aes(x = p, y = `Time ratio`, lty = n,
                  group = interaction(n, m, name_displayed, s_displayed)),
    data = tmp,
    col = "black"
  ) +
  scale_x_log10() +
  scale_y_log10(breaks = c(outer(c(1, 2, 5), (-2:3), FUN = function(x, y) x*10^y))) +
  coord_cartesian(ylim = c(0.5, NA)) +
  ggh4x::facet_grid2(name_displayed ~ s_displayed)
dev.off()

# Plot time ratios with bands for n = 1e6+1 ----
tmp <-
  df_times %>%
  filter(
    n == "1e6+1",
    name_displayed %in% c("t(Viterbi) / t(QATS)"),
    s %in% c(1.0)
  ) %>%
  pivot_wider(names_from = quantile) %>%
  rename(`Time ratio` = q_5)

p1 <-
  ggplot() +
  theme_bw() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1) +
  geom_line(
    mapping = aes(x = p, y = `Time ratio`, col = m),
    data = tmp,
    lwd = 2,
    alpha = 0.7
  ) +
  geom_line(
    mapping = aes(x = p, y = `Time ratio`, group = m),
    data = tmp,
    col = "black"
  ) +
  geom_ribbon(
    mapping = aes(x = p, ymin = q_1, ymax = q_9, fill = m),
    alpha = 0.2,
    data = tmp
  ) +
  ylab("Time ratio = t(Viterbi) / t(QATS)") +
  scale_x_log10() +
  scale_y_log10(breaks = c(outer(c(1, 2, 5), (-2:3), FUN = function(x, y) x*10^y))) +
  coord_cartesian(ylim = c(0.5, NA))

# Plot times with bands for n = 1e6+1 ----
tmp <-
  df_times %>%
  filter(
    n == "1e6+1",
    name %in% c("Viterbi", "QATS", "PMAP"),
    s %in% c(1.0)
  ) %>%
  pivot_wider(names_from = quantile) %>%
  rename(Time = q_5, Method = name)

p2 <-
  ggplot() +
  theme_bw() +
  geom_line(
    mapping = aes(x = p, y = Time, col = m, group = interaction(m, Method)),
    data = tmp,
    lwd = 2,
    alpha = 0.7
  ) +
  geom_line(
    mapping = aes(x = p, y = Time, lty = Method, group = interaction(m, Method)),
    data = tmp,
    col = "black"
  ) +
  geom_ribbon(
    mapping = aes(x = p, ymin = q_1, ymax = q_9, fill = m, group = interaction(m, Method)),
    alpha = 0.2,
    data = tmp
  ) +
  ylab("Time [s]") +
  scale_x_log10() +
  scale_y_log10() + # breaks = c(outer(c(1, 2, 5), (-5:3), FUN = function(x, y) x*10^y))
  coord_cartesian()

cairo_pdf("Plot_Time_Ratios_Bands.pdf", width = 7, height = 3)
p1 + p2
dev.off()

#__________________________________________________________________________ ----
# D. Plot errors ----

# Make sure n is a factor
df_errors %<>% mutate(n = paste0("1e", log10(n - 1), "+1"))

# Make sure m is a factor
df_errors %<>%
  mutate(m_displayed = paste0("m = ", m),
         m_displayed =
           factor(m_displayed,
                  levels = sapply(c(2, 3, 5, 10), function(x) paste0("m = ", x))),
         m = factor(m, levels = c(2, 3, 5, 10)))

# Modify names
df_errors %<>%
  mutate(
    name_displayed = gsub("_", " - ", name),
    .after = name
  )

# Modify standard deviations
df_errors %<>%
  mutate(
    s_displayed =
      case_when(
        s < 1 ~ paste0("\u03C3 = ", s),
        TRUE ~ paste0("\u03C3 = ", s, ".0"),
      )
  )

# Plot l0 errors ----
tmp <-
  df_errors %>%
  filter(
    name_displayed %in% c("Viterbi", "QATS", "PMAP"),
    s %in% c(0.1, 1.0),
    measurement %in% c("l0"),
    n == "1e6+1"
  ) %>%
  pivot_wider(names_from = quantile) %>%
  rename(`Misclassification rate` = q_5, Method = name)

ggplot() +
  theme_bw() +
  geom_line(
    mapping = aes(x = p, y = `Misclassification rate`,
                  lty = Method),
    data = tmp
  ) +
  geom_ribbon(
    mapping = aes(x = p, ymin = q_1, ymax = q_9, group = Method),
    alpha = 0.2,
    data = tmp
  ) +
  scale_x_log10() +
  ggh4x::facet_grid2(s_displayed ~ m_displayed, scales = "free_y")

# Plot l2 errors ----
tmp <-
  df_errors %>%
  filter(
    name_displayed %in% c("Viterbi", "QATS", "PMAP"),
    s %in% c(0.1, 1.0),
    measurement %in% c("l2"),
    n == "1e6+1"
  ) %>%
  pivot_wider(names_from = quantile) %>%
  rename(`Root mean squared error` = q_5, Method = name)

ggplot() +
  theme_bw() +
  geom_line(
    mapping = aes(x = p, y = `Root mean squared error`,
                  lty = Method),
    data = tmp
  ) +
  geom_ribbon(
    mapping = aes(x = p, ymin = q_1, ymax = q_9, group = Method),
    alpha = 0.2,
    data = tmp
  ) +
  scale_x_log10() +
  ggh4x::facet_grid2(s_displayed ~ m_displayed, scales = "free_y")

# Plot l0 and l2 errors ----
tmp <-
  df_errors %>%
  filter(
    name_displayed %in% c("Viterbi", "QATS", "PMAP"),
    s %in% c(0.1, 1.0),
    measurement %in% c("l0", "l2"),
    n == "1e6+1"
  ) %>%
  pivot_wider(names_from = quantile) %>%
  rename(Error = q_5, Method = name)

cairo_pdf("Plot_Errors.pdf", width = 7, height = 5)
ggplot() +
  theme_bw() +
  geom_line(
    mapping = aes(x = p, y = Error, col = m, group = interaction(m, s, measurement, Method)),
    data = tmp,
    alpha = 0.7,
    lwd = 2
  ) +
  geom_line(
    mapping = aes(x = p, y = Error, lty = Method, group = interaction(m, s, measurement, Method)),
    data = tmp,
    col = "black"
  ) +
  geom_ribbon(
    mapping = aes(x = p, ymin = q_1, ymax = q_9, fill = m, group = interaction(m, s, measurement, Method)),
    alpha = 0.2,
    data = tmp
  ) +
  scale_x_log10() +
  scale_y_sqrt() +
  ggh4x::facet_grid2(error_type_display ~ s_displayed, scales = "free_y", independent = "y")
dev.off()

# Plot error differences ----
tmp <-
  df_errors %>%
  filter(
    name_displayed %in% c("QATS - PMAP"),
    s %in% c(0.1, 1.0),
    measurement %in% c("l0", "l2"),
    quantile == "q_5"
  ) %>%
  rename(`Error difference` = value)

cairo_pdf("Plot_Errors_Differences.pdf", width = 7, height = 5)
ggplot() +
  theme_bw() +
  geom_line(
    mapping = aes(x = p, y = `Error difference`,
                  col = m, group = interaction(n, m, s, error_type_display)),
    data = tmp,
    lwd = 2,
    alpha = 0.7
  ) +
  geom_line(
    mapping = aes(x = p, y = `Error difference`,
                  lty = n, group = interaction(n, m, s, error_type_display)),
    data = tmp,
    col = "black"
  ) +
  ylab("Error difference = QATS - PMAP") +
  scale_x_log10() +
  scale_y_sqrt() +
  ggh4x::facet_grid2(error_type_display ~ s_displayed, scales = "free_y", independent = "y")
dev.off()

#__________________________________________________________________________ ----
# E. Error vs Time ----

cairo_pdf("Plot_Time_Errors.pdf", width = 4, height = 3)
df %>%
  mutate(p = K / (n - 1), .after = K) %>%
  filter(n == 1e6 + 1, m == 3, s == 1.0, p == 1e-4, method == "QATS", measurement %in% c("time", "l2")) %>%
  pivot_wider(names_from = measurement) %>%
  ggplot(aes(x = l2, y = time)) +
  theme_bw() +
  # theme(legend.position = "none") +
  scale_x_log10(expand = c(0, 0)) +
  scale_y_log10(expand = c(0, 0)) +
  xlab("Root mean squared error") +
  ylab("Time [s]") +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette = 4, direction = 1) +
  labs(fill = "Density")
dev.off()
