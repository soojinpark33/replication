# Run with working directory = this folder (source + CSV I/O).

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

set.seed(123)
source("simulation_functions.R")

B <- 500
Ns <- c(500, 1000, 10000)

to_tbl <- function(x) {
  if (is.data.frame(x)) return(as_tibble(x))
  if (is.list(x)) return(as_tibble(x))
  if (is.atomic(x)) return(as_tibble(as.list(x)))
  dplyr::tibble(value = x)
}

poc_zeta_from_vec <- function(v) {
  tibble(
    zeta_static = as.numeric(v[["zeta_static_poc"]]),
    zeta_stoch = as.numeric(v[["zeta_stoch_poc"]])
  )
}

run_method_reps <- function(pop, n_draw, method_id) {
  map_dfr(seq_len(B), function(rep_i) {
    df_s <- sample_n(pop, n_draw)
    out <- if (method_id == "poc") {
      poc_zeta_from_vec(poc(df_s))
    } else if (method_id == "imp") {
      imp(df_s)
    } else if (method_id == "lasso") {
      lasso(df_s, lambda_rule = "min")
    } else if (method_id == "lasso_1se") {
      lasso(df_s, lambda_rule = "1se")
    } else {
      stop("Unknown method_id", call. = FALSE)
    }
    if (method_id != "poc") {
      out <- out |>
        to_tbl() |>
        rename(zeta_static = zeta_static_imp, zeta_stoch = zeta_stoch_imp)
    }
    transmute(out, n = n_draw, method = method_id, rep = rep_i, zeta_static, zeta_stoch)
  })
}

bootstrap_zeta <- function(pop) {
  methods <- c("poc", "imp", "lasso", "lasso_1se")
  map_dfr(Ns, function(n_draw) {
    bind_rows(lapply(methods, function(m) run_method_reps(pop, n_draw, m)))
  })
}

write_zeta_truth_csv <- function(path, zeta_static_true, zeta_stoch_true) {
  write.csv(
    data.frame(zeta_static_true = zeta_static_true, zeta_stoch_true = zeta_stoch_true),
    path,
    row.names = FALSE
  )
}

# Truth: E_z[ mu(1,1,z,w*) - mu(1,0,z,w*) - (mu(0,1,z,w*) - mu(0,0,z,w*)) ] with structural mean; stochastic uses sample means on R=0.
truth_zeta <- function(popdat, w_fix, mu_struct) {
  zeta_static <- 0
  for (z in sort(unique(popdat$Z))) {
    pz <- mean(popdat$Z == z)
    d_r1 <- mu_struct(1, 1, z, w_fix) - mu_struct(1, 0, z, w_fix)
    d_r0 <- mu_struct(0, 1, z, w_fix) - mu_struct(0, 0, z, w_fix)
    zeta_static <- zeta_static + pz * (d_r1 - d_r0)
  }
  zeta_stoch <- 0
  for (z in sort(unique(popdat$Z))) {
    pz <- mean(popdat$Z == z)
    idx0 <- popdat$R == 0 & popdat$Z == z
    d_r0 <- mean(popdat$Y[idx0 & popdat$X == 1]) - mean(popdat$Y[idx0 & popdat$X == 0])
    Wp0 <- popdat$Wprime[idx0]
    wp_tilde <- if (length(Wp0)) mean(Wp0) else NA_real_
    d_r1_G <- if (length(Wp0)) {
      mu_struct(1, 1, z, wp_tilde) - mu_struct(1, 0, z, wp_tilde)
    } else {
      NA_real_
    }
    zeta_stoch <- zeta_stoch + pz * (d_r1_G - d_r0)
  }
  c(static = zeta_static, stoch = zeta_stoch)
}

# --- Correct specification ---
N <- 100000L
Z <- rbinom(N, 1, 0.5)
p_R <- plogis(-0.5 + 0.8 * Z)
R <- rbinom(N, 1, p_R)
W <- rnorm(N, mean = 0.2 + 0.5 * Z + 0.7 * R, sd = 1)
X <- rbinom(N, 1, 0.5)
Wprime <- rnorm(N, mean = -0.3 + 0.4 * Z + 0.6 * R + 0.3 * W, sd = 1)
mu_Y <- 0.7 + 1 * R + 0.5 * X + 0.3 * Z + 0.4 * W + 0.2 * Wprime + 0.1 * (X * Wprime) +
  0.5 * (R * X) + 0.2 * (R * X * Wprime) + 0.5 * (X * Z) - 0.7 * (R * Z) + 0.5 * (R * X * Z)
Y <- rnorm(N, mean = mu_Y, sd = 1)
popdat1 <- data.frame(Z, R, W, X, Wprime, Y)

w_fix <- 1
mu_Y_struct <- function(R, X, Z, wp) {
  mu_W <- 0.2 + 0.5 * Z + 0.7 * R
  0.7 + R + 0.5 * X + 0.3 * Z + 0.4 * mu_W + 0.2 * wp + 0.1 * (X * wp) +
    0.5 * (R * X) + 0.2 * (R * X * wp) + 0.5 * (X * Z) - 0.7 * (R * Z) + 0.5 * (R * X * Z)
}
tz1 <- truth_zeta(popdat1, w_fix, mu_Y_struct)
zeta_static_true_popdat1 <- tz1["static"]
zeta_stoch_true_popdat1 <- tz1["stoch"]
write_zeta_truth_csv("zeta_truth_popdat1.csv", zeta_static_true_popdat1, zeta_stoch_true_popdat1)

set.seed(123)
results_popdat1 <- bootstrap_zeta(popdat1)
write.csv(results_popdat1, "zeta_boot_popdat1_all.csv", row.names = FALSE)

# --- Misspecified outcome ---
Z <- rbinom(N, 1, 0.5)
p_R <- plogis(-0.5 + 0.8 * Z)
R <- rbinom(N, 1, p_R)
W <- rnorm(N, mean = 0.2 + 0.5 * Z + 0.7 * R, sd = 1)
X <- rbinom(N, 1, 0.5)
Wprime <- rnorm(N, mean = -0.3 + 0.4 * Z + 0.6 * R + 0.3 * W, sd = 1)
mu_Y1 <- 0.7 + R + 0.5 * X + 0.3 * Z + 0.4 * W + 0.2 * Wprime + 0.1 * (X * Wprime) +
  0.5 * (R * X) + 0.2 * (R * X * Wprime) + 0.5 * (X * Z) - 0.7 * (R * Z) + 0.5 * (R * X * Z) -
  0.5 * (X * W) - 0.5 * (X * W * Z)
Y <- rnorm(N, mean = mu_Y1, sd = 1)
popdat2 <- data.frame(Z, R, W, X, Wprime, Y)

mu_Y_struct_inc <- function(R, X, Z, wp) {
  mu_W <- 0.2 + 0.5 * Z + 0.7 * R
  0.7 + R + 0.5 * X + 0.3 * Z + 0.4 * mu_W + 0.2 * wp + 0.1 * (X * wp) +
    0.5 * (R * X) + 0.2 * (R * X * wp) + 0.5 * (X * Z) - 0.7 * (R * Z) + 0.5 * (R * X * Z) -
    0.5 * (X * mu_W) - 0.5 * (X * mu_W * Z)
}
tz2 <- truth_zeta(popdat2, w_fix, mu_Y_struct_inc)
zeta_static_true_popdat2 <- tz2["static"]
zeta_stoch_true_popdat2 <- tz2["stoch"]
write_zeta_truth_csv("zeta_truth_popdat2.csv", zeta_static_true_popdat2, zeta_stoch_true_popdat2)

set.seed(123)
results_popdat2 <- bootstrap_zeta(popdat2)
write.csv(results_popdat2, "zeta_boot_popdat2_all.csv", row.names = FALSE)

# --- Plots (can re-run from CSV if objects absent) ---
load_truth_or_stop <- function(csv, label) {
  if (!file.exists(csv)) {
    stop("Missing ", csv, ". Run the DGP block for ", label, " first.", call. = FALSE)
  }
  tr <- read.csv(csv)
  list(static = tr$zeta_static_true[[1]], stoch = tr$zeta_stoch_true[[1]])
}

if (!exists("zeta_static_true_popdat1", inherits = TRUE) ||
    !exists("zeta_stoch_true_popdat1", inherits = TRUE)) {
  t1 <- load_truth_or_stop("zeta_truth_popdat1.csv", "popdat1")
  zeta_static_true_popdat1 <- t1$static
  zeta_stoch_true_popdat1 <- t1$stoch
}
if (!exists("zeta_static_true_popdat2", inherits = TRUE) ||
    !exists("zeta_stoch_true_popdat2", inherits = TRUE)) {
  t2 <- load_truth_or_stop("zeta_truth_popdat2.csv", "popdat2")
  zeta_static_true_popdat2 <- t2$static
  zeta_stoch_true_popdat2 <- t2$stoch
}
if (!exists("results_popdat1", inherits = TRUE)) {
  results_popdat1 <- read.csv("zeta_boot_popdat1_all.csv")
}
if (!exists("results_popdat2", inherits = TRUE)) {
  results_popdat2 <- read.csv("zeta_boot_popdat2_all.csv")
}

center_results <- function(results, zeta_static_true, zeta_stoch_true) {
  results |>
    pivot_longer(cols = c(zeta_static, zeta_stoch), names_to = "parameter", values_to = "est") |>
    mutate(
      truth = case_when(
        parameter == "zeta_static" ~ zeta_static_true,
        parameter == "zeta_stoch" ~ zeta_stoch_true,
        TRUE ~ NA_real_
      ),
      diff = est - truth,
      parameter = factor(
        parameter,
        levels = c("zeta_static", "zeta_stoch"),
        labels = c("zeta[static]", "zeta[stoch]")
      )
    )
}

res1_centered <- center_results(results_popdat1, zeta_static_true_popdat1, zeta_stoch_true_popdat1)
res2_centered <- center_results(results_popdat2, zeta_static_true_popdat2, zeta_stoch_true_popdat2)

methods_plot <- c("poc", "imp", "lasso")
method_labels <- c("POC", "Imp", "Lasso")
res1_centered <- res1_centered |>
  filter(method %in% methods_plot) |>
  mutate(method = factor(method, levels = methods_plot, labels = method_labels)) |>
  droplevels()
res2_centered <- res2_centered |>
  filter(method %in% methods_plot) |>
  mutate(method = factor(method, levels = methods_plot, labels = method_labels)) |>
  droplevels()

yrng <- range(c(res1_centered$diff, res2_centered$diff), na.rm = TRUE, finite = TRUE)
pad <- 0.05 * diff(yrng)
if (!is.finite(pad) || pad == 0) pad <- 0.1
ylims <- c(yrng[1] - pad, yrng[2] + pad)

plot_per_n <- function(dat, n_val, title, ylims, n_title_size = 10) {
  ggplot(filter(dat, n == n_val), aes(x = method, y = diff)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "red", linewidth = 0.7) +
    geom_boxplot(outlier.alpha = 0.25, width = 0.7) +
    facet_grid(. ~ parameter, labeller = label_parsed) +
    labs(title = title, x = "Method", y = "Estimate - truth") +
    coord_cartesian(ylim = ylims) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = n_title_size),
      panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey85"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text.x = element_text(size = 11),
      axis.text.x = element_text(size = 10)
    )
}

p1_500 <- plot_per_n(res1_centered, 500L, "n = 500", ylims, n_title_size = 13)
p1_1000 <- plot_per_n(res1_centered, 1000L, "n = 1,000", ylims, n_title_size = 13)
p1_10000 <- plot_per_n(res1_centered, 10000L, "n = 10,000", ylims, n_title_size = 13)
p2_500 <- plot_per_n(res2_centered, 500L, "n = 500", ylims, n_title_size = 13)
p2_1000 <- plot_per_n(res2_centered, 1000L, "n = 1,000", ylims, n_title_size = 13)
p2_10000 <- plot_per_n(res2_centered, 10000L, "n = 10,000", ylims, n_title_size = 13)

row_title <- function(text, title_size = 16) {
  ggplot() +
    theme_void() +
    labs(title = text) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = title_size))
}

final_fig <-
  row_title("Correct specification", title_size = 18) /
  (p1_500 | p1_1000 | p1_10000) /
  row_title("Incorrect specification", title_size = 18) /
  (p2_500 | p2_1000 | p2_10000) +
  plot_layout(heights = c(0.07, 1, 0.07, 1))

print(final_fig)
ggsave("simulation_replication_suyeon_fig.png", final_fig, width = 14, height = 9, dpi = 120)
