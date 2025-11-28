# ============================================================
#      SCOUT-FDR: Diagnostics for the Calibrated CQR Model
# ============================================================

library(scoutFDR)
library(ggplot2)
library(dplyr)
library(purrr)
library(quantregForest)

# Calibration data
cal_df <- readRDS("inst/calibration/cal_df.rds")

# Trained CQR model
cqr_fit <- readRDS("inst/calibration/cqr_fit.rds")

# Features used for training (from the model)
feature_names <- cqr_fit$feature_names

# Conformal miscoverage level (tau = 1 - alpha_conf)
alpha_conf <- 1 - cqr_fit$tau   # since tau = 0.7 → alpha_conf = 0.3

# True π0 from simulation
pi0_true <- cal_df$pi0

pi0_pred <- apply(
  cal_df[, feature_names, drop = FALSE],
  1,
  function(x) scoutFDR:::predict_qtau(cqr_fit, setNames(as.numeric(x), feature_names))
)

# Create X_target_all (all feature rows)
X_target_all <- cal_df[, feature_names, drop = FALSE]

# Vectorized conformal upper bounds
source("inst/calibration/diagnostics/compute_upper_bounds.R")

prep <- prepare_conformal_cqr(
    cqr_fit       = cqr_fit,
    cal_df        = cal_df,
    feature_names = feature_names,
    X_target_all  = cal_df[, feature_names],
    alpha         = alpha_conf
)

pi0_upper <- pbapply::pbsapply(
  seq_len(nrow(cal_df)),
  function(i) predict_conformal_cqr(prep, cal_df[i, feature_names])
)

# COVERAGE DIAGNOSTICS

coverage <- mean(pi0_true <= pi0_upper)

message("Overall empirical coverage: ", round(coverage, 3),
        " (target: ", 1 - alpha_conf, ")")

df_cov_ind <- data.frame(
  covered = factor(pi0_true <= pi0_upper, levels = c(FALSE, TRUE))
)

p1 <- ggplot(df_cov_ind, aes(x = covered)) +
  geom_bar(fill = "steelblue") +
  geom_hline(yintercept = nrow(df_cov_ind) * coverage,
             linetype = "dashed", color = "red") +
  annotate("text",
           x = 1.5,
           y = nrow(df_cov_ind) * coverage,
           label = paste0("Coverage = ", round(coverage, 3)),
           vjust = -1,
           color = "red",
           size = 4) +
  labs(title = "Empirical Coverage Indicator",
       subtitle = paste0("Overall coverage: ", round(coverage, 3),
                         " (target ", 1 - alpha_conf, ")"),
       x = "pi0_true ≤ pi0_upper?",
       y = "Count") +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 10))

ggsave("inst/calibration/diagnostics/coverage_barplot.png", p1, width = 6, height = 6)


# Coverage by pi0 bins
df_cov <- data.frame(pi0_true, pi0_upper) %>%
  mutate(bin = cut(
    pi0_true,
    breaks = seq(0, 1, by = 0.1),
    include.lowest = TRUE
  ))

df_cov_plot <- df_cov %>%
  group_by(bin) %>%
  summarise(cov = mean(pi0_true <= pi0_upper), .groups = "drop")

# Add global coverage for reference
overall_cov <- mean(pi0_true <= pi0_upper)

p2 <- ggplot(df_cov_plot, aes(x = bin, y = cov)) +
  geom_col(fill = "darkorange") +
  geom_hline(yintercept = 1 - alpha_conf,
             linetype = 2, color = "red") +
  annotate("text",
           x = Inf, y = 1 - alpha_conf,
           label = paste0("Target = ", 1 - alpha_conf),
           hjust = 1.1, vjust = -0.5,
           color = "red", size = 3.5) +
  ylim(0, 1) +
  labs(
    title = "Coverage by true π₀ bins",
    subtitle = paste0(
      "Overall empirical coverage = ",
      round(overall_cov, 3),
      " (target ", 1 - alpha_conf, ")"
    ),
    x = "π₀ bin",
    y = "Empirical coverage"
  ) +
  theme_minimal()

ggsave("inst/calibration/diagnostics/coverage_by_pi0.png",
       p2, width = 6, height = 4)


# POINT PREDICTION ACCURACY (RAW CQR QUANTILE)

rmse <- sqrt(mean((pi0_pred - pi0_true)^2))
bias <- mean(pi0_pred - pi0_true)

message("RMSE(π₀_pred) = ", round(rmse, 4))
message("Bias(π₀_pred) = ", round(bias, 4))

p4 <- ggplot(data.frame(pi0_true, pi0_pred),
             aes(x = pi0_true, y = pi0_pred)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(title = "Point prediction of π₀ (raw CQR quantile)",
       x = "True π₀", y = "Predicted τ-quantile (CQR)") +
  theme_minimal()

ggsave("inst/calibration/diagnostics/point_prediction.png",
       p4, width = 6, height = 4)


