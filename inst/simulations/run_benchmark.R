
library(scoutFDR)

## --------------------------- ##
##  Load calibration objects   ##
## --------------------------- ##

cal_df  <- readRDS("inst/calibration/cal_df.rds")
cqr_fit <- readRDS("inst/calibration/cqr_fit.rds")

feature_names <- cqr_fit$feature_names
stopifnot(all(feature_names %in% colnames(cal_df)))

## ------------------ ##
## Simulation design  ##
## ------------------ ##

alpha_grid  <- c(0.01, 0.05, 0.10)
pi0_grid <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
power_grid <- c(0.25, 0.35, 0.50, 0.65)
m_grid <- c(2000, 5000, 10000)
n_grid <- c(3, 5, 10)

n_rep <- 100 # number of replicates per scenario

set.seed(17)

sim_grid <- expand.grid(
  pi0 = pi0_grid,
  target_power = power_grid,
  m = m_grid,
  n = n_grid,
  stringsAsFactors = FALSE
)

effect_grid <- scoutFDR:::build_effect_grid(
  power_grid = power_grid,
  n_grid     = n_grid,
  alpha      = 0.05,   # reference alpha for power calibration
  sd_ratio   = 0.3
)

## ---------------------- ##
## One-scenario function  ##
## ---------------------- ##

run_one_scenario <- function(s_id) {
  row <- sim_grid[s_id, ]

  pi0_s <- row$pi0
  pow_s <- row$target_power
  m_s   <- row$m
  n_s   <- row$n

  message(
    "Scenario ", s_id, "/", nrow(sim_grid),
    ": pi0=", pi0_s,
    ", power=", pow_s,
    ", m=", m_s,
    ", n=", n_s
  )

  # Scenario-specific seed
  set.seed(1000 + s_id)

  res_list <- vector(
    mode = "list",
    length = n_rep * length(alpha_grid)
  )
  idx <- 1L

  for (rep in seq_len(n_rep)) {
    # 1) Simulate one dataset for this scenario
    dat <- scoutFDR:::simulate_dataset(
      pi0         = pi0_s,
      target_power = pow_s,
      m           = m_s,
      n           = n_s,
      effect_grid = effect_grid
    )

    # 2) Evaluate all FDR levels on the same dataset
    for (alpha in alpha_grid) {
      df <- scoutFDR:::run_all_methods(
        pvals        = dat$pvals,
        tstats       = dat$tstats,
        n_samples    = n_s,
        cal_df       = cal_df,
        cqr_fit      = cqr_fit,
        feature_names = feature_names,
        is_null      = dat$is_null,
        alpha        = alpha,
        true_pi0     = pi0_s
      )

      df$scenario_id  <- s_id
      df$rep          <- rep
      df$pi0_true     <- pi0_s
      df$power_target <- pow_s
      df$m            <- m_s
      df$n            <- n_s
      df$alpha        <- alpha

      res_list[[idx]] <- df
      idx <- idx + 1L
    }
  }

  do.call(rbind, res_list)
}

## -------------------------- ##
## Interate on all scenarios  ##
## -------------------------- ##

scenario_ids <- seq_len(nrow(sim_grid))
start_time <- Sys.time()
res_list <- lapply(scenario_ids, run_one_scenario)
end_time <- Sys.time()
message("Total simulation time: ", as.numeric(end_time - start_time), " ", attr(end_time - start_time, "units"))

## ---------------------- ##
## Save combined results  ##
## ---------------------- ##

# Bind all results into one big data.frames
res_benchmark <- do.call(rbind, res_list)
saveRDS(res_benchmark, "inst/simulations/res_benchmark.rds")

message("Saved combined results to inst/simulations/res_benchmark.rds")

#Analyse the results of the benchmark simulation
res_benchmark <- readRDS("inst/simulations/res_benchmark.rds")
library(dplyr)
summary_df <- res_benchmark %>%
  group_by(method, pi0_true, power_target, m, n, alpha) %>%
  summarize(
    avg_fdr = mean(FDR),
    avg_power = mean(power),
    avg_discoveries = mean(discoveries),
    avg_pi0_hat = mean(pi0_hat),
    .groups = 'drop'
  )
saveRDS(summary_df, "inst/simulations/res_benchmark_summary.rds")
message("Saved summary results to inst/simulations/res_benchmark_summary.rds")

# Make plots summarizing the benchmark results
library(ggplot2)
library(tidyr)
summary_df <- readRDS("inst/simulations/res_benchmark_summary.rds")
# Example plot: boxplot of average FDR by method,true pi0 and alpha level
fdr_plot <- ggplot(summary_df, aes(x = method, y = avg_fdr, fill = method)) +
  geom_boxplot() +
  facet_grid(pi0_true ~ alpha, labeller = labeller(
    pi0_true = function(x) paste("True π0 =", x),
    alpha = function(x) paste("FDR target α =", x)
  ), scales = "free_y") +
  geom_hline(data = unique(summary_df[c("alpha")]), 
           aes(yintercept = alpha),
           linetype = "dashed", color = "red") +
  theme_bw() +
  theme(legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Average FDR across simulation scenarios",
       x = "Method",
       y = "Average FDR")
ggsave("inst/simulations/fdr_plot.png", plot = fdr_plot, width = 10, height = 12)
message("Saved average FDR plot to inst/simulations/fdr_plot.png")

