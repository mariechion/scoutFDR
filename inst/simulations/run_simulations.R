## run_simulations.R
## Run full simulation study for scoutFDR paper

#library(scoutFDR)
#library(future.apply)

## ------------------------------------------------------------------
## 0. Load calibration objects
## ------------------------------------------------------------------

cal_df  <- readRDS("inst/calibration/cal_df.rds")
cqr_fit <- readRDS("inst/calibration/cqr_fit.rds")

feature_names <- cqr_fit$feature_names

## Optional: small sanity check
stopifnot(all(feature_names %in% colnames(cal_df)))

## ------------------------------------------------------------------
## 1. Simulation design
## ------------------------------------------------------------------

alpha_grid  <- c(0.01, 0.05, 0.10)
pi0_grid    <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
power_grid  <- c(0.25, 0.35, 0.50, 0.65)
m_grid      <- c(2000, 5000, 10000)
n_grid      <- c(3, 5, 10)

n_rep <- 10  # adjust if you want a lighter run first

set.seed(17)

sim_grid <- expand.grid(
  pi0          = pi0_grid,
  target_power = power_grid,
  m            = m_grid,
  n            = n_grid,
  stringsAsFactors = FALSE
)

## ------------------------------------------------------------------
## 2. Effect size lookup (same for all scenarios)
## ------------------------------------------------------------------

effect_grid <- scoutFDR:::build_effect_grid(
  power_grid = power_grid,
  n_grid     = n_grid,
  alpha      = 0.05,   # reference alpha for power calibration
  sd_ratio   = 0.3
)

## ------------------------------------------------------------------
## 3. One-scenario runner
## ------------------------------------------------------------------

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

  # Scenario-specific seed for reproducibility (still compatible with future.seed=TRUE)
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

## ------------------------------------------------------------------
## 4. Parallel execution with future.apply
## ------------------------------------------------------------------

scenario_ids <- seq_len(nrow(sim_grid))
start_time <- Sys.time()
res_list <- lapply(scenario_ids, run_one_scenario)
end_time <- Sys.time()
message("Total simulation time: ", as.numeric(end_time - start_time), " ", attr(end_time - start_time, "units"))

## ------------------------------------------------------------------
## 5. Save combined results
## ------------------------------------------------------------------

# Bind all results into one big data.frame
sim_results <- do.call(rbind, res_list)
saveRDS(sim_results, "inst/simulations/sim_results.rds")

message("Saved combined results to inst/simulations/sim_results.rds")
