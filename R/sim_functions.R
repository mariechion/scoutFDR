#' @importFrom stats qt pnorm dnorm integrate rnorm t.test uniroot


############################################################
## 1. Effect size calibration helpers
############################################################

# Compute power of a two-sided t-test for a given true effect delta
# under equal-variance, balanced two-sample design (n per group),
# using a normal approximation to the non-central t distribution.
power_for_delta <- function(delta, n, alpha = 0.05) {
  nu <- 2 * (n - 1)
  # two-sided test critical value
  tcrit <- qt(1 - alpha / 2, df = nu)
  # standard normal approximation
  z1 <- tcrit - delta * sqrt(n / 2)
  z2 <- -tcrit - delta * sqrt(n / 2)
  1 - pnorm(z1) + pnorm(z2)
}

# Average power when delta ~ N(mu, (sd_ratio * mu)^2)
avg_power_for_mu <- function(mu, n, alpha = 0.05, sd_ratio = 0.3) {
  sd <- sd_ratio * mu
  f <- function(delta) {
    dnorm(delta, mean = mu, sd = sd) * power_for_delta(delta, n, alpha = alpha)
  }
  integrate(f, lower = -Inf, upper = Inf)$value
}

# Solve for mu such that the average power matches target_power
# for a given n and sd_ratio.
calibrate_mu <- function(target_power, n, alpha = 0.05, sd_ratio = 0.3) {
  f_root <- function(mu) avg_power_for_mu(mu, n, alpha = alpha, sd_ratio = sd_ratio) - target_power
   # Rough bracket: near 0 effect => low power, larger effect => high power
  uniroot(f_root, interval = c(0.01, 5), tol = 1e-3)$root
}

# Build a lookup table of (target_power, n) -> (mu_delta, sigma_delta)
# This lets us reuse the same grid for all simulations.
build_effect_grid <- function(power_grid, n_grid, alpha = 0.05, sd_ratio = 0.3) {
  eg <- expand.grid(
    target_power = power_grid,
    n = n_grid,
    stringsAsFactors = FALSE
  )
  
  eg$mu_delta <- mapply(
    calibrate_mu,
    target_power = eg$target_power,
    n = eg$n,
    MoreArgs = list(alpha = alpha, sd_ratio = sd_ratio)
  )
  eg$sigma_delta <- sd_ratio * eg$mu_delta
  
  eg
}

############################################################
## 2. Data generating mechanism
############################################################

#' Simulate one dataset for a given (pi0, target_power, m, n) scenario.
#'
#' @param pi0 True proportion of nulls
#' @param target_power Target *average* power used for calibration
#' @param m Total number of hypotheses
#' @param n Sample size per group
#' @param effect_grid Output of build_effect_grid()
#'
#' @returns A list with:
#'   - pvals   - numeric vector of p-values (length m)
#'   - tstats  - numeric vector of t-statistics (length m)
#'   - is_null - logical vector of length m, TRUE if hypothesis is null
#'
#' @export
simulate_dataset <- function(pi0, target_power, m, n, effect_grid) {
  # look up mu_delta, sigma_delta for this (target_power, n)
  idx <- which(effect_grid$target_power == target_power &
               effect_grid$n == n)
  if (length(idx) != 1L) {
    stop("Could not uniquely match (target_power, n) in effect_grid.")
  }
  mu_delta <- effect_grid$mu_delta[idx]
  sd_delta <- effect_grid$sigma_delta[idx]

  # Null / non-null labels
  m0 <- round(pi0 * m)
  m1 <- m - m0
  is_null <- c(rep(TRUE, m0), rep(FALSE, m1))
  is_null <- sample(is_null)  # shuffle

  # Effect sizes
  deltas <- numeric(m)
  if (m1 > 0L) {
    deltas[!is_null] <- rnorm(m1, mean = mu_delta, sd = sd_delta)
  }

  # Generate data matrices: rows = hypotheses, cols = samples
  X1 <- matrix(rnorm(m * n, mean = 0, sd = 1), nrow = m, ncol = n)
  X2 <- matrix(rnorm(m * n, mean = 0, sd = 1), nrow = m, ncol = n) +
    matrix(deltas, nrow = m, ncol = n)

  # Row-wise means
  mean1 <- rowMeans(X1)
  mean2 <- rowMeans(X2)

  # Unbiased sample variances (df = n - 1)
  s1_sq <- rowSums((X1 - mean1)^2) / (n - 1)
  s2_sq <- rowSums((X2 - mean2)^2) / (n - 1)

  # Pooled variance and t-stats (equal variances)
  df <- 2 * n - 2
  sp_sq <- ((n - 1) * s1_sq + (n - 1) * s2_sq) / df
  se <- sqrt(sp_sq * (2 / n))
  tstats <- (mean1 - mean2) / se

  # Two-sided p-values from t-distribution
  pvals <- 2 * stats::pt(abs(tstats), df = df, lower.tail = FALSE)

  list(
    pvals   = pvals,
    tstats  = tstats,
    is_null = is_null
  )
}


############################################################
## 3. Multiple testing methods
############################################################

#' Benjamini–Hochberg procedure (BH).
#'
#' Internal helper to compute BH rejections at a given FDR level.
#'
#' @param pvals Numeric vector of p-values.
#' @param level Target FDR level in (0, 1].
#'
#' @return A logical vector of length \code{length(pvals)} indicating
#'   which hypotheses are rejected.
#'
#' @keywords internal
bh_reject <- function(pvals, level) {
  m <- length(pvals)
  o <- order(pvals)
  p_sorted <- pvals[o]
  thresh <- (seq_len(m) * level) / m
  which_ok <- which(p_sorted <= thresh)
  k <- if (length(which_ok) == 0L) 0L else max(which_ok)
  reject <- logical(m)
  if (k > 0L) {
    reject[o[seq_len(k)]] <- TRUE
  }
  reject
}


#' Run all FDR methods on a single vector of p-values.
#'
#' 
#' @param pvals numeric vector of p-values
#' @param tstats numeric vector of t-statistics (used for Cheng)
#' @param n_samples number of samples per group (used for Cheng)
#' @param cal_df calibration dataset (used for scoutFDR)
#' @param cqr_fit cqr model (used for scoutFDR)
#' @param feature_names feature names (used for scoutFDR) 
#' @param is_null logical vector indicating which hypotheses are null
#' @param alpha FDR level
#' @param true_pi0 true proportion of nulls (used for Oracle method)
#' 
#' @return data.frame with one row per method containing:
#'     method, pi0_hat, FDR, power, discoveries
#
run_all_methods <- function(pvals, tstats, n_samples,
  cal_df, cqr_fit, feature_names, is_null, alpha, true_pi0) {
  m  <- length(pvals)
  m0 <- sum(is_null)
  m1 <- m - m0
  
  methods_out <- list()
  
  ## Benjamini-Hochberg method
  rej_bh <- bh_reject(pvals, level = alpha)
  methods_out[["BH"]] <- list(
    pi0_hat    = 1,
    rejections = rej_bh
  )
  
  ## Storey (lambda fixed at 0.5 as a simple choice)
  lambda <- 0.5
  pi0_storey <- mean(pvals > lambda) / (1 - lambda)
  pi0_storey <- min(max(pi0_storey, 1e-6), 1)  # keep in (0, 1]
  rej_storey <- bh_reject(pvals, level = alpha / pi0_storey)
  methods_out[["Storey"]] <- list(
    pi0_hat    = pi0_storey,
    rejections = rej_storey
  )

  ## Pounds-Cheng: average Storey over lambda
  lambda_grid <- seq(0.5, 0.9, by = 0.05)
  pi0_pc_vec <- sapply(lambda_grid, function(lam) {
    est <- mean(pvals > lam) / (1 - lam)
    est
  })
  pi0_pc <- mean(pi0_pc_vec)
  pi0_pc <- min(max(pi0_pc, 1e-6), 1)
  rej_pc <- bh_reject(pvals, level = alpha / pi0_pc)
  methods_out[["PoundsCheng"]] <- list(
    pi0_hat    = pi0_pc,
    rejections = rej_pc
  )

  ## Simple Beta: fit Beta(β,1) on p <= 0.25, correct Storey at λ = 0.5
  p_low <- pvals[pvals <= 0.25]
  if (length(p_low) >= 50) {
    beta_hat <- -mean(log(p_low))
    eta_lam  <- (1 - lambda)^beta_hat   
    pi0_sb <- (mean(pvals > lambda) / (1 - lambda)) / eta_lam
  } else {
    pi0_sb <- pi0_storey                      # fallback
  }
  pi0_sb <- min(max(pi0_sb, 1e-6), 1)
  rej_sb <- bh_reject(pvals, level = alpha / pi0_sb)
  methods_out[["SimpleBeta"]] <- list(
    pi0_hat    = pi0_sb,
    rejections = rej_sb
  )

  ## Cheng parametric
  F_hat <- which(pvals <= lambda)
  deltas <- tstats[F_hat] * sqrt(2 / n_samples)
  df <- 2 * n_samples - 2
  t_thresh <- stats::qt(1 - lambda / 2, df = df)
  Q_vals <- suppressWarnings(
  1 - stats::pt(t_thresh, df = df, ncp = deltas)
  )
  Q_vals[!is.finite(Q_vals)] <- 1
  eta_lambda <- mean(Q_vals)
  pi0_cheng <- (pi0_storey - eta_lambda) / (1 - eta_lambda)
  pi0_cheng <- min(max(pi0_cheng, 1e-6), 1)
  rej_cheng <- bh_reject(pvals, level = alpha / pi0_cheng)

  methods_out[["Cheng"]] <- list(
    pi0_hat    = pi0_cheng,
    rejections = rej_cheng
  )

  ## BUM-direct
  bum_fit   <- fit_bum(pvals)
  pi0_bum <- as.numeric(bum_fit$w0)
  pi0_bum <- min(max(pi0_bum, 1e-6), 1)
  rej_bum <- bh_reject(pvals, level = alpha / pi0_bum)

  methods_out[["BUMdirect"]] <- list(
    pi0_hat    = pi0_bum,
    rejections = rej_bum
  )

  ## BUM + model
  feat_vec  <- extract_bum_features(bum_fit, pvals)
  eta_hat   <- feat_vec["x4"]
  if (eta_hat > 1) eta_hat <- 1
  if (eta_hat < 1e-6) eta_hat <- 1e-6
  if (!bum_fit$converged || !is.finite(eta_hat) || eta_hat <= 0) {
    pi0_bum_model <- bum_fit$w0
  } else {
    pi0_bum_model <- bum_fit$w0/ eta_hat
  }
  pi0_bum_model <- min(max(pi0_bum_model, 1e-6), 1)
  rej_bum_model <- bh_reject(pvals, level = alpha / pi0_bum_model)

  methods_out[["BUMmodel"]] <- list(
    pi0_hat    = pi0_bum_model,
    rejections = rej_bum_model
  )
  
  ## scoutFDR method (BUM + calibrated correction, cross-fitted)
  cf <- adaptive_fdr(
    pvals = pvals,
    alpha = alpha,
    use_conformal = T,
    cqr_fit = cqr_fit,
    cal_df = cal_df,
    feature_names = feature_names
  )

  # Convert indices to logical vector of length m
  rej_cf <- logical(m)
  if (!is.null(cf$rejections) && length(cf$rejections) > 0L) {
    rej_cf[cf$rejections] <- TRUE
  }

  # Define a scalar pi0_hat from the fold-wise upper bounds
  pi0_hat_cf <- if (is.null(cf$U_fold) || all(is.na(cf$U_fold))) {
    1
  } else {
    min(cf$U_fold, na.rm = TRUE)
  }
  pi0_hat_cf <- min(max(pi0_hat_cf, 1e-6), 1)

  methods_out[["scoutFDR"]] <- list(
    pi0_hat    = pi0_hat_cf,
    rejections = rej_cf
  )
  
  ## Oracle method: uses true π0 but still runs BH on p-values.
  rej_oracle <- bh_reject(pvals, level = alpha / true_pi0)
  methods_out[["Oracle"]] <- list(
    pi0_hat    = true_pi0,
    rejections = rej_oracle
  )
  
  ## Turn everything into a data.frame of performance metrics.
  res <- lapply(names(methods_out), function(mname) {
    obj <- methods_out[[mname]]
    rej <- obj$rejections
    
    V <- sum(rej & is_null)
    R <- sum(rej)
    
    FDP <- if (R == 0L) 0 else V / R
    power <- if (m1 == 0L) NA_real_ else sum(rej & !is_null) / m1
    
    data.frame(
      method      = mname,
      pi0_hat     = obj$pi0_hat,
      FDR         = FDP,
      power       = power,
      discoveries = R,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, res)
}
