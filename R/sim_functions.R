

power_for_delta <- function(delta, n, alpha = 0.05) {
  nu <- 2 * (n - 1)
  tcrit <- stats::qt(1 - alpha / 2, df = nu)
  z1 <- tcrit - delta * sqrt(n / 2)
  z2 <- -tcrit - delta * sqrt(n / 2)
  1 - stats::pnorm(z1) + stats::pnorm(z2)
}


avg_power_for_mu <- function(mu, n, alpha = 0.05, sd_ratio = 0.3) {
  sd <- sd_ratio * mu
  f <- function(delta) {
    stats::dnorm(delta, mean = mu, sd = sd) * 
      power_for_delta(delta, n, alpha = alpha)
  }
  stats::integrate(f, lower = -Inf, upper = Inf)$value
}


calibrate_mu <- function(target_power, n, alpha = 0.05, sd_ratio = 0.3) {
  f_root <- function(mu)
    avg_power_for_mu(mu, n, alpha = alpha, sd_ratio = sd_ratio) - target_power

  stats::uniroot(f_root, interval = c(0.01, 5), tol = 1e-3)$root
}

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

simulate_pvalues_ttest <- function(pi0, target_power, m, n) {

  # Split into nulls and non-nulls
  m0 <- round(m * pi0)
  m1 <- m - m0

  # POWER → effect size distribution parameters
  z_power <- stats::qnorm(pmin(pmax(target_power, 1e-4), 1 - 1e-4))
  z_alpha <- stats::qnorm(0.975)
  delta_mean <- (z_alpha + z_power) / sqrt(n / 2)
  delta_sd   <- abs(delta_mean) * 0.3

  # 1) NULLS (delta = 0)
  # Generate all null samples at once
  if (m0 > 0) {
    X1_0 <- matrix(stats::rnorm(m0 * n, 0, 1), nrow = m0)
    X2_0 <- matrix(stats::rnorm(m0 * n, 0, 1), nrow = m0)

    mean_diff_0 <- matrixStats::rowMeans2(X1_0) - matrixStats::rowMeans2(X2_0)

    var1_0 <- matrixStats::rowVars(X1_0)
    var2_0 <- matrixStats::rowVars(X2_0)

    pooled_sd_0 <- sqrt(((n - 1) * var1_0 + (n - 1) * var2_0) / (2*n - 2))

    tvals_0 <- mean_diff_0 / (pooled_sd_0 * sqrt(2/n))

    p0 <- 2 * (1 - stats::pt(abs(tvals_0), df = 2*n - 2))
  } else {
    p0 <- numeric(0)
  }

  # NON-NULLS (random delta ~ Normal(delta_mean, delta_sd))
  if (m1 > 0) {

    delta <- stats::rnorm(m1, delta_mean, delta_sd)

    # Generate alternative sample matrices
    X1_1 <- matrix(stats::rnorm(m1 * n, 0, 1), nrow = m1)
    X2_1 <- matrix(stats::rnorm(m1 * n, delta, 1), nrow = m1)

    mean_diff_1 <- matrixStats::rowMeans2(X1_1) - matrixStats::rowMeans2(X2_1)

    var1_1 <- matrixStats::rowVars(X1_1)
    var2_1 <- matrixStats::rowVars(X2_1)

    pooled_sd_1 <- sqrt(((n - 1) * var1_1 + (n - 1) * var2_1) / (2*n - 2))

    tvals_1 <- mean_diff_1 / (pooled_sd_1 * sqrt(2/n))

    p1 <- 2 * (1 - stats::pt(abs(tvals_1), df = 2*n - 2))
  } else {
    p1 <- numeric(0)
  }

  # Combine and shuffle
  p <- c(p0, p1)
  sample(p)
}

simulate_dataset <- function(pi0, target_power, m, n, effect_grid) {
  # Determine null / non-null split
  m0 <- round(m * pi0)
  m1 <- m - m0

  # Extract effect size parameters from effect_grid
  idx <- which(effect_grid$n == n & effect_grid$target_power == target_power)
  if (length(idx) != 1L)
    stop("Effect grid does not contain unique match for (n, target_power).")

  delta_mean <- effect_grid$delta_mean[idx]
  delta_sd   <- effect_grid$delta_sd[idx]

  # Generate NULL test statistics
  if (m0 > 0) {
    X1_0 <- matrix(stats::rnorm(m0 * n, 0, 1), nrow = m0)
    X2_0 <- matrix(stats::rnorm(m0 * n, 0, 1), nrow = m0)

    mean_diff_0 <- matrixStats::rowMeans2(X1_0) - matrixStats::rowMeans2(X2_0)
    var1_0 <- matrixStats::rowVars(X1_0)
    var2_0 <- matrixStats::rowVars(X2_0)

    pooled_sd_0 <- sqrt(((n - 1) * var1_0 + (n - 1) * var2_0) / (2*n - 2))
    tstat_0 <- mean_diff_0 / (pooled_sd_0 * sqrt(2/n))

    p0 <- 2 * (1 - stats::pt(abs(tstat_0), df = 2*n - 2))
  } else {
    tstat_0 <- numeric(0)
    p0 <- numeric(0)
  }

  # Generate NON-NULL test statistics
  if (m1 > 0) {
    delta_vec <- stats::rnorm(m1, delta_mean, delta_sd)

    X1_1 <- matrix(stats::rnorm(m1 * n, 0, 1), nrow = m1)
    X2_1 <- matrix(stats::rnorm(m1 * n, delta_vec, 1), nrow = m1)

    mean_diff_1 <- matrixStats::rowMeans2(X1_1) - matrixStats::rowMeans2(X2_1)
    var1_1 <- matrixStats::rowVars(X1_1)
    var2_1 <- matrixStats::rowVars(X2_1)

    pooled_sd_1 <- sqrt(((n - 1) * var1_1 + (n - 1) * var2_1) / (2*n - 2))
    tstat_1 <- mean_diff_1 / (pooled_sd_1 * sqrt(2/n))

    p1 <- 2 * (1 - stats::pt(abs(tstat_1), df = 2*n - 2))
  } else {
    tstat_1 <- numeric(0)
    p1 <- numeric(0)
  }

  # Combine results
  tstats <- c(tstat_0, tstat_1)
  pvals  <- c(p0, p1)
  is_null <- c(rep(TRUE, m0), rep(FALSE, m1))

  # Shuffle to avoid ordering artifacts
  perm <- sample.int(m)
  list(
    pvals   = pvals[perm],
    tstats  = tstats[perm],
    is_null = is_null[perm]
  )
}


#' Run multiple FDR estimators on a vector of p-values
#'
#' @description
#' This function is intended *only* for simulation studies.  
#' It runs:
#'
#' - BH  
#' - Storey (λ = 0.5)  
#' - Pounds–Cheng  
#' - Simple Beta estimator  
#' - Cheng parametric  
#' - BUM-direct  
#' - BUM + tail-corrected  
#' - SCOUT-FDR  
#' - Oracle (requires true π₀)
#'
#' @param pvals Numeric p-values.
#' @param tstats t-statistics (for Cheng).
#' @param n_samples Sample size per group.
#' @param cal_df Calibration data for SCOUT-FDR.
#' @param cqr_fit CQR model.
#' @param feature_names Feature names for SCOUT-FDR.
#' @param is_null Logical vector marking null hypotheses.
#' @param alpha FDR target.
#' @param true_pi0 True null proportion (Oracle).
#'
#' @return A data.frame summarizing:
#'   - method
#'   - pi0_hat
#'   - FDR
#'   - power
#'   - discoveries
#'
#' @export
run_all_methods <- function(
    pvals,
    tstats,
    n_samples,
    cal_df,
    cqr_fit,
    feature_names,
    is_null,
    alpha,
    true_pi0
) {

  m0 <- sum(is_null)
  m1 <- length(is_null) - m0
  methods_out <- list()

  # ----------------------------------------------
  # BH
  methods_out$BH <- list(
    pi0_hat = 1,
    rejections = bh_reject(pvals, alpha)
  )

  # ----------------------------------------------
  # Storey (λ=0.5)
  lambda <- 0.5
  pi0_storey <- min(max(mean(pvals > lambda) / (1 - lambda), 1e-6), 1)
  methods_out$Storey <- list(
    pi0_hat = pi0_storey,
    rejections = bh_reject(pvals, alpha / pi0_storey)
  )

  # ----------------------------------------------
  # Pounds–Cheng (average Storey λ)
  lambda_grid <- seq(0.5, 0.9, by = 0.05)
  pi0_pc_vec <- sapply(lambda_grid, function(l)
    mean(pvals > l) / (1 - l)
  )
  pi0_pc <- min(max(mean(pi0_pc_vec), 1e-6), 1)
  methods_out$PoundsCheng <- list(
    pi0_hat = pi0_pc,
    rejections = bh_reject(pvals, alpha / pi0_pc)
  )

  # ----------------------------------------------
  # Simple beta model
  p_low <- pvals[pvals <= 0.25]
  if (length(p_low) >= 50) {
    beta_hat <- -mean(log(p_low))
    eta_lam  <- (1 - lambda)^beta_hat
    pi0_sb <- pi0_storey / eta_lam
  } else {
    pi0_sb <- pi0_storey
  }
  pi0_sb <- min(max(pi0_sb, 1e-6), 1)
  methods_out$SimpleBeta <- list(
    pi0_hat = pi0_sb,
    rejections = bh_reject(pvals, alpha / pi0_sb)
  )

  # ----------------------------------------------
  # Cheng parametric
  F_hat <- which(pvals <= lambda)
  deltas <- tstats[F_hat] * sqrt(2 / n_samples)
  df <- 2 * n_samples - 2
  t_thresh <- stats::qt(1 - lambda / 2, df)
  Q_vals <- suppressWarnings(1 - stats::pt(t_thresh, df = df, ncp = deltas))
  Q_vals[!is.finite(Q_vals)] <- 1
  eta_lambda <- mean(Q_vals)
  pi0_cheng <- min(max((pi0_storey - eta_lambda) / (1 - eta_lambda), 1e-6), 1)
  methods_out$Cheng <- list(
    pi0_hat = pi0_cheng,
    rejections = bh_reject(pvals, alpha / pi0_cheng)
  )

  # ----------------------------------------------
  # BUM direct
  bum_fit <- fit_bum(pvals)
  pi0_bum <- min(max(bum_fit$w0 %||% 1, 1e-6), 1)
  methods_out$BUMdirect <- list(
    pi0_hat = pi0_bum,
    rejections = bh_reject(pvals, alpha / pi0_bum)
  )

  # ----------------------------------------------
  # BUM + correction
  bum_feats <- extract_bum_features(bum_fit, pvals)
  eta_hat <- bum_feats["bum_4"]
  if (!is.finite(eta_hat) || eta_hat <= 0) eta_hat <- 1
  pi0_bm <- min(max(pi0_bum / eta_hat, 1e-6), 1)
  methods_out$BUMmodel <- list(
    pi0_hat = pi0_bm,
    rejections = bh_reject(pvals, alpha / pi0_bm)
  )

  # ----------------------------------------------
  # SCOUT-FDR
  cf <- adaptive_fdr(
    pvals = pvals,
    alpha = alpha,
    use_conformal = TRUE,
    cqr_fit = cqr_fit,
    cal_df = cal_df,
    feature_names = feature_names
  )
  rej_cf <- logical(length(pvals))
  if (length(cf$rejections) > 0)
    rej_cf[cf$rejections] <- TRUE

  pi0_hat_cf <- min(max(min(cf$U_fold, na.rm = TRUE), 1e-6), 1)

  methods_out$scoutFDR <- list(
    pi0_hat = pi0_hat_cf,
    rejections = rej_cf
  )

  # ----------------------------------------------
  # Oracle
  methods_out$Oracle <- list(
    pi0_hat = true_pi0,
    rejections = bh_reject(pvals, alpha / true_pi0)
  )

  # ----------------------------------------------
  # Compile results
  df <- lapply(names(methods_out), function(mname) {
    obj <- methods_out[[mname]]
    rej <- obj$rejections

    V <- sum(rej & is_null)
    R <- sum(rej)
    FDP <- if (R == 0) 0 else V / R
    power <- if (m1 == 0) NA_real_ else sum(rej & !is_null) / m1

    data.frame(
      method = mname,
      pi0_hat = obj$pi0_hat,
      FDR = FDP,
      power = power,
      discoveries = R,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, df)
}
