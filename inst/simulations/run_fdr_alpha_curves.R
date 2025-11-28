

# ============================================================
# 6) Additional FDR–alpha curves (continuous, more scenarios)
# ============================================================
suppressPackageStartupMessages({
  have_gg <- requireNamespace("ggplot2", quietly = TRUE)
  have_qr_local <- requireNamespace("quantreg", quietly = TRUE)
})

set.seed(5252)

# Utilities (define if not already present)
if (!exists("bh", mode = "function")) {
  bh <- function(p, alpha) {
    m <- length(p); if (m == 0 || alpha <= 0) return(integer(0))
    o <- order(p); ps <- p[o]; thr <- (seq_len(m) * alpha) / m
    k <- max(which(ps <= thr), 0)
    if (k == 0) integer(0) else sort(o[seq_len(k)])
  }
}
if (!exists("fdp", mode = "function")) {
  fdp <- function(rej_idx, null_idx) {
    if (length(rej_idx) == 0) return(0)
    v <- sum(rej_idx %in% null_idx)
    v / length(rej_idx)
  }
}
if (!exists("exceed_prob", mode = "function")) {
  exceed_prob <- function(fdp_vec, alpha) {
    mean(fdp_vec > alpha, na.rm = TRUE)
  }
}
if (!exists("run_scout", mode = "function")) {
  run_scout <- function(p, alpha, K, delta_ecdf, H, mode = c("safe","both","conformal"),
                        cqr_fit = NULL, cal_df = NULL, feature_names = NULL,
                        ood_params = list(alpha_cqr = 0.05, max_weight = 4)) {
    mode <- match.arg(mode)
    res <- adaptive_fdr(
      p = p, alpha = alpha, K = K,
      delta_ecdf = delta_ecdf, H = H,
      mode = mode, cqr_fit = cqr_fit,
      cal_df = cal_df, feature_names = feature_names,
      ood_params = ood_params
    )
    list(rejections = res$rejections)
  }
}

# Additional scenarios: emphasize low n and low pi0
# (IDs L1..L15 so they do not collide with earlier S1..S6)
scenarios_extra <- data.frame(
  scen_id = paste0("L", 1:15),
  n    = c(3, 3, 3, 5, 5, 5, 5, 5, 5, 8, 8, 8, 10, 10, 10),
  pi0  = c(0.2, 0.4, 0.7, 0.2, 0.4, 0.7, 0.2, 0.4, 0.7, 0.2, 0.4, 0.7, 0.2, 0.4, 0.7),
  power= c(0.20,0.20,0.20,0.20,0.20,0.20,0.35,0.35,0.35,0.35,0.35,0.35,0.50,0.50,0.50),
  stringsAsFactors = FALSE
)

m_curve <- 5000
alphas_curve <- c(0.01, 0.05, 0.10)
B_reps_curve <- 100
K_curve <- 2

# Conformal training helper per scenario
prepare_conformal_cont <- function(n, pi0, power, m, tau, cal_n = 2000) {
  if (!have_qr_local) return(NULL)
  base_pi0 <- function(p) storey_median_lambda(p, c(0.6,0.7,0.8))
  cal_df <- generate_calibration_data(
    n_sim = cal_n,
    base_pi0 = base_pi0,
    feature_fun = function(p) as.list(features_ecdf(p))
  )
  # Candidate features = intersection of calibration features and a prototype target feature set
  proto_cols <- names(as.list(features_ecdf(runif(2000))))
  cand <- setdiff(intersect(colnames(cal_df), proto_cols), "pi0")

  # Build prototype target folds using simulated p-values for the scenario
  p_probe <- simulate_pvalues_ttest(m = m, n = n, pi0 = pi0, target_power = power)
  fold_id <- sample(rep(seq_len(K_curve), length.out = length(p_probe)))
  X_target_all <- matrix(NA_real_, nrow = K_curve, ncol = length(proto_cols))
  colnames(X_target_all) <- proto_cols
  for (k in seq_len(K_curve)) {
    fk <- try(features_ecdf(p_probe[fold_id == k]), silent = TRUE)
    if (!inherits(fk, "try-error")) {
      nm <- intersect(names(fk), proto_cols)
      X_target_all[k, nm] <- as.numeric(fk[nm])
    }
  }
  X_target_all <- as.data.frame(X_target_all, stringsAsFactors = FALSE)

  keep_cal <- vapply(cand, function(nm) {
    v <- cal_df[[nm]]
    is.finite(sum(v)) && sum(is.finite(v)) >= 0.99 * nrow(cal_df) &&
      stats::sd(v, na.rm = TRUE) > 1e-10
  }, logical(1))
  keep_tar <- vapply(cand, function(nm) {
    v <- X_target_all[[nm]]
    !is.null(v) && is.numeric(v) && all(is.finite(v))
  }, logical(1))
  feat_names <- cand[keep_cal & keep_tar]
  if (length(feat_names) < 5) warning("Few robust features; conformal may be less effective for this scenario.")
  cqr_fit <- train_cqr_pi0(cal_df, feature_names = feat_names, tau = tau, model = "qr")
  list(fit = cqr_fit, cal = cal_df, feats = feat_names)
}

# Evaluate one scenario at one alpha with given settings
eval_scenario_alpha <- function(n, pi0, power, alpha,
                                conf_obj = NULL,
                                H, delta_ecdf,
                                tau_label,
                                alpha_cqr, max_weight,
                                K = K_curve, B = B_reps_curve) {
  out <- replicate(B, {
    # Generate p-values (first m0 are nulls in simulate_pvalues_ttest)
    p <- simulate_pvalues_ttest(m = m_curve, n = n, pi0 = pi0, target_power = power)
    null_idx <- seq_len(round(m_curve * pi0))

    # Baselines
    r_bh <- bh(p, alpha)
    A_med <- storey_median_lambda(p, c(0.6,0.7,0.8))
    r_st <- bh(p, alpha / A_med)

    # SCOUT variants
    res_safe <- run_scout(p, alpha, K, delta_ecdf = delta_ecdf, H = H, mode = "safe")

    res_both <- if (!is.null(conf_obj)) run_scout(
      p, alpha, K, delta_ecdf = delta_ecdf, H = H, mode = "both",
      cqr_fit = conf_obj$fit, cal_df = conf_obj$cal, feature_names = conf_obj$feats,
      ood_params = list(alpha_cqr = alpha_cqr, max_weight = max_weight)
    ) else NULL

    res_conf <- if (!is.null(conf_obj)) run_scout(
      p, alpha, K, delta_ecdf = delta_ecdf, H = H, mode = "conformal",
      cqr_fit = conf_obj$fit, cal_df = conf_obj$cal, feature_names = conf_obj$feats,
      ood_params = list(alpha_cqr = alpha_cqr, max_weight = max_weight)
    ) else NULL

    c(
      FDP_BH    = fdp(r_bh, null_idx),
      FDP_ST    = fdp(r_st, null_idx),
      FDP_safe  = fdp(res_safe$rejections, null_idx),
      FDP_both  = if (is.null(res_both)) NA_real_ else fdp(res_both$rejections, null_idx),
      FDP_conf  = if (is.null(res_conf)) NA_real_ else fdp(res_conf$rejections, null_idx)
    )
  }, simplify = "matrix")

  data.frame(
    scen_id = NA_character_, n = n, pi0 = pi0, power = power,
    alpha = alpha, method = c("BH","Storey","SCOUT_safe","SCOUT_both","SCOUT_confOnly"),
    mean_fdp = c(mean(out["FDP_BH",]),
                 mean(out["FDP_ST",]),
                 mean(out["FDP_safe",]),
                 mean(out["FDP_both",], na.rm = TRUE),
                 mean(out["FDP_conf",], na.rm = TRUE)),
    exceed = c(exceed_prob(out["FDP_BH",], alpha),
               exceed_prob(out["FDP_ST",], alpha),
               exceed_prob(out["FDP_safe",], alpha),
               exceed_prob(out["FDP_both",], alpha),
               exceed_prob(out["FDP_conf",], alpha)),
    H = paste(H, collapse = ","),
    delta_ecdf = delta_ecdf,
    tau = tau_label,
    alpha_cqr = alpha_cqr,
    max_weight = max_weight,
    stringsAsFactors = FALSE
  )
}

# --------------------------
# A) Original settings (match your previous)
# --------------------------
H_orig <- c(0.10, 0.20, 0.30)
delta_ecdf_orig <- 0.02
tau_orig <- 0.65
alpha_cqr_orig <- 0.02
max_weight_orig <- 4

results_extra <- list()
for (i in seq_len(nrow(scenarios_extra))) {
  scen <- scenarios_extra[i, ]
  cat(sprintf("\n[EXTRA original] Preparing conformal for %s (n=%d, pi0=%.2f, power=%.2f, tau=%.2f)...\n",
              scen$scen_id, scen$n, scen$pi0, scen$power, tau_orig))
  conf_obj <- try(prepare_conformal_cont(scen$n, scen$pi0, scen$power, m_curve, tau = tau_orig, cal_n = 2000), silent = TRUE)
  if (inherits(conf_obj, "try-error")) conf_obj <- NULL

  for (a in alphas_curve) {
    cat(sprintf("  %s @ alpha=%.2f (original settings)\n", scen$scen_id, a))
    res <- eval_scenario_alpha(
      n = scen$n, pi0 = scen$pi0, power = scen$power, alpha = a,
      conf_obj = conf_obj, H = H_orig, delta_ecdf = delta_ecdf_orig,
      tau_label = tau_orig, alpha_cqr = alpha_cqr_orig, max_weight = max_weight_orig
    )
    res$scen_id <- scen$scen_id
    results_extra[[length(results_extra) + 1]] <- res
  }
}
fdr_curve_extra <- do.call(rbind, results_extra)
row.names(fdr_curve_extra) <- NULL
write.csv(fdr_curve_extra, "fdr_curve_continuous_extra.csv", row.names = FALSE)
cat("\nSaved fdr_curve_continuous_extra.csv\n")

if (have_gg) {
  library(ggplot2)
  p_extra <- ggplot(fdr_curve_extra, aes(x = alpha, y = mean_fdp, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~ scen_id + paste0("n=", n, ", pi0=", pi0, ", pow=", power), scales = "free_y") +
    scale_x_continuous(breaks = alphas_curve, limits = c(min(alphas_curve), max(alphas_curve))) +
    labs(title = "Observed mean FDP vs target FDR (continuous, extra scenarios, original settings)",
         subtitle = "H={0.10,0.20,0.30}, delta_ecdf=0.02; conformal tau=0.65, alpha_cqr=0.02, max_weight=4",
         x = "Target FDR (alpha)", y = "Observed mean FDP") +
    theme_bw() + theme(legend.position = "bottom")
  ggsave("fdr_curve_continuous_extra.png", p_extra, width = 13, height = 9, dpi = 300)
  cat("Saved fdr_curve_continuous_extra.png\n")
} else {
  cat("ggplot2 not installed; skipping extra plot (original settings).\n")
}

# --------------------------
# B) Conservative settings (optional, for comparison)
# --------------------------
H_cons <- c(0.02, 0.05, 0.10)
delta_ecdf_cons <- 0.005
tau_cons <- 0.70
alpha_cqr_cons <- 0.03
max_weight_cons <- 3

results_extra_cons <- list()
for (i in seq_len(nrow(scenarios_extra))) {
  scen <- scenarios_extra[i, ]
  cat(sprintf("\n[EXTRA conservative] Preparing conformal for %s (n=%d, pi0=%.2f, power=%.2f, tau=%.2f)...\n",
              scen$scen_id, scen$n, scen$pi0, scen$power, tau_cons))
  conf_obj <- try(prepare_conformal_cont(scen$n, scen$pi0, scen$power, m_curve, tau = tau_cons, cal_n = 2000), silent = TRUE)
  if (inherits(conf_obj, "try-error")) conf_obj <- NULL

  for (a in alphas_curve) {
    cat(sprintf("  %s @ alpha=%.2f (conservative settings)\n", scen$scen_id, a))
    res <- eval_scenario_alpha(
      n = scen$n, pi0 = scen$pi0, power = scen$power, alpha = a,
      conf_obj = conf_obj, H = H_cons, delta_ecdf = delta_ecdf_cons,
      tau_label = tau_cons, alpha_cqr = alpha_cqr_cons, max_weight = max_weight_cons
    )
    res$scen_id <- scen$scen_id
    results_extra_cons[[length(results_extra_cons) + 1]] <- res
  }
}
fdr_curve_extra_cons <- do.call(rbind, results_extra_cons)
row.names(fdr_curve_extra_cons) <- NULL
write.csv(fdr_curve_extra_cons, "fdr_curve_continuous_extra_conservative.csv", row.names = FALSE)
cat("\nSaved fdr_curve_continuous_extra_conservative.csv\n")

if (have_gg) {
  library(ggplot2)
  p_extra_cons <- ggplot(fdr_curve_extra_cons, aes(x = alpha, y = mean_fdp, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~ scen_id + paste0("n=", n, ", pi0=", pi0, ", pow=", power), scales = "free_y") +
    scale_x_continuous(breaks = alphas_curve, limits = c(min(alphas_curve), max(alphas_curve))) +
    labs(title = "Observed mean FDP vs target FDR (continuous, extra scenarios, conservative settings)",
         subtitle = "H={0.02,0.05,0.10}, delta_ecdf=0.005; conformal tau=0.70, alpha_cqr=0.03, max_weight=3",
         x = "Target FDR (alpha)", y = "Observed mean FDP") +
    theme_bw() + theme(legend.position = "bottom")
  ggsave("fdr_curve_continuous_extra_conservative.png", p_extra_cons, width = 13, height = 9, dpi = 300)
  cat("Saved fdr_curve_continuous_extra_conservative.png\n")
} else {
  cat("ggplot2 not installed; skipping extra plot (conservative settings).\n")
}

# Notes
# - Leaves earlier results intact; appends two new sets (original vs conservative) for the extra scenarios.
# - Uses simulate_pvalues_ttest to generate p-values (first m*pi0 indices are nulls).
# - Conformal runs only if quantreg is installed; otherwise the 'both'/'conformal' modes gracefully fall back to safe-only when conf_obj is NULL.
# - Increase B_reps_curve for tighter Monte Carlo error, at the cost of runtime.
