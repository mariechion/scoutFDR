# ---------- Load required packages and functions ----------
library(scoutFDR)
library(quantreg) # Required for conformal training
library(ggplot2)
library(dplyr)
library(purrr)

# ---------- Utility functions ----------
# Benjamini-Hochberg procedure
bh <- function(p, alpha) {
  m <- length(p); o <- order(p); ps <- p[o]; k <- seq_len(m)
  thr <- k * alpha / m
  kmax <- max(which(ps <= thr), 0)
  if (kmax == 0) integer(0) else sort(o[ps <= ps[kmax]])
}
# FDP calculation
fdp <- function(rej, null_idx) if (length(rej) == 0) 0 else sum(rej %in% null_idx) / length(rej)
# FDP exceedance probability i.e. P(FDP > alpha)
exceed_prob <- function(fdp_vec, alpha) mean(fdp_vec > alpha, na.rm = TRUE)
# Pinball loss function allowing for weights
pinball_loss <- function(y, yhat, tau, weights = NULL) {
  u <- y - yhat
  loss <- ifelse(u >= 0, tau * u, (1 - tau) * (-u))
  ifelse(is.null(weights), mean(loss), mean(loss * weights))
}

# ---------- Scenario and simulation setup ----------
m_curve <- 5000 # total number of hypotheses
alphas_curve <- c(0.01, 0.05, 0.10) # target FDR levels
B_reps_curve <- 100 # number of repetitions per scenario
K_curve <- 2 # number of folds for fold-safe SCOUT

# Construct scenario grid
scenarios_extra <- data.frame(
  scen_id = paste0("L", 1:15),
  n    = c(3, 3, 3, 5, 5, 5, 5, 5, 5, 8, 8, 8, 10, 10, 10),
  pi0  = c(0.2, 0.4, 0.7, 0.2, 0.4, 0.7, 0.2, 0.4, 0.7, 0.2, 0.4, 0.7, 0.2, 0.4, 0.7),
  power= c(0.20,0.20,0.20,0.20,0.20,0.20,0.35,0.35,0.35,0.35,0.35,0.35,0.50,0.50,0.50),
  stringsAsFactors = FALSE
)

set.seed(17)

# ---------- Conformal training (global model) ----------
# Initialise objects
cqr_fit <- NULL
cal_df <- NULL
feature_names <- NULL
base_pi0 <- function(pp) scoutFDR:::storey_median_lambda(pp, c(0.6,0.7,0.8))

# Generate calibration data
cal_df <- generate_calibration_data(n_sim = 5000,
  base_pi0 = base_pi0,
  feature_fun = function(pp) as.list(features_ecdf(pp))
)

# Train QR at a moderately conservative tau
cqr_fit <- train_cqr_pi0(cal_df, feature_names = "auto", tau = 0.70, model = "qr")


# ---------- SCOUT runner ----------
run_scout_once <- function(p, alpha, K, H, delta_ecdf, use_conformal, cqr_fit, cal_df, feature_names) {
  scoutFDR::adaptive_fdr(
    p, alpha = alpha, K = K,
    delta_ecdf = delta_ecdf, H = H,
    use_conformal = use_conformal,
    cqr_fit = if (use_conformal) cqr_fit else NULL,
    cal_df   = if (use_conformal) cal_df else NULL,
    feature_names = if (use_conformal) feature_names else NULL
  )
}

# ---------- Evaluate one scenario at one alpha ----------
eval_scenario_alpha <- function(n, pi0, power, alpha,
  conf_obj = NULL,
  H, delta_ecdf,
  tau_label,
  alpha_cqr, max_weight,
  K = K_curve, B = B_reps_curve) {
    out <- replicate(B, {
      # Simulate p-values
      p <- scoutFDR:::simulate_pvalues_ttest(m = m_curve, n = n, pi0 = pi0, target_power = power)
      # Identify nulls
      null_idx <- seq_len(round(m_curve * pi0))
      # Apply methods
      r_bh <- bh(p, alpha)
      A_med <- scoutFDR:::storey_median_lambda(p, c(0.6,0.7,0.8))
      r_st  <- bh(p, alpha / A_med)
      # SCOUT safe
      res_safe <- run_scout_once(p, alpha, K, H, delta_ecdf, FALSE, NULL, NULL, NULL)
      # SCOUT both 
      res_both <- if (!is.null(conf_obj)) scoutFDR::adaptive_fdr(
        p, alpha = alpha, K = K, delta_ecdf = delta_ecdf, H = H,
        mode = "both", use_conformal = TRUE,
        cqr_fit = conf_obj$fit, cal_df = conf_obj$cal, feature_names = conf_obj$feats,
        ood_params = list(alpha_cqr = alpha_cqr, max_weight = max_weight)
      ) else NULL
      # SCOUT conformal only
      res_conf <- if (!is.null(conf_obj)) scoutFDR::adaptive_fdr(
        p, alpha = alpha, K = K, delta_ecdf = delta_ecdf, H = H,
        mode = "conformal", use_conformal = TRUE,
        cqr_fit = conf_obj$fit, cal_df = conf_obj$cal, feature_names = conf_obj$feats,
        ood_params = list(alpha_cqr = alpha_cqr, max_weight = max_weight)
      ) else NULL
      # Compute FDPs
      c(
        FDP_BH    = fdp(r_bh, null_idx),
        FDP_ST    = fdp(r_st, null_idx),
        FDP_safe  = fdp(res_safe$rejections, null_idx),
        FDP_both  = if (is.null(res_both)) NA_real_ else fdp(res_both$rejections, null_idx),
        FDP_conf  = if (is.null(res_conf)) NA_real_ else fdp(res_conf$rejections, null_idx)
      )
    }, simplify = "matrix")
    # Summarise results
    data.frame(
      scen_id = NA_character_, 
      n = n, 
      pi0 = pi0, 
      power = power,
      alpha = alpha, 
      method = c("BH","Storey","SCOUT_safe","SCOUT_both","SCOUT_confOnly"),
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
  
# ---------- Feature importance (global + scenario-weighted) ----------
# Permutation feature importance based on increase in pinball loss
compute_feature_importance_perm <- function(cqr_fit, cal_df, weights = NULL,
  repeats = 5, seed = 123) {
  # Sanity checks
  stopifnot(!is.null(cqr_fit), !is.null(cal_df))
      
  # Prepare data
  feat_names <- cqr_fit$feature_names
  keep_rows <- stats::complete.cases(cal_df[, c("pi0", feat_names), drop = FALSE])
  df <- cal_df[keep_rows, , drop = FALSE]
  y <- as.numeric(df$pi0)
  X <- as.matrix(df[, feat_names, drop = FALSE])
      
  mu <- cqr_fit$mu
  sd <- cqr_fit$sd
  if (is.null(mu) || is.null(sd)) { 
     mu <- rep(0, length(feat_names)); 
    sd <- rep(1, length(feat_names)) 
  }
  sd[sd < 1e-8] <- 1
      
  Xs <- scale(X, center = mu, scale = sd)
  if (!is.null(cqr_fit$vars)) colnames(Xs) <- cqr_fit$vars
      
  tau <- cqr_fit$tau %||% 0.5
      
  # Identify CQR model
  model_type <- dplyr::case_when(
    !is.null(cqr_fit$rqfit) ~ "rq",
    !is.null(cqr_fit$qrf)   ~ "qrf",
    !is.null(cqr_fit$lmfit) ~ "lm",
    TRUE                    ~ "custom"
  )
      
  # Baseline prediction based on type of CQR model
  switch(model_type,
    rq = {
      b <- coef(cqr_fit$rqfit)
      yhat0 <- as.numeric(b[1] + Xs %*% as.numeric(b[-1]))
    },
    qrf = {
      yhat0 <- as.numeric(predict(cqr_fit$qrf, X, what = tau))
    },
    lm = {
      yhat0 <- as.numeric(predict(cqr_fit$lmfit, as.data.frame(X)))
    },
    custom = {
      yhat0 <- apply(X, 1, function(xx) predict_qtau(cqr_fit, setNames(xx, feat_names)))
    }
  ) 

  # Validate and align weights if provided
  if (!is.null(weights)) {
    weights <- as.numeric(weights)[keep_rows]
    if (length(weights) != nrow(X))
      stop("Weights length does not match filtered calibration rows.")
  }
  # Baseline unified pinball loss (weighted or unweighted)
  L0 <- pinball_loss(y, yhat0, tau, weights = weights)
      
  # Permutation importance calculation
  set.seed(seed)
  perm_importance <- purrr::map_dbl(seq_along(feat_names), function(j) {
    deltas <- purrr::map_dbl(seq_len(repeats), function(r) {
      Xs_perm <- Xs
      Xs_perm[, j] <- sample(Xs_perm[, j])
      switch(model_type,
       rq = {
          b <- coef(cqr_fit$rqfit)
          yhat_p <- as.numeric(b[1] + Xs_perm %*% as.numeric(b[-1]))
        },
        qrf = {
          Xp <- X
          Xp[, j] <- sample(Xp[, j])
          yhat_p <- as.numeric(predict(cqr_fit$qrf, Xp, what = tau))
        },
        lm = {
          Xp <- X
          Xp[, j] <- sample(Xp[, j])
          yhat_p <- as.numeric(predict(cqr_fit$lmfit, as.data.frame(Xp)))
        },
        custom = {
          Xp <- X
          Xp[, j] <- sample(Xp[, j])
          yhat_p <- apply(Xp, 1, function(xx) predict_qtau(cqr_fit, setNames(xx, feat_names)))
        }
      )
      Lp <- pinball_loss(y, yhat_p, tau, weights)
      Lp - L0
    })
  mean(deltas)
  })
      
  data.frame(
    feature = feat_names,
    importance = perm_importance,
    type = if (is.null(weights)) "perm_global" else "perm_weighted",
    stringsAsFactors = FALSE
  )
}
    
# Estimate density-ratio weights: logistic regression cal vs target features
estimate_dr_weights <- function(cal_df, feat_for_weight, X_target_all, clip = c(1/3, 3)) {
  # Sanity checks
  feat_for_weight <- intersect(feat_for_weight, intersect(colnames(cal_df), colnames(X_target_all)))
  if (length(feat_for_weight) == 0) return(rep(1, nrow(cal_df)))
  
  # Prepare data
  X_cal <- as.matrix(cal_df[, feat_for_weight, drop = FALSE])
  X_tar <- as.matrix(X_target_all[, feat_for_weight, drop = FALSE])
  n_c <- nrow(X_cal)
  n_t <- nrow(X_tar)
  if (n_c == 0 || n_t == 0) return(rep(1, nrow(cal_df)))
  if (n_c > n_t) X_cal <- X_cal[sample.int(n_c, n_t), , drop = FALSE]
  if (n_t > n_c) X_tar <- X_tar[sample.int(n_t, n_c), , drop = FALSE]
     
  # Standardize features
  X_all <- rbind(X_cal, X_tar)
  mu <- colMeans(X_all)
  sdv <- apply(X_all, 2, sd) 
  sdv[sdv < 1e-8] <- 1
  Xc <- scale(X_cal, center = mu, scale = sdv)
  Xt <- scale(X_tar, center = mu, scale = sdv)
  
  # Fit logistic regression
  df_lr <- data.frame(z = c(rep(0L, nrow(Xc)), rep(1L, nrow(Xt))), rbind(Xc, Xt))
  colnames(df_lr) <- c("z", feat_for_weight)
  gl <- suppressWarnings(try(stats::glm(z ~ ., data = df_lr, family = stats::binomial()), silent = TRUE))
  if (inherits(gl, "try-error")) return(rep(1, nrow(cal_df)))
  
  # Predict on calibration set
  Xc_full <- as.matrix(cal_df[, feat_for_weight, drop = FALSE])
  Xc_full_s <- scale(Xc_full, center = mu, scale = sdv)
  p_cal <- suppressWarnings(try(stats::predict(gl, newdata = as.data.frame(Xc_full_s), type = "response"), silent = TRUE))
  if (inherits(p_cal, "try-error")) return(rep(1, nrow(cal_df)))
  
  # Compute weights and clip
  w_raw <- p_cal / pmax(1 - p_cal, 1e-12)
  pmin(pmax(w_raw, clip[1], na.rm = TRUE), clip[2], na.rm = TRUE)
}
    
# Build per-scenario target features (fold-safe) for importance weighting
build_target_features <- function(p, K, feature_names) {
  fold_id <- sample(rep(seq_len(K), length.out = length(p)))
  X_target_all <- matrix(NA_real_, nrow = K, ncol = length(feature_names))
  colnames(X_target_all) <- feature_names
  for (k in seq_len(K)) {
    fk <- try(features_ecdf(p[fold_id == k]), silent = TRUE)
    if (!inherits(fk, "try-error")) {
      nm <- intersect(names(fk), feature_names)
      X_target_all[k, nm] <- as.numeric(fk[nm])
    }
  }
  as.data.frame(X_target_all, stringsAsFactors = FALSE)
}
    
# ---------- Run extra scenarios: original and conservative ----------
    
run_block <- function(scenarios, H, delta_ecdf, tau_label, alpha_cqr, max_weight, fname_prefix) {
  results <- list()
  conf_obj <- if (!is.null(cqr_fit)) list(fit = cqr_fit, cal = cal_df, feats = feature_names) else NULL
      
  for (i in seq_len(nrow(scenarios))) {
    scen <- scenarios[i, ]
    for (a in alphas_curve) {
      res <- eval_scenario_alpha(
        n = scen$n, pi0 = scen$pi0, power = scen$power, alpha = a,
        conf_obj = conf_obj,
        H = H, delta_ecdf = delta_ecdf,
        tau_label = tau_label, alpha_cqr = alpha_cqr, max_weight = max_weight
      )
      res$scen_id <- scen$scen_id
      results[[length(results) + 1]] <- res
    }
  }
  df <- do.call(rbind, results)
  write.csv(df, file = sprintf("inst/simulations/%s.csv", fname_prefix), row.names = FALSE)

  p_plot <- ggplot2::ggplot(df, aes(x = alpha, y = mean_fdp, color = method, group = method)) +
    geom_line(size = 1) + 
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~ scen_id + paste0("n=", n, ", pi0=", pi0, ", pow=", power), scales = "free_y") +
    scale_x_continuous(breaks = alphas_curve, limits = c(min(alphas_curve), max(alphas_curve))) +
    labs(title = sprintf("Observed mean FDP vs target FDR (%s)", fname_prefix),
      subtitle = sprintf("H={%s}, delta_ecdf=%.3f; conformal tau=%.2f, alpha_cqr=%.2f, max_weight=%d",
      paste(H, collapse = ","), delta_ecdf, tau_label, alpha_cqr, max_weight),
      x = "Target FDR (alpha)", y = "Observed mean FDP") +
    theme_bw() + 
    theme(legend.position = "bottom")
  ggsave(sprintf("inst/simulations/%s.png", fname_prefix), p_plot, width = 13, height = 9, dpi = 300)
  invisible(df)
  }
    
# Original settings
df_extra_orig <- run_block(
  scenarios = scenarios_extra,
  H = c(0.10, 0.20, 0.30),
  delta_ecdf = 0.02,
  tau_label = 0.65,
  alpha_cqr = 0.02,
  max_weight = 4,
  fname_prefix = "fdr_curve_continuous_extra"
)
    
# Conservative settings
df_extra_cons <- run_block(
  scenarios = scenarios_extra,
  H = c(0.02, 0.05, 0.10),
  delta_ecdf = 0.005,
  tau_label = 0.70,
  alpha_cqr = 0.03,
  max_weight = 3,
  fname_prefix = "fdr_curve_continuous_extra_conservative"
)
    
# ---------- Feature importance outputs ----------
    
# Global permutation importance on calibration
imp_global <- compute_feature_importance_perm(
  cqr_fit = cqr_fit,
  cal_df  = cal_df,
  weights = NULL,
  repeats = 100
)
write.csv(imp_global, "inst/simulations/feature_importance_cont_perm_global.csv", row.names = FALSE)
    
ig <- imp_global[order(-imp_global$importance), ]
topK <- min(30, nrow(ig))
ig_top <- ig[1:topK, , drop = FALSE]
ig_top$feature <- factor(ig_top$feature, levels = rev(ig_top$feature))
g <- ggplot2::ggplot(ig_top, aes(x = feature, y = importance)) +
    geom_col(fill = "#2c7fb8") +
    coord_flip() + 
    theme_minimal(base_size = 12) +
    labs(title = "Permutation importance (global calibration)",
          subtitle = sprintf("tau=%.2f; Δ pinball loss (top %d)", cqr_fit$tau %||% NA_real_, topK),
          x = NULL, 
          y = "Δ pinball loss")
ggsave("inst/simulations/feature_importance_cont_perm_global.png", g, width = 8, height = 10, dpi = 300)
    
    
# Scenario-weighted permutation importance
for (i in seq_len(nrow(scenarios_extra))) {
  scen <- scenarios_extra[i, ]
      
  p_tar <- scoutFDR:::simulate_pvalues_ttest(m = m_curve, n = scen$n, pi0 = scen$pi0, target_power = scen$power)
  X_tar <- build_target_features(p_tar, K = K_curve, feature_names = cqr_fit$feature_names)
      
  keep_cols <- vapply(X_tar, function(v) all(is.finite(v)), logical(1))
  feat_ok <- intersect(cqr_fit$feature_names, names(X_tar)[keep_cols])
  if (length(feat_ok) < 5) next
      
  w <- estimate_dr_weights(
    cal_df         = cal_df,
    feat_for_weight= feat_ok,
    X_target_all   = X_tar[, feat_ok, drop = FALSE],
    clip           = c(1/3, 3)
  )
      
  imp_w <- compute_feature_importance_perm(
    cqr_fit = cqr_fit,
    cal_df  = cal_df,
    weights = w,
    repeats = 100
  )
      
  ofile <- sprintf("inst/simulations/feature_importance_cont_perm_weighted_%s.csv", scen$scen_id)
  write.csv(imp_w, ofile, row.names = FALSE)
      
  # Plot weighted importance
  iw <- imp_w[order(-imp_w$importance), ]
  topK <- min(25, nrow(iw))
  iw_top <- iw[1:topK, , drop = FALSE]
  iw_top$feature <- factor(iw_top$feature, levels = rev(iw_top$feature))
      
  gw <- ggplot2::ggplot(iw_top, aes(x = feature, y = importance)) +
      geom_col(fill = "#41ab5d") +
      coord_flip() + 
      theme_minimal(base_size = 12) +
      labs(title = sprintf("Permutation importance (weighted) – %s", scen$scen_id),
          subtitle = sprintf("n=%d, pi0=%.2f, power=%.2f; τ=%.2f", scen$n, scen$pi0, scen$power, cqr_fit$tau %||% NA_real_),
          x = NULL, 
          y = "Δ pinball loss (weighted)")
  ggsave(sprintf("inst/simulations/feature_importance_cont_perm_weighted_%s.png", scen$scen_id), gw, width = 8, height = 9, dpi = 300)
}
    
    
    # ---------- Done ----------
    cat("\nSimulations and feature-importance exports written to 'inst/simulations/'.\n\n")
    