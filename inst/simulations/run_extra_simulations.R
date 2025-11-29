# ============================================================
# SCOUT-FDR — Extra FDR Curve Simulations
# ============================================================

library(scoutFDR)
library(quantreg)
library(ggplot2)
library(dplyr)
library(purrr)

# ---------- Utility functions ----------
bh <- function(p, alpha) {
  m <- length(p); o <- order(p); ps <- p[o]; k <- seq_len(m)
  thr <- k * alpha / m
  kmax <- max(which(ps <= thr), 0)
  if (kmax == 0) integer(0) else sort(o[ps <= ps[kmax]])
}

fdp <- function(rej, null_idx) if (length(rej)==0) 0 else sum(rej %in% null_idx)/length(rej)
exceed_prob <- function(fdp_vec, alpha) mean(fdp_vec > alpha, na.rm=TRUE)

# ---------- Scenario setup ----------
m_curve <- 5000
alphas_curve <- c(0.01,0.05,0.10)
B_reps_curve <- 100
K_curve <- 2

scenarios_extra <- data.frame(
  scen_id = paste0("L",1:15),
  n    = c(3,3,3, 5,5,5, 5,5,5, 8,8,8, 10,10,10),
  pi0  = c(0.2,0.4,0.7, 0.2,0.4,0.7, 0.2,0.4,0.7, 0.2,0.4,0.7, 0.2,0.4,0.7),
  power= c(0.20,0.20,0.20, 0.20,0.20,0.20, 0.35,0.35,0.35, 0.35,0.35,0.35, 0.50,0.50,0.50),
  stringsAsFactors=FALSE
)

set.seed(17)

# ---------- Load existing calibration data + CQR model ----------
cal_df  <- readRDS("inst/calibration/cal_df.rds")
cqr_fit <- readRDS("inst/calibration/cqr_fit.rds")
feature_names <- cqr_fit$feature_names

# ---------- SCOUT runner ----------
run_scout_once <- function(p, alpha, K, H, delta_ecdf, use_conformal,
                           cqr_fit, cal_df, feature_names) {
  scoutFDR::adaptive_fdr(
    p, alpha=alpha, K=K,
    delta_ecdf=delta_ecdf, H=H,
    use_conformal = use_conformal,
    cqr_fit = if (use_conformal) cqr_fit else NULL,
    cal_df   = if (use_conformal) cal_df else NULL,
    feature_names = if (use_conformal) feature_names else NULL
  )
}

# ---------- Evaluate scenario ----------
eval_scenario_alpha <- function(n, pi0, power, alpha,
                                conf_obj=NULL,
                                H, delta_ecdf,
                                tau_label,
                                alpha_cqr, max_weight,
                                K = K_curve, B = B_reps_curve) {

  out <- replicate(B, {
    p <- scoutFDR:::simulate_pvalues_ttest(
      m=m_curve, n=n, pi0=pi0, target_power=power
    )
    null_idx <- seq_len(round(m_curve*pi0))

    r_bh <- bh(p, alpha)
    A_med <- scoutFDR:::storey_median_lambda(p, c(0.6,0.7,0.8))
    r_st <- bh(p, alpha/A_med)

    res_safe <- run_scout_once(p, alpha, K, H, delta_ecdf,
                               FALSE, NULL, NULL, NULL)

    res_both <- if (!is.null(conf_obj)) scoutFDR::adaptive_fdr(
      p, alpha=alpha, K=K, delta_ecdf=delta_ecdf, H=H,
      mode="both", use_conformal=TRUE,
      cqr_fit=conf_obj$fit, cal_df=conf_obj$cal,
      feature_names=conf_obj$feats,
      ood_params=list(alpha_cqr=alpha_cqr, max_weight=max_weight)
    ) else NULL

    res_conf <- if (!is.null(conf_obj)) scoutFDR::adaptive_fdr(
      p, alpha=alpha, K=K, delta_ecdf=delta_ecdf, H=H,
      mode="conformal", use_conformal=TRUE,
      cqr_fit=conf_obj$fit, cal_df=conf_obj$cal,
      feature_names=conf_obj$feats,
      ood_params=list(alpha_cqr=alpha_cqr, max_weight=max_weight)
    ) else NULL

    c(
      FDP_BH   = fdp(r_bh, null_idx),
      FDP_ST   = fdp(r_st, null_idx),
      FDP_safe = fdp(res_safe$rejections, null_idx),
      FDP_both = if (is.null(res_both)) NA_real_
                 else fdp(res_both$rejections, null_idx),
      FDP_conf = if (is.null(res_conf)) NA_real_
                 else fdp(res_conf$rejections, null_idx)
    )
  }, simplify="matrix")

  data.frame(
    scen_id = NA_character_,
    n=n, pi0=pi0, power=power, alpha=alpha,
    method=c("BH","Storey","SCOUT_safe","SCOUT_both","SCOUT_confOnly"),
    mean_fdp=c(
      mean(out["FDP_BH",]),
      mean(out["FDP_ST",]),
      mean(out["FDP_safe",]),
      mean(out["FDP_both",], na.rm=TRUE),
      mean(out["FDP_conf",], na.rm=TRUE)
    ),
    exceed=c(
      exceed_prob(out["FDP_BH",], alpha),
      exceed_prob(out["FDP_ST",], alpha),
      exceed_prob(out["FDP_safe",], alpha),
      exceed_prob(out["FDP_both",], alpha),
      exceed_prob(out["FDP_conf",], alpha)
    ),
    H=paste(H, collapse=","),
    delta_ecdf=delta_ecdf,
    tau=tau_label,
    alpha_cqr=alpha_cqr,
    max_weight=max_weight,
    stringsAsFactors=FALSE
  )
}

# ---------- Block runner ----------
run_block <- function(scenarios, H, delta_ecdf, tau_label,
  alpha_cqr, max_weight, fname_prefix,
  cqr_fit, cal_df) {

  results <- list()
  conf_obj <- if (!is.null(cqr_fit)) 
  list(
    fit  = cqr_fit,
    cal  = cal_df,
    feats = cqr_fit$feature_names   # <-- always correct
  ) else NULL

  for (i in seq_len(nrow(scenarios))) {
    scen <- scenarios[i, ]
    for (a in alphas_curve) {
      res <- eval_scenario_alpha(
        n=scen$n, pi0=scen$pi0, power=scen$power, alpha=a,
        conf_obj=conf_obj,
        H=H, delta_ecdf=delta_ecdf,
        tau_label=tau_label, alpha_cqr=alpha_cqr, max_weight=max_weight
      )
      res$scen_id <- scen$scen_id
      results[[length(results)+1]] <- res
    }
  }

  df <- do.call(rbind, results)
  write.csv(df,
            sprintf("inst/simulations/%s.csv", fname_prefix),
            row.names=FALSE)

  p <- ggplot(df, aes(x=alpha, y=mean_fdp, color=method, group=method)) +
    geom_line(size=1) +
    geom_point(size=2) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    facet_wrap(~ scen_id, scales="free_y") +
    theme_bw() + theme(legend.position="bottom")

  ggsave(sprintf("inst/simulations/%s.png", fname_prefix),
         p, width=13, height=9, dpi=300)

  invisible(df)
}

# ---------- Run both settings ----------
message("Starting FDR simulations... Time:", Sys.time())
message("Starting fdr_curve_continuous_extra... Time:", Sys.time())
run_block(
  scenarios=scenarios_extra,
  H=c(0.10,0.20,0.30),
  delta_ecdf=0.02,
  tau_label=0.65,
  alpha_cqr=0.02,
  max_weight=4,
  fname_prefix="fdr_curve_continuous_extra",
  cqr_fit=cqr_fit,
  cal_df=cal_df
)
message("Starting fdr_curve_continuous_extra_conservative... Time:", Sys.time())
run_block(
  scenarios=scenarios_extra,
  H=c(0.02,0.05,0.10),
  delta_ecdf=0.005,
  tau_label=0.70,
  alpha_cqr=0.03,
  max_weight=3,
  fname_prefix="fdr_curve_continuous_extra_conservative",
  cqr_fit=cqr_fit,
  cal_df=cal_df
)
message("FDR simulations complete. Time:", Sys.time())