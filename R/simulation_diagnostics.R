#' Simulate labeled t-test experiment (null vs alternative)
#'
#' Generates p-values from the same two-sample t-test mixture model used
#' in [simulate_pvalues_ttest()], but additionally returns the indices
#' of null and alternative hypotheses. This is useful for computing
#' realised FDP/FDR in simulations.
#'
#' @inheritParams simulate_pvalues_ttest
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{p}: numeric vector of p-values of length \code{m},
#'   \item \code{null_idx}: integer indices of true null hypotheses,
#'   \item \code{alt_idx}: integer indices of true alternatives.
#' }
#' @export
simulate_ttest_labeled <- function(m, n, pi0, target_power, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  m0 <- round(m * pi0)
  m1 <- m - m0

  p <- simulate_pvalues_ttest(m = m, n = n, pi0 = pi0, target_power = target_power)

  null_idx <- if (m0 > 0) seq_len(m0) else integer(0)
  alt_idx  <- if (m1 > 0) (m0 + 1L):m else integer(0)

  list(p = p, null_idx = null_idx, alt_idx = alt_idx)
}

# Internal helper: false discovery proportion
fdp <- function(rej, null_idx) {
  if (length(rej) == 0) {
    return(0)
  }
  sum(rej %in% null_idx) / length(rej)
}

# Internal helper: benchmark a single replicate
bench_one <- function(cal_df, cqr_fit, feature_names,
                      m = 2000, n = 5, pi0 = 0.2, target_power = 0.35,
                      alpha = 0.10, seed = NULL,
                      K = 8,
                      base_pi0 = make_base_pi0("ecdf_storey_median"),
                      delta_ecdf = 0.02,
                      H = c(0.05, 0.10, 0.20)) {
  sim <- simulate_ttest_labeled(m = m, n = n, pi0 = pi0,
                                target_power = target_power, seed = seed)
  p <- sim$p
  null_idx <- sim$null_idx

  # Safe-U only
  r1 <- adaptive_fdr_safeU(
    p = p,
    alpha = alpha,
    K = K,
    base_pi0 = base_pi0,
    delta_ecdf = delta_ecdf,
    H = H,
    use_conformal = FALSE
  )

  # Safe-U + conformal (CQR)
  r2 <- adaptive_fdr_safeU(
    p = p,
    alpha = alpha,
    K = K,
    base_pi0 = base_pi0,
    delta_ecdf = delta_ecdf,
    H = H,
    use_conformal = TRUE,
    cqr_fit = cqr_fit,
    cal_df = cal_df,
    feature_names = feature_names
  )

  c(
    Umin_safe = min(r1$U_fold, na.rm = TRUE),
    Umin_conf = min(r2$U_fold, na.rm = TRUE),
    disc_safe = length(r1$rejections),
    disc_conf = length(r2$rejections),
    fdp_safe = fdp(r1$rejections, null_idx),
    fdp_conf = fdp(r2$rejections, null_idx)
  )
}

#' Run Monte Carlo benchmark for Safe-U and conformal Safe-U
#'
#' Runs a Monte Carlo simulation comparing the baseline Safe-U procedure
#' to its conformal refinement using a pre-trained CQR model for \code{pi0}.
#' Internally, this function:
#'
#' \enumerate{
#'   \item Generates calibration data using [generate_calibration_data()],
#'   \item Trains a CQR model via [train_cqr_pi0()],
#'   \item Repeats \code{B} independent simulation runs and for each,
#'         computes summary statistics (minimum per-fold \code{U},
#'         number of discoveries, and FDP for each method).
#' }
#'
#' @param B Integer, number of Monte Carlo replicates.
#' @param m Number of hypotheses per replicate.
#' @param n Sample size per group for each t-test.
#' @param pi0 Proportion of nulls in the simulator.
#' @param target_power Target power used to calibrate non-null effects.
#' @param alpha Target FDR level.
#'
#' @return A data frame with one row per replicate and columns:
#'   \code{Umin_safe}, \code{Umin_conf}, \code{disc_safe},
#'   \code{disc_conf}, \code{fdp_safe}, \code{fdp_conf}. The function
#'   also prints summary statistics to the console.
#'
#' @examples
#' \dontrun{
#' # Small benchmark (B = 20) just to check everything runs
#' tab <- run_benchmark(B = 20, m = 1000, n = 5, pi0 = 0.2,
#'                      target_power = 0.35, alpha = 0.10)
#' head(tab)
#' }
#'
#' @export
run_benchmark <- function(B = 200,
                          m = 2000, n = 5, pi0 = 0.2, target_power = 0.35,
                          alpha = 0.10) {
  if (!isTRUE(have_qr)) {
    stop("quantreg not installed; cannot run conformal benchmark.")
  }

  # 1) Train calibration model once
  cal_df <- generate_calibration_data(
    n_sim = 1000,
    base_pi0 = make_base_pi0("ecdf_storey_median"),
    feature_fun = function(p) as.list(features_ecdf(p))
  )

  feature_names_ecdfA <- intersect(
    c(
      "A",
      "x9", "x10", "x9_05", "x9_10", "x9_20", "x9_curv",
      "A_lam60", "A_lam70", "A_lam80", "A_lam90", "A_slope", "A_range",
      "F050", "ER050",
      "R0", "t0", "Simes_min",
      "KS", "HC",
      "lnp_med", "lnp_mean", "lnp_p90",
      "m", "uniq_ratio", "boundary_mass"
    ),
    colnames(cal_df)
  )

  cqr_fit <- train_cqr_pi0(
    cal_df = cal_df,
    feature_names = feature_names_ecdfA,
    tau = 0.8,
    model = "qr"
  )

  # 2) Monte Carlo loop
  set.seed(1)
  tab <- t(sapply(seq_len(B), function(b) {
    bench_one(
      cal_df = cal_df,
      cqr_fit = cqr_fit,
      feature_names = feature_names_ecdfA,
      m = m,
      n = n,
      pi0 = pi0,
      target_power = target_power,
      alpha = alpha,
      seed = 1000 + b
    )
  }))

  tab <- as.data.frame(tab)

  # 3) Print summaries
  cat("Mean Umin (safe, conf):",
      mean(tab$Umin_safe, na.rm = TRUE),
      mean(tab$Umin_conf, na.rm = TRUE), "\n")
  cat("Mean discoveries (safe, conf):",
      mean(tab$disc_safe, na.rm = TRUE),
      mean(tab$disc_conf, na.rm = TRUE), "\n")
  cat("Coverage P(Umin_conf >= pi0):",
      mean(tab$Umin_conf >= pi0, na.rm = TRUE), "\n")
  cat("Mean FDR (safe, conf):",
      mean(tab$fdp_safe, na.rm = TRUE),
      mean(tab$fdp_conf, na.rm = TRUE), "\n")

  tab
}

#' One-shot conformal example with Safe-U comparison
#'
#' Runs a single simulated experiment (similar to the examples used in
#' the development script) to compare the Safe-U-only procedure to its
#' conformal refinement, printing basic diagnostics to the console.
#'
#' This helper is mainly for interactive exploration, not for use in
#' production analyses.
#'
#' @return Invisibly, a list with components \code{res_safeU} and
#'   \code{res_conf}, the outputs of [adaptive_fdr_safeU()] without and
#'   with conformal tightening respectively.
#'
#' @examples
#' \dontrun{
#' out <- run_conformal_with_diagnostics()
#' }
#'
#' @export
run_conformal_with_diagnostics <- function() {
  set.seed(123)
  p <- simulate_pvalues_ttest(
    m = 2000,
    n = 5,
    pi0 = 0.2,
    target_power = 0.35
  )

  # Safe-U only
  res_safeU <- adaptive_fdr_safeU(
    p = p,
    alpha = 0.10,
    K = 8,
    base_pi0 = make_base_pi0("ecdf_storey_median"),
    delta_ecdf = 0.02,
    H = c(0.05, 0.10, 0.20),
    use_conformal = FALSE
  )

  cat("Safe-U only: discoveries =",
      length(res_safeU$rejections),
      " | mean(U_fold) =",
      round(mean(res_safeU$U_fold, na.rm = TRUE), 3), "\n")

  if (!isTRUE(have_qr)) {
    cat("quantreg not installed; skipping conformal example.\n")
    return(invisible(list(res_safeU = res_safeU, res_conf = NULL)))
  }

  # Calibration for conformal run
  cal_df <- generate_calibration_data(
    n_sim = 1000,
    base_pi0 = make_base_pi0("ecdf_storey_median"),
    feature_fun = function(p) as.list(features_ecdf(p))
  )

  feature_names_ecdfA <- intersect(
    c(
      "A",
      "x9", "x10", "x9_05", "x9_10", "x9_20", "x9_curv",
      "A_lam60", "A_lam70", "A_lam80", "A_lam90", "A_slope", "A_range",
      "F050", "ER050",
      "R0", "t0", "Simes_min",
      "KS", "HC",
      "lnp_med", "lnp_mean", "lnp_p90",
      "m", "uniq_ratio", "boundary_mass"
    ),
    colnames(cal_df)
  )

  cqr_fit <- train_cqr_pi0(
    cal_df = cal_df,
    feature_names = feature_names_ecdfA,
    tau = 0.8,
    model = "qr"
  )

  res_conf <- adaptive_fdr_safeU(
    p = p,
    alpha = 0.10,
    K = 8,
    base_pi0 = make_base_pi0("ecdf_storey_median"),
    delta_ecdf = 0.02,
    H = c(0.05, 0.10, 0.20),
    use_conformal = TRUE,
    cqr_fit = cqr_fit,
    cal_df = cal_df,
    feature_names = feature_names_ecdfA,
    feature_fun = function(pp) as.list(features_ecdf(pp)),
    ood_params = list(max_weight = 2, maha_alpha = 0.01, alpha_cqr = 0.10)
  )

  cat("Safe-U + conformal: discoveries =",
      length(res_conf$rejections),
      " | mean(U_fold) =",
      round(mean(res_conf$U_fold, na.rm = TRUE), 3), "\n")

  invisible(list(res_safeU = res_safeU, res_conf = res_conf))
}
