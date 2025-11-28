#' ECDF-based feature extraction for p-values
#'
#' This function computes a rich collection of features derived from the
#' empirical CDF (ECDF) of a vector of p-values. These features are used
#' by \code{scoutFDR} to model local π₀ behavior through calibrated
#' regression (CQR).
#'
#' The feature set includes:
#'   * Tail mass ratios (A_h) for several right-tail intervals.
#'   * Storey-type estimates across multiple λ values.
#'   * Slopes and curvature diagnostics for ECDF near 1.
#'   * Early-tail enrichment features (F₀₀₁, F₀₁₀, F₀₅₀, etc.).
#'   * Simes, KS, and Higher Criticism (HC) summaries.
#'   * Log-p summaries, quantiles, uniqueness ratios.
#'   * Right-tail densities, enrichment ratios, pilot BH stats.
#'   * Optional stability diagnostics via subsampling.
#'
#' @param pvals Numeric vector of p-values between 0 and 1.
#' @param H Vector of right-tail widths for computing A_h ratios.
#' @param ss_lambda0 Start of the tail slope regression segment (default 0.7).
#' @param alpha_pilot Pilot BH level for early rejection statistics.
#' @param n_df Optional degrees-of-freedom metadata.
#' @param compute_stability Whether to compute subsample-based stability features.
#'
#' @return A named numeric vector of ECDF-based features.
#' 
#' @importFrom graphics hist
#' @export
features_ecdf <- function(pvals,
                          H = c(0.02, 0.05, 0.10),
                          ss_lambda0 = 0.7,
                          alpha_pilot = 0.10,
                          n_df = NULL,
                          compute_stability = TRUE) {

  stopifnot(is.numeric(pvals), length(pvals) > 0)

  # Clean p-values
  pvals <- pmax(pmin(as.numeric(pvals), 1), 0)
  m <- length(pvals)
  Fhat <- stats::ecdf(pvals)

  # Sort once
  o <- order(pvals)
  ps <- pvals[o]
  k  <- seq_len(m)

  # ---- Tail mass A_h ratios ----
  Ah_fun <- function(h) (1 - Fhat(1 - h)) / h
  Ah_H   <- sapply(H, Ah_fun)
  x9     <- min(Ah_H)

  # ---- Right-tail slope via robust regression (or fallback) ----
  t_grid <- seq(ss_lambda0, 1, length.out = 41)
  y_grid <- sapply(t_grid, Fhat)

  slope_main <- NA_real_
  if (exists("have_MASS") && isTRUE(have_MASS)) {
    slope_main <- tryCatch({
      fit <- MASS::rlm(y_grid ~ t_grid, psi = MASS::psi.huber)
      as.numeric(stats::coef(fit)[2])
    }, error = function(e) NA_real_)
  }
  if (!is.finite(slope_main)) {
    slope_main <- tryCatch({
      as.numeric(stats::coef(stats::lm(y_grid ~ t_grid))[2])
    }, error = function(e) NA_real_)
  }
  x10 <- slope_main

  # ---- Curvature of A_h over standard h ----
  x9_05 <- Ah_fun(0.05)
  x9_10 <- Ah_fun(0.10)
  x9_20 <- Ah_fun(0.20)
  x9_curv <- x9_20 - 2 * x9_10 + x9_05

  # ---- Storey-type A_lambda values ----
  A_lambda <- function(l) (1 - Fhat(l)) / (1 - l)
  lamv     <- c(0.6, 0.7, 0.8, 0.9)
  A_vals   <- sapply(lamv, A_lambda)
  names(A_vals) <- paste0("A_lam", c("60","70","80","90"))

  A_slope <- (A_vals["A_lam90"] - A_vals["A_lam60"]) / (0.9 - 0.6)
  A_range <- max(A_vals) - min(A_vals)

  # ---- Early-tail CDF + enrichment ----
  gridL <- c(1e-3, 1e-2, 0.05, 0.10)
  Fg    <- sapply(gridL, Fhat)
  names(Fg) <- c("F001","F010","F050","F100")
  ER    <- Fg / gridL
  names(ER) <- c("ER001","ER010","ER050","ER100")

  # ---- Pilot BH at alpha_pilot ----
  thr0  <- (k * alpha_pilot) / m
  pass0 <- which(ps <= thr0)
  if (length(pass0) > 0) {
    R0 <- max(pass0)
    t0 <- ps[R0]
  } else {
    R0 <- 0L
    t0 <- NA_real_
  }

  # ---- Simes minimum ----
  Simes_min <- min(m * ps / k)

  # ---- KS + Higher Criticism ----
  ks_stat <- suppressWarnings(
    as.numeric(stats::ks.test(ps, "punif")$statistic)
  )
  HC <- {
    gridHC <- c(1e-3, 1e-2, 0.05, 0.10)
    vals   <- sapply(gridHC, function(u) {
      sqrt(m) * abs(Fhat(u) - u) / sqrt(u * (1 - u))
    })
    max(vals)
  }

  # ---- Log-p and quantiles ----
  lp <- -log(pmax(ps, .Machine$double.eps))
  lnp_mean <- mean(lp)
  lnp_med  <- stats::median(lp)
  lnp_p90  <- stats::quantile(lp, 0.90, names = FALSE)

  q01 <- ps[max(1L, ceiling(0.01 * m))]
  q05 <- ps[max(1L, ceiling(0.05 * m))]

  # ---- Data quality ----
  uniq_ratio <- length(unique(ps)) / m
  boundary_mass <- mean(ps == 0 | ps == 1)
  df <- if (is.null(n_df)) NA_real_ else as.numeric(n_df)

  # ---- Tail proxies (simple) ----
  left_tail_01 <- mean(ps <= 0.01)
  left_tail_05 <- mean(ps <= 0.05)
  right_tail_95 <- mean(ps > 0.95)
  enrichment_ratio <- {
    enrich_early <- if (left_tail_05 > 0) left_tail_05 / 0.05 else 0
    enrich_late  <- if (right_tail_95 > 0) right_tail_95 / 0.05 else 1e-8
    enrich_early / enrich_late
  }

  # ---- Frequency diagnostics ----
  p_tab <- table(ps)
  max_freq <- max(p_tab) / m
  tail_unique_ratio <- if (sum(ps > 0.9) > 0) {
    length(unique(ps[ps > 0.9])) / sum(ps > 0.9)
  } else 1

  # ---- Optional subsampling stability ----
  A_07_subsample_sd <- NA_real_
  if (compute_stability && m >= 100) {
    A_07_subsample_sd <- tryCatch({
      n_sub <- max(20L, floor(0.8 * m))
      A_subs <- replicate(20, {
        idx <- sample.int(m, n_sub, replace = FALSE)
        Ah_fun_small <- function(h) {
          (1 - stats::ecdf(ps[idx])(1 - h)) / h
        }
        Ah_fun_small(0.07)
      })
      stats::sd(A_subs)
    }, error = function(e) NA_real_)
  }

  # ---- Right-tail bins ----
  t80_90  <- mean(ps > 0.80 & ps <= 0.90)
  t90_95  <- mean(ps > 0.90 & ps <= 0.95)
  t95_100 <- mean(ps > 0.95)

  d80_90 <- t80_90 / 0.10
  d90_95 <- t90_95 / 0.05
  d95_100 <- t95_100 / 0.05

  # ---- Clopper-Pearson upper bounds ----
  cp_upper <- function(x, n, h, alpha = 0.05) {
    stats::qbeta(1 - alpha, x + 1, n - x) / h
  }
  CP05 <- cp_upper(sum(ps > 0.95), m, 0.05)
  CP10 <- cp_upper(sum(ps > 0.90), m, 0.10)
  CP20 <- cp_upper(sum(ps > 0.80), m, 0.20)

  # ---- Center / flatness ----
  mean_p <- mean(ps)
  median_p <- stats::median(ps)
  br <- seq(0, 1, length.out = 101)
  histc <- hist(ps, breaks = br, plot = FALSE)$counts / m
  F_emp <- cumsum(histc)
  F_uni <- seq(0.01, 1, by = 0.01)
  CvM   <- mean((F_emp - F_uni)^2)

  # ---- Additional right-tail Storey ratios ----
  A_lam95 <- (1 - Fhat(0.95)) / 0.05
  A_lam98 <- (1 - Fhat(0.98)) / 0.02

  # ---- Return all features ----
  c(
    x9 = x9, x10 = x10,
    x9_05 = x9_05, x9_10 = x9_10, x9_20 = x9_20, x9_curv = x9_curv,
    A_lam60 = A_vals["A_lam60"], A_lam70 = A_vals["A_lam70"],
    A_lam80 = A_vals["A_lam80"], A_lam90 = A_vals["A_lam90"],
    A_slope = A_slope, A_range = A_range,
    F001 = Fg["F001"], F010 = Fg["F010"], F050 = Fg["F050"], F100 = Fg["F100"],
    ER001 = ER["ER001"], ER010 = ER["ER010"], ER050 = ER["ER050"], ER100 = ER["ER100"],
    R0 = as.numeric(R0), t0 = as.numeric(t0), Simes_min = as.numeric(Simes_min),
    KS = ks_stat, HC = HC,
    q01 = as.numeric(q01), q05 = as.numeric(q05),
    lnp_mean = lnp_mean, lnp_med = lnp_med, lnp_p90 = lnp_p90,
    m = m, df = df,
    uniq_ratio = uniq_ratio, boundary_mass = boundary_mass,
    left_tail_01 = left_tail_01, left_tail_05 = left_tail_05,
    right_tail_95 = right_tail_95, enrichment_ratio = enrichment_ratio,
    max_freq = max_freq, tail_unique_ratio = tail_unique_ratio,
    A_07_subsample_sd = A_07_subsample_sd,
    tail_80_90 = t80_90, tail_90_95 = t90_95, tail_95_100 = t95_100,
    dens_80_90 = d80_90, dens_90_95 = d90_95, dens_95_100 = d95_100,
    CP05 = CP05, CP10 = CP10, CP20 = CP20,
    mean_p = mean_p, median_p = median_p, CvM = CvM,
    A_lam95 = A_lam95, A_lam98 = A_lam98
  )
}
