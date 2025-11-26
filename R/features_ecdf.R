utils::globalVariables(c("have_MASS"))

#' @importFrom stats coef lm median setNames

# --- Small internal helpers --- #

.safe_quant <- function(x, prob, default = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default)
  stats::quantile(x, prob, na.rm = TRUE, names = FALSE)
}

.safe_stat <- function(x, fun, default = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default)
  fun(x)
}

# --- Main function --- #

#' ECDF-based feature extraction for p-values
#'
#' Computes a rich set of summary features of the empirical distribution
#' of p-values. These are intended as inputs to pi0 regression or other
#' meta-learning components.
#'
#' @param pvals Numeric vector of p-values.
#' @param H Tail bandwidths for Storey-like ratios.
#' @param ss_lambda0 Start of the upper-region for slope fitting.
#' @param alpha_pilot Pilot level used to define BH-derived features.
#' @param n_df Optional degrees-of-freedom metadata feature.
#'
#' @return Named numeric vector of features.
#' @export
features_ecdf <- function(pvals,
                          H = c(0.02, 0.05, 0.10),
                          ss_lambda0 = 0.7,
                          alpha_pilot = 0.10,
                          n_df = NULL) {

  stopifnot(is.numeric(pvals), length(pvals) > 0)

  ## --- Clean input --- ##
  pvals <- as.numeric(pvals)
  pvals <- pmin(pmax(pvals, .Machine$double.eps), 1 - .Machine$double.eps)
  pvals <- pvals[is.finite(pvals)]
  m <- length(pvals)

  ## --- Early return for very small m --- ##
  if (m < 20) {
    feature_names <- c(
      "x9","x10","x9_05","x9_10","x9_20","x9_curv",
      "A_lam60","A_lam70","A_lam80","A_lam90","A_slope","A_range",
      "F001","F010","F050","F100",
      "ER001","ER010","ER050","ER100",
      "R0","t0","Simes_min","KS","HC",
      "q01","q05","lnp_mean","lnp_med","lnp_p90",
      "m","df","uniq_ratio","boundary_mass"
    )
    out <- setNames(rep(0, length(feature_names)), feature_names)
    out["m"]  <- m
    out["df"] <- if (is.null(n_df)) NA_real_ else n_df
    return(out)
  }

  ## --- Main ECDF body --- ##
  Fhat <- stats::ecdf(pvals)
  ps <- sort(pvals)

  # Tail ratio function (Storey-like)
  Ah_fun <- function(h) {
    val <- (1 - Fhat(1 - h)) / h
    ifelse(is.finite(val), val, 0)
  }

  ## Tail-based ratios
  Ah_H  <- sapply(H, Ah_fun)
  x9    <- min(Ah_H)
  x9_05 <- Ah_fun(0.05)
  x9_10 <- Ah_fun(0.10)
  x9_20 <- Ah_fun(0.20)
  x9_curv <- x9_20 - 2*x9_10 + x9_05

  ## Slope estimation on upper portion
  t_grid <- seq(ss_lambda0, 1, length.out = 41)
  y_grid <- sapply(t_grid, Fhat)

  slope_main <- NA_real_
  if (isTRUE(getOption("have_MASS", FALSE)) || exists("have_MASS")) {
    slope_main <- tryCatch({
    fit <- suppressWarnings(
    MASS::rlm(y ~ t_grid, psi = MASS::psi.huber, k = 1.345, maxit = 50)
    )
    as.numeric(coef(fit)[2])
}, error = function(e) NA_real_)
  }
  if (!is.finite(slope_main)) {
    slope_main <- tryCatch(coef(lm(y_grid ~ t_grid))[2],
                           error = function(e) 0)
  }
  x10 <- slope_main

  ## A_lambda family
  A_lambda <- function(lam) {
    val <- (1 - Fhat(lam)) / (1 - lam)
    ifelse(is.finite(val), val, 0)
  }
  lamv   <- c(0.6, 0.7, 0.8, 0.9)
  A_vals <- sapply(lamv, A_lambda)
  A_lam60 <- A_vals[1]
  A_lam70 <- A_vals[2]
  A_lam80 <- A_vals[3]
  A_lam90 <- A_vals[4]
  A_slope <- A_lam90 - A_lam60
  A_range <- max(A_vals) - min(A_vals)

  ## Early rejections and enrichment
  gridL <- c(1e-3, 1e-2, 0.05, 0.10)
  Fg <- sapply(gridL, Fhat)
  ER <- Fg / gridL

  ## Pilot BH threshold
  k <- seq_len(m)
  thr0 <- (k * alpha_pilot) / m
  pass <- which(ps <= thr0)
  if (length(pass) > 0) {
    kmax <- max(pass)
    t0   <- ps[kmax]
    R0   <- kmax
  } else {
    t0 <- 0
    R0 <- 0
  }

  ## Simes
  Simes_min <- .safe_stat(m * ps / k, min, default = 1)

  ## KS
  ks_stat <- suppressWarnings(
    as.numeric(stats::ks.test(ps, "punif")$statistic)
  )
  ks_stat <- ifelse(is.finite(ks_stat), ks_stat, 0)

  ## Higher Criticism
  gridHC <- c(1e-3, 1e-2, 0.05, 0.10)
  HC <- max(sapply(gridHC, function(u) {
    sqrt(m) * abs(Fhat(u) - u) / sqrt(u*(1-u))
  }))

  ## Log-p transforms
  lp <- -log(ps)
  lnp_mean <- .safe_stat(lp, mean)
  lnp_med  <- .safe_stat(lp, median)
  lnp_p90  <- .safe_quant(lp, 0.9)

  uniq_ratio    <- length(unique(ps)) / m
  boundary_mass <- mean(ps <= .Machine$double.eps | ps >= 1 - .Machine$double.eps)

  ## --- Return all features --- ##

  c(
    x9 = x9,
    x10 = x10,
    x9_05 = x9_05, x9_10 = x9_10, x9_20 = x9_20, x9_curv = x9_curv,
    A_lam60 = A_lam60, A_lam70 = A_lam70,
    A_lam80 = A_lam80, A_lam90 = A_lam90,
    A_slope = A_slope, A_range = A_range,
    F001 = Fg[1], F010 = Fg[2], F050 = Fg[3], F100 = Fg[4],
    ER001 = ER[1], ER010 = ER[2], ER050 = ER[3], ER100 = ER[4],
    R0 = R0, t0 = t0, Simes_min = Simes_min,
    KS = ks_stat, HC = HC,
    q01 = .safe_quant(ps, 0.01),
    q05 = .safe_quant(ps, 0.05),
    lnp_mean = lnp_mean, lnp_med = lnp_med, lnp_p90 = lnp_p90,
    m = m,
    df = if (is.null(n_df)) NA_real_ else n_df,
    uniq_ratio = uniq_ratio,
    boundary_mass = boundary_mass
  )
}
