utils::globalVariables(c("have_MASS"))

#' ECDF-based feature extraction for p-values
#'
#' Computes a rich set of summary features of the empirical distribution
#' of p-values. These are intended as inputs to pi0 regression or other
#' meta-learning components.
#'
#' @param p Numeric vector of p-values.
#' @param H Tail bandwidths used for some Storey-type ratios.
#' @param ss_lambda0 Starting point of the linear region for trend fitting.
#' @param alpha_pilot Pilot FDR level used to define BH-type features.
#' @param n_df Optional degrees-of-freedom feature (e.g. from a t-test).
#'
#' @return A named numeric vector of features.
#' @export
features_ecdf <- function(p,
                          H = c(0.02, 0.05, 0.10),
                          ss_lambda0 = 0.7,
                          alpha_pilot = 0.10,
                          n_df = NULL) {
  stopifnot(is.numeric(p), length(p) > 0)
  m <- length(p)
  p <- pmax(pmin(as.numeric(p), 1), 0)
  Fhat <- stats::ecdf(p)
  o <- order(p); ps <- p[o]

  # Tail-based ratios A_h = (1 - Fhat(1 - h)) / h
  Ah_fun <- function(h) (1 - Fhat(1 - h)) / h
  Ah_H <- sapply(H, Ah_fun)
  x9 <- min(Ah_H)

  # Robust slope of Fhat in the upper region [ss_lambda0, 1]
  t_grid <- seq(ss_lambda0, 1, length.out = 41)
  y <- sapply(t_grid, Fhat)
  slope_main <- NA_real_
  if (exists("have_MASS") && isTRUE(have_MASS)) {
    slope_main <- tryCatch({
      fit <- MASS::rlm(y ~ t_grid, psi = MASS::psi.huber, k = 1.345, maxit = 50)
      as.numeric(stats::coef(fit)[2])
    }, warning = function(w) NA_real_, error = function(e) NA_real_)
  }
  if (!is.finite(slope_main)) {
    slope_main <- tryCatch({
      as.numeric(stats::coef(stats::lm(y ~ t_grid))[2])
    }, error = function(e) NA_real_)
  }
  x10 <- slope_main

  # Curvature of A_h at specific h values
  x9_05 <- Ah_fun(0.05)
  x9_10 <- Ah_fun(0.10)
  x9_20 <- Ah_fun(0.20)
  Ah_vec <- c(x9_05, x9_10, x9_20)
  x9_curv <- Ah_vec[3] - 2 * Ah_vec[2] + Ah_vec[1]

  # A_lambda type Storey estimators at a grid of lambdas
  A_lambda <- function(lam) (1 - Fhat(lam)) / (1 - lam)
  lamv <- c(0.6, 0.7, 0.8, 0.9)
  A_vals <- sapply(lamv, A_lambda)
  A_lam60 <- A_vals[1]; A_lam70 <- A_vals[2]; A_lam80 <- A_vals[3]; A_lam90 <- A_vals[4]
  A_slope <- (A_lam90 - A_lam60) / (0.9 - 0.6)
  A_range <- max(A_vals) - min(A_vals)

  # Early rejection and enrichment ratios
  gridL <- c(1e-3, 1e-2, 0.05, 0.10)
  Fg <- sapply(gridL, Fhat)
  names(Fg) <- c("F001", "F010", "F050", "F100")
  ER <- Fg / gridL
  names(ER) <- c("ER001", "ER010", "ER050", "ER100")

  # BH pilot threshold and Simes statistic
  k <- seq_len(m)
  thr0 <- (k * alpha_pilot) / m
  pass <- which(ps <= thr0)
  if (length(pass) > 0) {
    kmax <- max(pass); t0 <- ps[kmax]; R0 <- kmax
  } else {
    t0 <- NA_real_; R0 <- 0L
  }
  Simes_min <- min(m * ps / k)

  # Goodness-of-fit statistics
  ks_stat <- suppressWarnings(as.numeric(stats::ks.test(ps, "punif")$statistic))
  gridHC <- c(1e-3, 1e-2, 0.05, 0.10)
  HC <- {
    vals <- sapply(gridHC, function(u) {
      sqrt(m) * abs(Fhat(u) - u) / sqrt(u * (1 - u))
    })
    max(vals)
  }

  # Quantiles and log-p features
  q01 <- ps[max(1L, ceiling(0.01 * m))]
  q05 <- ps[max(1L, ceiling(0.05 * m))]
  lp <- -log(pmax(ps, .Machine$double.eps))
  lnp_mean <- mean(lp)
  lnp_med  <- stats::median(lp)
  lnp_p90  <- stats::quantile(lp, 0.90, names = FALSE)

  # Misc structural features
  df <- if (is.null(n_df)) NA_real_ else as.numeric(n_df)
  uniq_ratio <- length(unique(ps)) / m
  boundary_mass <- mean(ps == 0 | ps == 1)

  c(
    x9 = x9, x10 = x10,
    x9_05 = x9_05, x9_10 = x9_10, x9_20 = x9_20, x9_curv = x9_curv,
    A_lam60 = A_lam60, A_lam70 = A_lam70, A_lam80 = A_lam80, A_lam90 = A_lam90,
    A_slope = A_slope, A_range = A_range,
    F001 = Fg[["F001"]], F010 = Fg[["F010"]], F050 = Fg[["F050"]], F100 = Fg[["F100"]],
    ER001 = ER[["ER001"]], ER010 = ER[["ER010"]], ER050 = ER[["ER050"]], ER100 = ER[["ER100"]],
    R0 = as.numeric(R0), t0 = as.numeric(t0), Simes_min = as.numeric(Simes_min),
    KS = ks_stat, HC = as.numeric(HC),
    q01 = as.numeric(q01), q05 = as.numeric(q05),
    lnp_mean = as.numeric(lnp_mean), lnp_med = as.numeric(lnp_med), lnp_p90 = as.numeric(lnp_p90),
    m = as.numeric(m), df = df, uniq_ratio = as.numeric(uniq_ratio),
    boundary_mass = as.numeric(boundary_mass)
  )
}

# =========================
# BUM (Beta-Uniform Mixture) fitting (internal)
# =========================

fit_bum <- function(p_values, K = 2, a_min = 1.05, a_max = 50,
                    max_iter = 200, tol = 1e-6, damping = 0.3) {
  init <- initialize_bum(p_values, K, a_min, a_max)
  if (is.null(init)) return(list(converged = FALSE))
  w0 <- init$w0; omega <- init$omega; a <- init$a
  loglik_history <- numeric(max_iter)

  for (iter in 1:max_iter) {
    f_curr <- compute_bum_density(p_values, w0, omega, a)
    f_curr <- pmax(f_curr, 1e-12)

    gamma_unif <- w0 / f_curr
    gamma_beta <- sapply(seq_len(K), function(k) {
      gk <- a[k] * (1 - p_values)^(a[k] - 1)
      (1 - w0) * omega[k] * gk / f_curr
    })

    w0_new <- (1 - damping) * mean(gamma_unif) + damping * w0

    sum_beta <- rowSums(gamma_beta)
    denom <- sum(sum_beta)
    if (denom <= 0) return(list(converged = FALSE))
    omega_new <- colSums(gamma_beta) / denom
    omega_new <- omega_new / sum(omega_new)

    a_new <- numeric(K)
    log_q <- log(1 - p_values + 1e-12)
    for (k in 1:K) {
      num <- sum(gamma_beta[, k])
      den <- -sum(gamma_beta[, k] * log_q)
      a_hat <- if (den > 0) num / den else a[k]
      a_hat <- min(max(a_hat, a_min), a_max)
      a_new[k] <- (1 - damping) * a_hat + damping * a[k]
    }

    f_new <- compute_bum_density(p_values, w0_new, omega_new, a_new)
    f_new <- pmax(f_new, 1e-12)
    loglik <- sum(log(f_new))
    loglik_history[iter] <- loglik

    if (iter > 1) {
      delta <- loglik_history[iter] - loglik_history[iter - 1]
      if (abs(delta) < tol) {
        return(list(converged = TRUE, w0 = w0_new, omega = omega_new, a = a_new))
      }
      if (delta < -100 * tol) {
        return(list(converged = FALSE))
      }
    }

    w0 <- w0_new
    omega <- omega_new
    a <- a_new
  }
  list(converged = TRUE, w0 = w0, omega = omega, a = a)
}

initialize_bum <- function(p_values, K, a_min, a_max) {
  lambda_grid <- seq(0.5, 0.9, by = 0.1)
  pi0_est <- sapply(lambda_grid, function(l) {
    sum(p_values > l) / ((1 - l) * length(p_values))
  })
  w0_init <- min(max(stats::median(pi0_est, na.rm = TRUE), 0.01), 0.99)

  p_lower <- p_values[p_values <= 0.3]
  a_init <- if (length(p_lower) >= 20) {
    1 / mean(-log(1 - p_lower + 1e-12))
  } else {
    2
  }
  a_init <- min(max(a_init, a_min + 0.1), a_max - 0.1)

  if (K == 1) {
    omega_init <- 1
    a_vec <- a_init
  } else {
    a_vec <- c(
      max(a_init, a_min + 0.1),
      min(max(a_init * 1.5, a_min + 0.2), a_max - 0.2)
    )
    omega_init <- c(0.7, 0.3)
  }
  list(w0 = w0_init, omega = omega_init, a = a_vec)
}

compute_bum_density <- function(p_values, w0, omega, a) {
  f <- rep(w0, length(p_values))
  for (k in seq_along(a)) {
    f <- f + (1 - w0) * omega[k] * a[k] * (1 - p_values)^(a[k] - 1)
  }
  f
}

extract_bum_features <- function(bum_fit, p_values, lambda = 0.5) {
  if (!isTRUE(bum_fit$converged)) {
    out <- rep(NA_real_, 8)
    names(out) <- paste0("x", 1:8)
    return(out)
  }
  w0 <- bum_fit$w0
  omega <- bum_fit$omega
  a <- bum_fit$a
  K <- length(a)

  x1 <- w0
  x2 <- sum(omega * a)
  x3 <- sum(omega * a^2) - x2^2
  x4 <- sum(omega * (1 - lambda)^a)
  x5 <- mean(p_values > 0.9) / 0.1
  x6 <- max(omega)
  x7 <- K
  x8 <- w0 - x5

  out <- c(x1, x2, x3, x4, x5, x6, x7, x8)
  names(out) <- paste0("x", 1:8)
  out
}

# =========================
# Combined features (ECDF + optional BUM)
# =========================

#' Combined ECDF and BUM-based features
#'
#' Computes ECDF-based features and, when the BUM fit converges, additional
#' mixture-model-based features. If the BUM fit fails, only ECDF features
#' are returned.
#'
#' @param p Numeric vector of p-values.
#' @param lambda Threshold used in some BUM-based features.
#'
#' @return A named list of features.
#' @export
combined_features <- function(p, lambda = 0.5) {
  ecdf_feats <- features_ecdf(p)
  bf <- try(fit_bum(p, K = 2, max_iter = 300), silent = TRUE)
  if (!inherits(bf, "try-error") && isTRUE(bf$converged)) {
    bum_feats <- extract_bum_features(bf, p, lambda)
    as.list(c(ecdf_feats, bum_feats))
  } else {
    as.list(ecdf_feats)
  }
}
