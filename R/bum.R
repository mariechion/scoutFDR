
# NOTE:
# ECDF features are the *default* for SCOUT-FDR.
# BUM is optional and only used when users explicitly request
# additional tail-shape features via combined_features().


#' Fit a Beta-Uniform Mixture (BUM) model
#'
#' @description
#' Fits a mixture of:
#' - a Uniform(0,1) component (nulls)
#' - K Beta(a_k, 1) components (alternatives)
#'
#' via an EM algorithm with damping and simple constraints.
#'
#' @param pvals Numeric vector of p-values.
#' @param K Number of Beta components (default 2).
#' @param a_min Minimum shape parameter.
#' @param a_max Maximum shape parameter.
#' @param max_iter Maximum EM iterations.
#' @param tol Convergence tolerance for log-likelihood.
#' @param damping Damping parameter in (0,1).
#'
#' @return A list containing:
#'   - converged: logical
#'   - w0: uniform mixing proportion
#'   - omega: mixture weights over Beta components
#'   - a: shape parameters
#'
#' @export
fit_bum <- function(
    pvals,
    K = 2,
    a_min = 1.05,
    a_max = 50,
    max_iter = 200,
    tol = 1e-6,
    damping = 0.3
) {

  init <- initialize_bum(pvals, K, a_min, a_max)
  if (is.null(init)) return(list(converged = FALSE))

  w0 <- init$w0
  omega <- init$omega
  a <- init$a

  loglik_history <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {

    f_curr <- compute_bum_density(pvals, w0, omega, a)
    f_curr <- pmax(f_curr, 1e-12)

    # responsibilities
    gamma_unif <- w0 / f_curr

    gamma_beta <- sapply(seq_len(K), function(k) {
      gk <- a[k] * (1 - pvals)^(a[k] - 1)
      (1 - w0) * omega[k] * gk / f_curr
    })

    # M-step
    w0_new <- (1 - damping) * mean(gamma_unif) + damping * w0

    sum_beta <- rowSums(gamma_beta)
    denom <- sum(sum_beta)
    if (denom <= 0) return(list(converged = FALSE))

    omega_new <- colSums(gamma_beta) / denom
    omega_new <- omega_new / sum(omega_new)

    # update a
    a_new <- numeric(K)
    log_q <- log(1 - pvals + 1e-12)

    for (k in seq_len(K)) {
      num <- sum(gamma_beta[, k])
      den <- -sum(gamma_beta[, k] * log_q)
      a_hat <- if (den > 0) num / den else a[k]
      a_hat <- min(max(a_hat, a_min), a_max)
      a_new[k] <- (1 - damping) * a_hat + damping * a[k]
    }

    # log-likelihood
    f_new <- compute_bum_density(pvals, w0_new, omega_new, a_new)
    f_new <- pmax(f_new, 1e-12)
    loglik <- sum(log(f_new))
    loglik_history[iter] <- loglik

    if (iter > 1) {
      delta <- loglik - loglik_history[iter - 1]

      if (abs(delta) < tol)
        return(list(converged = TRUE, w0 = w0_new, omega = omega_new, a = a_new))

      if (delta < -100 * tol)
        return(list(converged = FALSE))
    }

    w0 <- w0_new
    omega <- omega_new
    a <- a_new
  }

  list(converged = TRUE, w0 = w0, omega = omega, a = a)
}


initialize_bum <- function(pvals, K, a_min, a_max) {

  lambda_grid <- seq(0.5, 0.9, by = 0.1)
  pi0_est <- sapply(lambda_grid, function(l)
    sum(pvals > l) / ((1 - l) * length(pvals))
  )

  w0_init <- min(max(stats::median(pi0_est, na.rm = TRUE), 0.01), 0.99)

  p_lower <- pvals[pvals <= 0.3]
  a_init <- if (length(p_lower) >= 20)
    1 / mean(-log(1 - p_lower + 1e-12))
  else
    2

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

compute_bum_density <- function(pvals, w0, omega, a) {
  f <- rep(w0, length(pvals))
  for (k in seq_along(a)) {
    f <- f + (1 - w0) * omega[k] * a[k] * (1 - pvals)^(a[k] - 1)
  }
  f
}

extract_bum_features <- function(bum_fit, pvals, lambda = 0.5) {
  if (!isTRUE(bum_fit$converged)) {
    out <- rep(NA_real_, 8)
    names(out) <- paste0("bum_", seq_len(8))
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
  x5 <- mean(pvals > 0.9) / 0.1
  x6 <- max(omega)
  x7 <- K
  x8 <- w0 - x5

  out <- c(x1, x2, x3, x4, x5, x6, x7, x8)
  names(out) <- paste0("bum_", seq_len(8))
  out
}


#' Combined ECDF + optional BUM features
#'
#' @description
#' Extracts ECDF features and (optionally) appends BUM features.
#'
#' BUM is fitted only if it converges; otherwise, only ECDF features
#' are returned.
#'
#' @param pvals Numeric vector of p-values.
#' @param lambda Threshold for BUM feature x4.
#' @param compute_stability Whether to compute expensive stability features.
#'
#' @return A named list of features.
#'
#' @export
combined_features <- function(
    pvals,
    lambda = 0.5,
    compute_stability = TRUE
) {
  ecdf_feats <- features_ecdf(pvals, compute_stability = compute_stability)

  bf <- try(fit_bum(pvals, K = 2, max_iter = 300), silent = TRUE)

  if (!inherits(bf, "try-error") && isTRUE(bf$converged)) {
    bum_feats <- extract_bum_features(bf, pvals, lambda)
    return(as.list(c(ecdf_feats, bum_feats)))
  }

  # fallback: ECDF only
  as.list(ecdf_feats)
}
