

bh_reject <- function(pvals, level) {
  m <- length(pvals)
  o <- order(pvals)
  p_sorted <- pvals[o]
  thresh <- (seq_len(m) * level) / m

  k <- max(which(p_sorted <= thresh), 0)
  rej <- logical(m)

  if (k > 0) rej[o[seq_len(k)]] <- TRUE
  rej
}


#' Adaptive FDR control (Safe-U + optional Conformal)
#'
#' @description
#' This implements the SCOUT-FDR procedure:
#'
#' **Safe-U component:**  
#' - Computes an upper bound \eqn{U_k} on \eqn{\pi_0} using the one-sided
#'   DKW inequality on ECDF tail masses.
#'
#' **Optional conformal component:**  
#' - Cross-fitted quantile regression gives an additional, conditional
#'   upper bound on \eqn{\pi_0}, producing a sharper \eqn{U_k}.
#'
#' Final bound:  
#' - If `mode="both"` (default): \eqn{U_k = \min(U_k^{safe}, U_k^{conf})}
#'
#' @param pvals Numeric p-values.
#' @param alpha Target FDR level.
#' @param K Number of folds (default 2).
#' @param delta_ecdf Total miscoverage for Safe-U component.
#' @param H Tail thresholds (default: c(0.1, 0.2, 0.3)).
#' @param use_conformal Logical: whether to use CQR-based conformal bound.
#' @param cqr_fit Optional CQR model (from \code{train_cqr_pi0()}).
#' @param cal_df Calibration dataset used for CQR.
#' @param feature_names Feature columns used in CQR.
#' @param feature_fun Function producing feature vectors (default ECDF).
#' @param ood_params List with:
#'   - alpha_cqr: miscoverage level for conformal bound
#'   - max_weight: weight clipping for density ratio
#' @param mode One of "safe", "conformal", "both".
#'
#' @return A list with:
#'   - rejections: vector of rejected indices
#'   - U_fold: per-fold \eqn{\pi_0} upper bounds
#'   - safe_U: Safe-U bounds only
#'   - conf_U: conformal bounds only
#'   - fold_id: fold assignments
#'
#' @export
adaptive_fdr <- function(
    pvals,
    alpha        = 0.10,
    K            = 2,
    delta_ecdf   = 0.02,
    H            = c(0.10, 0.20, 0.30),
    use_conformal = FALSE,
    cqr_fit       = NULL,
    cal_df        = NULL,
    feature_names = NULL,
    feature_fun   = function(pp) as.list(features_ecdf(pp)),
    ood_params    = list(alpha_cqr = 0.05, max_weight = 4),
    mode          = c("both", "safe", "conformal")
) {

  mode <- match.arg(mode)
  if (missing(mode) && isTRUE(use_conformal)) mode <- "both"

  stopifnot(is.numeric(pvals), all(is.finite(pvals)), length(pvals) > 0)

  m <- length(pvals)
  if (K < 1 || K > m) stop("K must be in [1, m].")

  # Random fold assignment
  fold_id <- sample(rep(seq_len(K), length.out = m))

  # ---- Safe-U helper
  compute_safe <- function(train_p) {
    u <- try(safety_pi0_upper(train_p, H = H, delta_total = delta_ecdf, K = K),
             silent = TRUE)
    if (inherits(u, "try-error") || !is.finite(u)) u <- 1
    max(min(u, 1), 1e-8)
  }

  # ---- Prepare conformal features (if enabled)
  use_conf <- (mode %in% c("conformal", "both"))

  if (use_conf) {
    if (is.null(cqr_fit) || is.null(cal_df) || is.null(feature_names)) {
      if (mode == "conformal")
        stop("Conformal mode requested but cqr_fit/cal_df/feature_names missing.")
      warning("Conformal inputs missing; falling back to Safe-U only.")
      use_conf <- FALSE
    }
  }

  # Build target fold features
  X_target_all <- NULL
  if (use_conf) {

    feature_names <- intersect(feature_names, cqr_fit$feature_names %||% feature_names)

    X_target_all <- matrix(NA_real_, nrow = K, ncol = length(feature_names))
    colnames(X_target_all) <- feature_names

    for (k in seq_len(K)) {
      p_tar <- pvals[fold_id == k]
      feats_k <- try(feature_fun(p_tar), silent = TRUE)
      if (!inherits(feats_k, "try-error")) {
        nm <- intersect(names(feats_k), feature_names)
        X_target_all[k, nm] <- as.numeric(feats_k[nm])
      }
    }

    X_target_all <- as.data.frame(X_target_all)
    keep <- vapply(X_target_all, function(v) all(is.finite(v)), logical(1))
    X_target_all <- X_target_all[, keep, drop = FALSE]
    feature_names <- colnames(X_target_all)

    if (length(feature_names) == 0) {
      warning("No usable features for conformal CQR; falling back to Safe-U only.")
      use_conf <- FALSE
    }
  }

  # Conformal hyperparameters
  alpha_cqr <- ood_params$alpha_cqr %||% 0.05
  max_weight <- max(1, as.numeric(ood_params$max_weight %||% 4))
  weight_clip <- c(1 / max_weight, max_weight)

  # Storage
  safe_U <- numeric(K)
  conf_U <- rep(NA_real_, K)
  rejections <- integer(0)

  # ============================================================
  # Cross-fitted loop over folds
  # ============================================================
  for (k in seq_len(K)) {

    idx_tar <- which(fold_id == k)
    idx_trn <- which(fold_id != k)

    p_trn <- pvals[idx_trn]
    p_tar <- pvals[idx_tar]

    # ---- Safe-U bound
    safe_U[k] <- compute_safe(p_trn)

    # ---- Conformal bound (optional)
    if (use_conf) {

      x0 <- as.numeric(X_target_all[k, feature_names, drop = TRUE])
      names(x0) <- feature_names

      if (!all(is.finite(x0))) {
        conf_U[k] <- 1
      } else {
        u_conf <- try(
          cqr_upper_bound(
            cqr_fit      = cqr_fit,
            cal_df       = cal_df,
            feature_names = feature_names,
            x0           = x0,
            X_target_all = X_target_all,
            alpha        = alpha_cqr,
            weight_clip  = weight_clip
          ),
          silent = TRUE
        )
        if (inherits(u_conf, "try-error") || !is.finite(u_conf)) u_conf <- 1
        conf_U[k] <- max(min(u_conf, 1), 1e-8)
      }
    }

    # ---- Select which U to use
    U_k <- switch(
      mode,
      safe = safe_U[k],
      conformal = conf_U[k],
      both = {
        if (is.finite(conf_U[k])) min(safe_U[k], conf_U[k]) else safe_U[k]
      }
    )
    if (!is.finite(U_k) || U_k <= 0) U_k <- 1

    alpha_k <- alpha / U_k

    # ---- Reject in this fold
    rej_k <- bh_reject(p_tar, level = alpha_k)
    if (any(rej_k)) {
      rejections <- c(rejections, idx_tar[which(rej_k)])
    }
  }

  rejections <- sort(unique(rejections))

  list(
    rejections = rejections,
    U_fold     = switch(
      mode,
      safe = safe_U,
      conformal = conf_U,
      both = pmin(safe_U, ifelse(is.finite(conf_U), conf_U, 1))
    ),
    safe_U   = safe_U,
    conf_U   = conf_U,
    fold_id  = fold_id,
    params   = list(
      alpha = alpha, K = K, delta_ecdf = delta_ecdf,
      H = H, mode = mode, alpha_cqr = alpha_cqr,
      max_weight = max_weight
    )
  )
}


