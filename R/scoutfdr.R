utils::globalVariables(c("have_cp4p"))

#' Base pi0 estimators and Safe-U ECDF upper bound
#'
#' Storey-type ECDF estimators and cp4p-based estimators that act as
#' base features for the adaptive procedures.
#'
#' @param type Character string specifying the estimator type.
#' @param Lambda Optional vector of lambda thresholds for Storey estimators.
#' @param eps Numerical stability parameter.
#' @param nbins,pz,censor,cp4p_fallback Arguments passed to cp4p-based estimators.
#'
#' @return A function that maps a numeric p-value vector to a scalar estimate of pi0.
#' @export
make_base_pi0 <- function(type = c("ecdf_storey_median", "ecdf_storey_weighted",
                                   "cp4p_st.boot","cp4p_st.spline","cp4p_langaas","cp4p_jiang",
                                   "cp4p_histo","cp4p_pounds","cp4p_abh","cp4p_slim"),
                          Lambda = NULL, eps = 1e-8,
                          nbins = 20, pz = 0.05, censor = FALSE,
                          cp4p_fallback = c("storey_median","storey_weighted")) {
  type <- match.arg(type)
  if (grepl("^ecdf_", type)) {
    if (is.null(Lambda)) {
      Lambda <- if (type == "ecdf_storey_median") c(0.6, 0.7, 0.8) else c(0.6, 0.7, 0.8, 0.9)
    }
    if (type == "ecdf_storey_median") {
      return(function(p) storey_median_lambda(p, Lambda, eps))
    } else {
      return(function(p) storey_weighted_multi(p, Lambda, eps))
    }
  } else {
    method <- sub("^cp4p_", "", type)
    return(make_cp4p_base(method, nbins, pz, censor, match.arg(cp4p_fallback), Lambda, eps))
  }
}

storey_fixed <- function(p, lambda = 0.5, eps = 1e-8) {
  m <- length(p)
  clip01(sum(p > lambda) / ((1 - lambda) * m), eps)
}

storey_median_lambda <- function(p, Lambda = c(0.6, 0.7, 0.8), eps = 1e-8) {
  ests <- vapply(Lambda, function(l) storey_fixed(p, l, eps), numeric(1))
  clip01(stats::median(ests), eps)
}

storey_weighted_multi <- function(p, Lambda = c(0.6, 0.7, 0.8, 0.9), eps = 1e-8) {
  w <- (1 - Lambda); w <- w / sum(w)
  ests <- vapply(Lambda, function(l) storey_fixed(p, l, eps), numeric(1))
  clip01(sum(w * ests), eps)
}

make_cp4p_base <- function(method = c("st.boot","st.spline","langaas","jiang","histo","pounds","abh","slim"),
                           nbins = 20, pz = 0.05, censor = FALSE,
                           fallback = c("storey_median","storey_weighted"),
                           Lambda = NULL, eps = 1e-8) {
  method <- match.arg(method)
  fallback <- match.arg(fallback)
  if (is.null(Lambda)) {
    Lambda <- if (fallback == "storey_median") c(0.6, 0.7, 0.8) else c(0.6, 0.7, 0.8, 0.9)
  }
  fb_fun <- if (fallback == "storey_median") {
    function(p) storey_median_lambda(p, Lambda = Lambda, eps = eps)
  } else {
    function(p) storey_weighted_multi(p, Lambda = Lambda, eps = eps)
  }
  if (!isTRUE(getOption("have_cp4p", FALSE)) && !exists("have_cp4p")) {
  warning("cp4p not installed; using ECDF fallback base estimator.")
  return(fb_fun)
}
  function(p) {
    p_in <- if (isTRUE(censor)) p else c(p, 1)
    res <- try(cp4p::estim.pi0(p_in, pi0.method = method, nbins = nbins, pz = pz), silent = TRUE)
    if (inherits(res, "try-error") || is.null(res$pi0) || !is.finite(res$pi0)) {
      fb_fun(p)
    } else {
      clip01(as.numeric(res$pi0), eps)
    }
  }
}

dkw_epsilon_one_sided <- function(delta, n) {
  sqrt(log(1 / delta) / (2 * n))
}

#' ECDF + one-sided DKW upper bound on pi0
#'
#' @param p Numeric vector of p-values.
#' @param H Tail bandwidths for the grid.
#' @param delta_total Total miscoverage budget across folds and H.
#' @param K Number of folds.
#'
#' @return Scalar upper bound on pi0.
#' @export
safety_pi0_upper <- function(p, H = c(0.05, 0.10, 0.20), delta_total = 0.02, K = 10) {
  n <- length(p); if (n <= 0) return(1)
  Fhat <- stats::ecdf(p)
  delta_per_h <- delta_total / (K * length(H))
  vals <- sapply(H, function(h) {
    num <- 1 - Fhat(1 - h) + dkw_epsilon_one_sided(delta_per_h, n)
    clip01(num / h, 1e-12)
  })
  clip01(min(vals), 1e-12)
}


#' Cross-fitted adaptive FDR with Safe-U denominator
#'
#' Implements a cross-fitted adaptive FDR procedure in which each fold
#' uses an upper bound on the null proportion \code{pi0} estimated from
#' the complementary folds. The denominator in the BH threshold within
#' each fold is taken to be an upper bound \code{U_use} on \code{pi0},
#' combining:
#'
#' \itemize{
#'   \item an ECDF-based one-sided DKW upper bound (\code{U_safe}),
#'   \item and, optionally, a conformal upper bound derived from a
#'         conditional quantile regression fit (\code{U_conf}).
#' }
#'
#' The overall denominator is \code{U_use = min(U_safe, U_conf)}, clipped
#' to \code{[eps, 1]}. When \code{use_conformal = FALSE}, only
#' \code{U_safe} is used.
#'
#' @param pvals Numeric vector of p-values.
#' @param alpha Target FDR level (e.g. 0.1).
#' @param K Number of folds for cross-fitting.
#' @param base_pi0 Base estimator factory (function) used only to compute
#'   a scalar feature \code{A} in each fold. Defaults to
#'   \code{make_base_pi0("ecdf_storey_median")}.
#' @param delta_ecdf Total miscoverage budget for the ECDF-based upper
#'   bound, spread across folds and bandwidths.
#' @param H Numeric vector of tail bandwidths for the ECDF upper bound.
#' @param use_conformal Logical; if \code{TRUE}, use \code{cqr_fit} and
#'   \code{cal_df} to compute an additional conformal upper bound.
#' @param cqr_fit Fitted object from [scoutFDR::train_cqr_pi0()], or \code{NULL}
#'   if conformal tightening is not used.
#' @param cal_df Calibration data frame used to fit the CQR model.
#' @param feature_names Character vector of feature names to feed into
#'   the CQR model (must be present in \code{cal_df} and constructed by
#'   \code{feature_fun}).
#' @param feature_fun Function that maps a numeric vector of p-values
#'   to a named list of features. Defaults to [scoutFDR::features_ecdf()].
#' @param ood_params List of parameters controlling out-of-distribution
#'   checks and conformal correction. Recognised elements:
#'   \itemize{
#'     \item \code{max_weight}: maximum density-ratio weight (default 2),
#'     \item \code{maha_alpha}: significance level for Mahalanobis
#'           outlier detection (default 0.01),
#'     \item \code{alpha_cqr}: miscoverage level for the conformal
#'           correction (default 0.05).
#'   }
#' @param folds Optional integer vector of fold labels of the same length
#'   as \code{p}. If \code{NULL}, a random assignment into \code{K}
#'   approximately equal folds is used.
#' @param use_BY Logical; if \code{TRUE}, replace BH thresholds with
#'   Benjamini--Yekutieli thresholds within each fold.
#' @param eps Numerical stability parameter for bounding probabilities.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{rejections}: integer vector of indices of rejected
#'         hypotheses,
#'   \item \code{U_fold}: numeric vector of length \code{K} giving the
#'         per-fold \code{U_use},
#'   \item \code{folds}: the fold assignment used,
#'   \item \code{mode}: a list summarising whether Safe-U and conformal
#'         tightening were used.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' pvals <- simulate_pvalues_ttest(m = 2000, n = 5, pi0 = 0.2, target_power = 0.35)
#' res <- adaptive_fdr(pvals, alpha = 0.10, K = 8)
#' length(res$rejections)
#' mean(res$U_fold, na.rm = TRUE)
#' }
#'
#' @export
adaptive_fdr <- function(
    pvals,
    alpha = 0.10,
    K = 10,
    base_pi0 = make_base_pi0("ecdf_storey_median"),
    delta_ecdf = 0.02,
    H = c(0.05, 0.10, 0.20),
    use_conformal = FALSE,
    cqr_fit = NULL,
    cal_df = NULL,
    feature_names = c("A", "x9", "x10"),
    feature_fun = function(pp) as.list(features_ecdf(pp)),
    ood_params = list(max_weight = 2, maha_alpha = 0.01, alpha_cqr = 0.05),
    folds = NULL,
    use_BY = FALSE,
    eps = 1e-8
) {
  stopifnot(is.numeric(pvals), length(pvals) > 0)
  m <- length(pvals)

  # Fold assignment
  if (is.null(folds)) {
    folds <- sample(rep(seq_len(K), length.out = m))
  } else {
    K <- max(folds)
  }

  # Precompute target features for density-ratio (including A per fold)
  X_target_all <- NULL
  if (use_conformal && !is.null(cqr_fit) && !is.null(cal_df)) {
    X_target_all <- do.call(rbind, lapply(seq_len(K), function(ell) {
      idx_comp <- which(folds != ell)
      feats_list <- feature_fun(pvals[idx_comp])
      A_tmp <- base_pi0(pvals[idx_comp])
      df <- as.data.frame(as.list(c(A = A_tmp, unlist(feats_list))))
      df[, feature_names, drop = FALSE]
    }))
  }

  # Mahalanobis out-of-distribution check (internal)
  mahalanobis_ood <- function(x, Xref, alpha = 0.01) {
    mu <- colMeans(Xref, na.rm = TRUE)
    S <- stats::cov(Xref, use = "pairwise.complete.obs")

    Sinv <- try(solve(S + diag(1e-6, nrow(S))), silent = TRUE)
    if (inherits(Sinv, "try-error")) {
      return(list(is_ood = TRUE))
    }

    d2 <- as.numeric(t(x - mu) %*% Sinv %*% (x - mu))
    crit <- stats::qchisq(1 - alpha, df = ncol(Xref))

    list(is_ood = (!is.finite(d2)) || (d2 > crit))
  }

  rejections <- integer(0)
  U_fold <- rep(NA_real_, K)

  for (ell in seq_len(K)) {
    idx_comp <- which(folds != ell)
    idx_fold <- which(folds == ell)
    p_comp <- pvals[idx_comp]
    p_fold <- pvals[idx_fold]

    if (length(p_fold) == 0) {
      next
    }

    # Safe-U from ECDF + DKW
    U_safe <- safety_pi0_upper(
      p_comp,
      H = H,
      delta_total = delta_ecdf,
      K = K
    )

    # Optional conformal tightening
    U_conf <- 1
    if (use_conformal && !is.null(cqr_fit) && !is.null(cal_df) && !is.null(X_target_all)) {
      feats0 <- feature_fun(p_comp)
      A_minus <- base_pi0(p_comp)
      x0_full <- c(A = A_minus, unlist(feats0))
      x0 <- as.numeric(x0_full[feature_names])
      names(x0) <- feature_names

      ood <- mahalanobis_ood(
        x = x0,
        Xref = as.matrix(cal_df[, feature_names, drop = FALSE]),
        alpha = ood_params$maha_alpha %||% 0.01
      )

      if (!ood$is_ood) {
        U_conf <- cqr_upper_bound(
          cqr_fit = cqr_fit,
          cal_df = cal_df,
          feature_names = feature_names,
          x0 = x0,
          X_target_all = X_target_all,
          alpha = ood_params$alpha_cqr %||% 0.05,
          weight_clip = c(0.2, ood_params$max_weight %||% 2)
        )
      } else {
        U_conf <- 1
      }
    }

    # Combined bound and storage
    U_use <- min(U_safe, U_conf)
    U_use <- clip01(U_use, eps)
    U_fold[ell] <- U_use

    # BH (or BY) within fold with denominator U_use
    ord <- order(p_fold)
    p_sorted <- p_fold[ord]
    k_seq <- seq_along(p_sorted)

    if (!use_BY) {
      thr <- (k_seq * alpha) / (m * U_use)
    } else {
      c_m <- sum(1 / (1:length(p_fold)))
      thr <- (k_seq * alpha) / (m * U_use * c_m)
    }

    pass <- which(p_sorted <= thr)
    if (length(pass) > 0) {
      kmax <- max(pass)
      tval <- p_sorted[kmax]
      rejections <- c(rejections, idx_fold[p_fold <= tval])
    }

    if (getOption("scoutFDR.debug", FALSE)) {
      message(
        sprintf(
          "Fold %d: U_safe=%.3f, U_conf=%.3f, U_use=%.3f",
          ell, U_safe, U_conf, U_use
        )
      )
    }
  }

  list(
    rejections = sort(unique(rejections)),
    U_fold = U_fold,
    folds = folds,
    mode = list(
      safeU = TRUE,
      conformal = use_conformal && !is.null(cqr_fit) && !is.null(cal_df)
    )
  )
}
