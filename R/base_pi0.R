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
  if (!have_cp4p) {
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
