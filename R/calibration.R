# ============================================================
# Calibration dataset generation, base pi0 estimators,
# CQR training, conformal upper bounds, adaptive FDR
# ============================================================

utils::globalVariables(c("have_cp4p"))

# ============================================================
# Base pi0 estimators (Storey-type and CP4P)
# ============================================================

storey_fixed <- function(pvals, lambda = 0.5, eps = 1e-8) {
  m <- length(pvals)
  val <- sum(pvals > lambda) / ((1 - lambda) * m)
  pmin(pmax(val, eps), 1)
}

storey_median_lambda <- function(pvals, Lambda = c(0.6, 0.7, 0.8), eps = 1e-8) {
  ests <- vapply(Lambda, function(l) storey_fixed(pvals, l, eps), numeric(1))
  pmin(pmax(stats::median(ests), eps), 1)
}

storey_weighted_multi <- function(pvals, Lambda = c(0.6, 0.7, 0.8, 0.9), eps = 1e-8) {
  w <- (1 - Lambda); w <- w / sum(w)
  ests <- vapply(Lambda, function(l) storey_fixed(pvals, l, eps), numeric(1))
  pmin(pmax(sum(w * ests), eps), 1)
}

make_cp4p_base <- function(
    method = c("st.boot","st.spline","langaas","jiang","histo","pounds","abh","slim"),
    nbins = 20, pz = 0.05, censor = FALSE,
    fallback = c("storey_median","storey_weighted"),
    Lambda = NULL, eps = 1e-8
) {
  method <- match.arg(method)
  fallback <- match.arg(fallback)

  if (is.null(Lambda)) {
    Lambda <- if (fallback == "storey_median") c(0.6, 0.7, 0.8)
    else c(0.6, 0.7, 0.8, 0.9)
  }

  fb_fun <- if (fallback == "storey_median") {
    function(pvals) storey_median_lambda(pvals, Lambda, eps)
  } else {
    function(pvals) storey_weighted_multi(pvals, Lambda, eps)
  }

  if (!have_cp4p) {
    warning("cp4p not installed; using ECDF fallback base estimator.")
    return(fb_fun)
  }

  function(pvals) {
    p_in <- if (isTRUE(censor)) pvals else c(pvals, 1)
    res <- try(cp4p::estim.pi0(p_in, pi0.method = method, nbins = nbins, pz = pz),
               silent = TRUE)
    if (inherits(res, "try-error") || is.null(res$pi0) || !is.finite(res$pi0)) {
      fb_fun(pvals)
    } else {
      pmin(pmax(as.numeric(res$pi0), eps), 1)
    }
  }
}

#' Construct a base \eqn{\pi_0} estimator function
#'
#' @description
#' Creates and returns a function that estimates the null proportion π₀
#' from a vector of p-values.  
#'
#' There are two types of estimators:
#' \itemize{
#'   \item \strong{ECDF-based}: Storey-style estimators using multiple
#'         λ thresholds.
#'   \item \strong{cp4p-based}: Flexible density-estimation methods from the
#'         \pkg{cp4p} package, with automatic fallback if that package is not
#'         available or if the estimate fails.
#' }
#' 
#' @param type Character string specifying the estimator. One of:
#'   \itemize{
#'     \item `"ecdf_storey_median"` — median of Storey(\eqn{\lambda}) over \eqn{\lambda \in \{0.6, 0.7, 0.8\}}
#'     \item `"ecdf_storey_weighted"` — weighted Storey(\eqn{\lambda}) over \eqn{\lambda \in \{0.6, 0.7, 0.8, 0.9\}}
#'     \item `"cp4p_st.boot"`, `"cp4p_st.spline"`, `"cp4p_langaas"`,
#'           `"cp4p_jiang"`, `"cp4p_histo"`, `"cp4p_pounds"`,
#'           `"cp4p_abh"`, `"cp4p_slim"`
#'   }
#' @param Lambda Optional numeric vector of λ values used for ECDF-based
#'   Storey estimators. If `NULL`, a default set is chosen depending on `type`.
#' @param eps Numeric stability value for clipping the estimator to the interval
#'   \eqn{[eps, 1]}. Default: `1e-8`.
#' @param nbins Number of histogram bins for cp4p estimators. Default: 20.
#' @param pz Zero threshold (argument passed to \code{cp4p}). Default: 0.05.
#' @param censor Logical. If `TRUE`, censor p-values at 1 before passing to cp4p.
#'   Default: `FALSE`.
#' @param cp4p_fallback Character string specifying which ECDF estimator to use
#'   if cp4p fails (`"storey_median"` or `"storey_weighted"`).
#'
#' @return A function with signature `function(pvals)` returning a scalar
#'   \eqn{\pi_0} estimate in the interval \eqn{[eps, 1]}.
#'
#' @export
make_base_pi0 <- function(
    type = c("ecdf_storey_median", "ecdf_storey_weighted",
             "cp4p_st.boot","cp4p_st.spline","cp4p_langaas","cp4p_jiang",
             "cp4p_histo","cp4p_pounds","cp4p_abh","cp4p_slim"),
    Lambda = NULL, eps = 1e-8,
    nbins = 20, pz = 0.05, censor = FALSE,
    cp4p_fallback = c("storey_median","storey_weighted")
) {
  type <- match.arg(type)

  if (grepl("^ecdf_", type)) {
    if (is.null(Lambda))
      Lambda <- if (type == "ecdf_storey_median")
        c(0.6,0.7,0.8) else c(0.6,0.7,0.8,0.9)

    if (type == "ecdf_storey_median")
      return(function(pvals) storey_median_lambda(pvals, Lambda, eps))
    else
      return(function(pvals) storey_weighted_multi(pvals, Lambda, eps))
  }

  method <- sub("^cp4p_", "", type)
  make_cp4p_base(method, nbins, pz, censor,
                 match.arg(cp4p_fallback), Lambda, eps)
}


# ============================================================
# Generate calibration dataset
# ============================================================

#' Generate calibration data for pi0 regression
#'
#' @description
#' Repeatedly simulates p-values under a t-test model, computes a base
#' pi0 estimate, extracts features, and returns a calibration data frame.
#'
#' @param n_sim Number of simulated datasets.
#' @param grid Optional data.frame of (pi0, m, n, target_power).
#' @param base_pi0 Function: base pi0 estimator.
#' @param feature_fun Function: maps pvals -> feature list.
#' @param seed Random seed.
#'
#' @return A numeric data frame, one row per simulation.
#'
#' @export
generate_calibration_data <- function(
    n_sim = 2000,
    grid = NULL,
    base_pi0 = make_base_pi0("ecdf_storey_median"),
    feature_fun = function(pvals) as.list(features_ecdf(pvals)),
    seed = 1
) {
  if (is.null(grid)) {
    grid <- expand.grid(
      pi0 = seq(0.1, 0.9, by = 0.05),
      m = seq(500, 12000, by = 100),
      n = 3:12,
      target_power = seq(0.05, 0.7, by = 0.05)
    )
  }

  set.seed(seed)
  out <- vector("list", n_sim)
  nG <- nrow(grid)

  pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
  on.exit(close(pb), add = TRUE)

  for (i in seq_len(n_sim)) {
    g <- grid[sample.int(nG, 1), ]

    pvals <- simulate_pvalues_ttest(
      m = g$m, n = g$n,
      pi0 = g$pi0, target_power = g$target_power
    )

    A <- base_pi0(pvals)
    feats <- feature_fun(pvals)

    row <- c(pi0 = g$pi0, A = as.numeric(A), unlist(feats))
    storage.mode(row) <- "double"
    out[[i]] <- row

    utils::setTxtProgressBar(pb, i)
  }

  df <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
  for (nm in names(df)) df[[nm]] <- as.numeric(df[[nm]])
  df
}


# ============================================================
# CQR training + prediction
# ============================================================

scout_select_stable_features <- function(cal_df, probe = TRUE, K = 2) {
  proto_cols <- names(as.list(features_ecdf(stats::runif(2000))))
  cand <- setdiff(intersect(colnames(cal_df), proto_cols), "pi0")

  # Calibration stability filter
  keep_cal <- vapply(cand, function(nm) {
    v <- cal_df[[nm]]
    all(is.finite(v)) && stats::sd(v) > 1e-10
  }, logical(1))

  if (!probe) return(cand[keep_cal])

  # Small probe scenario (tests stability in "real" SCOUT usage)
  p_probe <- simulate_pvalues_ttest(m = 5000, n = 6, pi0 = 0.4, target_power = 0.35)
  fold_id <- sample(rep(seq_len(K), length.out = length(p_probe)))

  keep_tar <- vapply(cand, function(nm) {
    for (k in seq_len(K)) {
      fk <- try(features_ecdf(p_probe[fold_id == k]), silent = TRUE)
      if (inherits(fk, "try-error") || !(nm %in% names(fk))) return(FALSE)
      if (!is.finite(fk[[nm]])) return(FALSE)
    }
    TRUE
  }, logical(1))

  cand[keep_cal & keep_tar]
}


#' Train a conditional quantile regression model for \eqn{\pi_0}
#'
#' @description
#' Fits a model predicting the null proportion \eqn{\pi_0} from calibration
#' features, using either:
#'
#' - **quantile regression** via the `quantreg` package (`model = "qr"`)
#' - **quantile regression forest** via the `quantregForest` package (`model = "qrf"`)
#'
#' The function automatically:
#' - keeps only available feature columns,
#' - removes rows with missing values,
#' - removes near-zero-variance predictors,
#' - enforces full column rank (QR pivoting),
#' - standardizes predictors when using quantile regression,
#' - falls back to quantile forests or linear regression if needed.
#'
#' This model is later used inside [`adaptive_fdr()`] to produce conformal
#' upper bounds for \eqn{\pi_0}.
#'
#' @param cal_df A data.frame containing at least:
#'   - a column named `"pi0"` with the true simulation π₀,
#'   - feature columns named in `feature_names`.
#'
#' @param feature_names Character vector naming the predictor columns from
#'   `cal_df` to use in the model. Default `"auto"` allows for automatic selection of stable features.
#'
#' @param tau Numeric scalar in \eqn{(0,1)}; the quantile level to be
#'   estimated. Defaults to `0.8`, following conformal upper-bound practice.
#'
#' @param model Character string, `"qr"` (default) for quantile regression
#'   or `"qrf"` for quantile regression forest. If `"qr"` fails, the function
#'   will attempt a fallback to `"qrf"` (if available), and otherwise to `lm()`.
#'
#' @return A list describing the fitted model, containing:
#' \describe{
#'   \item{type}{One of `"qr"`, `"qrf"`, or `"lm"`.}
#'   \item{feature_names}{Character vector of features actually used.}
#'   \item{tau}{The quantile level used.}
#'
#'   \item{rqfit}{(for `"qr"`) fitted `quantreg::rq` object.}
#'   \item{mu}{(for `"qr"`) column means used for scaling.}
#'   \item{sd}{(for `"qr"`) column standard deviations used for scaling.}
#'   \item{vars}{(for `"qr"`) names of standardized predictor columns.}
#'
#'   \item{qrf}{(for `"qrf"`) fitted `quantregForest` model.}
#'
#'   \item{lmfit}{(fallback) linear regression fit if both QR and QRF are unavailable.}
#' }
#'
#' @export
train_cqr_pi0 <- function(cal_df, feature_names = "auto",
                          tau = 0.8,
                          model = c("qr", "qrf")) {
  model <- match.arg(model)

  # Automatic stable feature selection
  if (identical(feature_names, "auto") ||
      is.null(feature_names) ||
      identical(feature_names, "")) {

    feature_names <- scout_select_stable_features(cal_df, probe = TRUE, K = 2)
    message("train_cqr_pi0(): selected ", length(feature_names), " stable features.")
  }

  keep <- intersect(feature_names, colnames(cal_df))
  df <- cal_df[, c("pi0", keep), drop = FALSE]
  df <- df[stats::complete.cases(df),, drop = FALSE]

  if (nrow(df) < 100)
    stop("Too few complete calibration rows.")

  for (nm in names(df)) df[[nm]] <- as.numeric(df[[nm]])

  y <- df$pi0
  X <- as.matrix(df[, keep, drop = FALSE])

  vars <- apply(X, 2, function(xx) stats::var(xx, na.rm = TRUE))
  keep_cols <- names(vars)[is.finite(vars) & vars > 1e-10]
  X <- X[, keep_cols, drop = FALSE]

  QR <- qr(X)
  rank <- QR$rank
  piv <- QR$pivot[seq_len(rank)]
  X <- X[, piv, drop = FALSE]

  kept_names <- colnames(X)

  # Drop low-variance surplus predictors
  if (nrow(X) < ncol(X) + 20) {
    v <- apply(X, 2, stats::var)
    ord <- order(v, decreasing = TRUE)
    keep_k <- max(1L, nrow(X) - 20L)
    X <- X[, ord[seq_len(keep_k)], drop = FALSE]
    kept_names <- colnames(X)
  }

  X_mu <- colMeans(X)
  X_sd <- apply(X, 2, sd)
  X_sd[X_sd < 1e-8] <- 1
  Xs <- scale(X, center = X_mu, scale = X_sd)

  if (model == "qr") {
    if (!requireNamespace("quantreg", quietly = TRUE))
      stop("quantreg not installed.")

    df_fit <- data.frame(y = y, Xs)
    vars <- colnames(Xs)
    form <- stats::as.formula(paste("y ~", paste(vars, collapse = " + ")))

    rqfit <- try(
      quantreg::rq(form, data = df_fit, tau = tau, method = "fn"),
      silent = TRUE
    )

    if (!inherits(rqfit, "try-error")) {
      return(list(
        type = "qr",
        rqfit = rqfit,
        mu = X_mu,
        sd = X_sd,
        feature_names = kept_names,
        vars = vars,
        tau = tau
      ))
    }

    # fallback to QRF
    if (requireNamespace("quantregForest", quietly = TRUE)) {
      qrf <- quantregForest::quantregForest(X, y, ntree = 500)
      return(list(type = "qrf", qrf = qrf,
                  feature_names = kept_names, tau = tau))
    }

    lmfit <- stats::lm(y ~ ., data = data.frame(y = y, X))
    return(list(
      type = "lm", lmfit = lmfit, feature_names = kept_names,
      mu = rep(0, length(kept_names)), sd = rep(1, length(kept_names)),
      tau = tau
    ))
  }

  # quantile forest
  if (!requireNamespace("quantregForest", quietly = TRUE))
    stop("quantregForest not installed.")

  qrf <- quantregForest::quantregForest(X, y, ntree = 500)
  list(type = "qrf", qrf = qrf, feature_names = kept_names, tau = tau)
}

predict_qtau <- function(cqr_fit, x_vec) {
  x <- as.numeric(x_vec[cqr_fit$feature_names])

  if (cqr_fit$type == "qr") {
    xs <- (x - cqr_fit$mu) / cqr_fit$sd
    newd <- as.data.frame(matrix(xs, nrow = 1))
    colnames(newd) <- cqr_fit$vars
    return(as.numeric(stats::predict(cqr_fit$rqfit, newd)))
  }

  if (cqr_fit$type == "qrf") {
    X_new <- matrix(x, nrow = 1)
    colnames(X_new) <- cqr_fit$feature_names
    return(as.numeric(predict(cqr_fit$qrf, X_new, what = cqr_fit$tau)))
  }

  # fallback linear model
  newd <- as.data.frame(matrix(x, nrow = 1))
  colnames(newd) <- cqr_fit$feature_names
  as.numeric(stats::predict(cqr_fit$lmfit, newd))
}


# ============================================================
# Conformal upper bound
# ============================================================

#' Conformal Upper Bound for \eqn{\pi_0} Using CQR Residuals
#'
#' @description
#' Computes a conformal upper bound on the proportion of true null hypotheses
#' (\code{pi0}) at a new feature vector \code{x0}.  
#'
#' This function implements the weighted conformal prediction method described
#' in the SCOUT-FDR framework:
#'
#' 1. Balance calibration and target samples.  
#' 2. Estimate a density ratio via logistic regression to obtain importance
#'    weights.  
#' 3. Compute conformal residuals from the CQR model on the calibration sample.  
#' 4. Form a weighted empirical quantile of the residuals at level \code{1 - alpha}.  
#' 5. Add this quantile to the CQR prediction at \code{x0} to obtain the final
#'    upper bound.  
#'
#' @param cqr_fit A fitted object returned by
#'   \code{\link{train_cqr_pi0}}. Contains either a quantile regression
#'   model or a quantile regression forest.
#'
#' @param cal_df Calibration data frame. Must contain a column \code{pi0}
#'   and all feature columns listed in \code{feature_names}.
#'
#' @param feature_names Character vector of feature names used by the CQR model.
#'   Must appear as column names in both \code{cal_df} and \code{X_target_all}.
#'
#' @param x0 Numeric vector of features for a single target point (a single
#'   p-value distribution). Should be named and match \code{feature_names}.
#'
#' @param X_target_all Data frame containing one or more new feature vectors at
#'   which upper bounds may be needed. Only the rows are used to match the
#'   calibration sample size.
#'
#' @param alpha Miscoverage level. The returned upper bound satisfies
#'   (approximately) \code{P(pi0 <= U) >= 1 - alpha}.
#'
#' @param weight_clip Numeric vector of length two specifying the lower and
#'   upper bounds used to truncate density-ratio weights. Defaults to
#'   \code{c(0.2, 2)} to avoid extreme weights.
#'
#' @details
#' If logistic regression fails or yields unstable probabilities, the method
#' automatically falls back to uniform weighting.  
#'
#' If the conformal residual sample is degenerate (too few finite values), the
#' function defaults to returning \code{1}.
#'
#' @return
#' A single numeric value in \code{[0, 1]}, representing a conformal upper bound
#' on \eqn{\pi_0(x0)}.
#'
#' @export
#'
#' @importFrom stats glm predict approx sd setNames
cqr_upper_bound <- function(
    cqr_fit, cal_df, feature_names,
    x0, X_target_all,
    alpha = 0.05,
    weight_clip = c(0.2, 2)
) {

  X_cal_full <- as.matrix(cal_df[, feature_names, drop = FALSE])
  X_tar_full <- as.matrix(X_target_all[, feature_names, drop = FALSE])

  n_c <- nrow(X_cal_full)
  n_t <- nrow(X_tar_full)
  if (n_c == 0 || n_t == 0) return(1)

  # balance sample sizes
  if (n_c > n_t) {
    sel_c <- sample.int(n_c, n_t)
    X_cal <- X_cal_full[sel_c,, drop = FALSE]
    X_tar <- X_tar_full
  } else {
    sel_c <- seq_len(n_c)
    X_cal <- X_cal_full
    X_tar <- X_tar_full[sample.int(n_t, n_c),, drop = FALSE]
  }

  # standardize
  X_all <- rbind(X_cal, X_tar)
  mu <- colMeans(X_all, na.rm = TRUE)
  sdv <- apply(X_all, 2, sd)
  sdv[sdv < 1e-8] <- 1

  X_cal_s <- scale(X_cal, center = mu, scale = sdv)
  X_tar_s <- scale(X_tar, center = mu, scale = sdv)

  # logistic regression density ratio
  df_lr <- data.frame(z = c(rep(0L, nrow(X_cal_s)),
                            rep(1L, nrow(X_tar_s))),
                      rbind(X_cal_s, X_tar_s))
  colnames(df_lr) <- c("z", feature_names)

  fit_ok <- TRUE
  p_cal <- NULL
  suppressWarnings({
    gl <- try(stats::glm(z ~ ., data = df_lr, family = stats::binomial(),
                         control = list(maxit = 200)), silent = TRUE)
    if (inherits(gl, "try-error")) fit_ok <- FALSE else {
      newd <- as.data.frame(X_cal_s); colnames(newd) <- feature_names
      p_cal <- try(stats::predict(gl, newdata = newd, type = "response"),
                   silent = TRUE)
      if (inherits(p_cal, "try-error") || any(!is.finite(p_cal)))
        fit_ok <- FALSE
    }
  })

  if (fit_ok &&
      (min(p_cal, na.rm = TRUE) < 1e-6 ||
       max(p_cal, na.rm = TRUE) > 1 - 1e-6)) {
    fit_ok <- FALSE
  }

  if (fit_ok) {
    w_raw <- p_cal / pmax(1 - p_cal, 1e-12)
    w <- pmin(pmax(w_raw, weight_clip[1]), weight_clip[2])
  } else {
    w <- rep(1, nrow(X_cal_s))
  }

  # Conformal residuals
  qhat_cal <- apply(
    X_cal,
    1,
    function(xx) predict_qtau(cqr_fit, setNames(as.numeric(xx), feature_names))
  )
  y_bal <- cal_df$pi0[sel_c]
  r <- y_bal - qhat_cal

  ok <- is.finite(r) & is.finite(w)
  r <- r[ok]; w <- w[ok]

  if (length(r) < 2) return(1)

  ord <- order(r)
  r_ord <- r[ord]
  w_ord <- w[ord]
  cw <- cumsum(w_ord) / sum(w_ord)

  x_vec <- c(0, cw)
  y_vec <- c(r_ord[1], r_ord)

  ok2 <- is.finite(x_vec) & is.finite(y_vec)
  x_vec <- x_vec[ok2]; y_vec <- y_vec[ok2]
  if (length(x_vec) < 2) return(1)

  c_hat <- stats::approx(
    x_vec, y_vec,
    xout = 1 - alpha,
    ties = "ordered",
    rule = 2
  )$y

  c_hat <- max(c_hat, 0)

  q_x0 <- predict_qtau(cqr_fit, setNames(as.numeric(x0), feature_names))
  if (!is.finite(q_x0)) return(1)

  pmin(pmax(q_x0 + c_hat, 0), 1)
}

# ============================================================
# SafeU upper bound (ECDF)
# ============================================================

dkw_epsilon_one_sided <- function(delta, n) {
  sqrt(log(1 / delta) / (2 * n))
}

safety_pi0_upper <- function(
    pvals,
    H = c(0.05, 0.10, 0.20),
    delta_total = 0.02,
    K = 10
) {
  n <- length(pvals)
  Fhat <- stats::ecdf(pvals)

  delta_per_h <- delta_total / (K * length(H))

  vals <- sapply(H, function(h) {
    num <- 1 - Fhat(1 - h) + dkw_epsilon_one_sided(delta_per_h, n)
    pmin(pmax(num / h, 1e-12), 1)
  })

  min(vals)
}
