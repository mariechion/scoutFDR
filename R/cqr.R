utils::globalVariables(c("have_qr", "have_qrf"))

#' Train conditional quantile regression for pi0
#'
#' Fits a conditional quantile model for \code{pi0} given features, using
#' linear quantile regression (via \pkg{quantreg}) when available, and
#' falling back to quantile regression forests (\pkg{quantregForest})
#' or an ordinary linear model.
#'
#' This function assumes that \code{cal_df} contains a column named
#' \code{pi0} and the specified \code{feature_names}. It performs basic
#' preprocessing (dropping constant columns, QR-based column selection,
#' and standardisation for linear quantile regression).
#'
#' @importFrom stats as.formula approx
#' @param cal_df Data frame of calibration samples, typically produced by
#'   [generate_calibration_data()]. Must contain column \code{pi0} and the
#'   specified features.
#' @param feature_names Character vector giving the names of features in
#'   \code{cal_df} to be considered as predictors.
#' @param tau Numeric in \eqn{[0,1]}, the quantile level to model (e.g. 0.8).
#' @param model Character string, either \code{"qr"} for linear quantile
#'   regression or \code{"qrf"} for quantile regression forests.
#'
#' @return A list describing the fitted model, with at least the elements
#'   \code{type} (one of \code{"qr"}, \code{"qrf"}, \code{"lm"}),
#'   \code{feature_names}, and model-specific components.
#' @export
train_cqr_pi0 <- function(cal_df, feature_names, tau = 0.8, model = c("qr", "qrf")) {
  model <- match.arg(model)

  # Keep only pi0 and requested features
  keep <- intersect(feature_names, colnames(cal_df))
  df <- cal_df[, c("pi0", keep), drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) < 100) {
    stop("Too few complete calibration rows after filtering.")
  }

  # Coerce all columns to numeric
  for (nm in colnames(df)) {
    df[[nm]] <- as.numeric(df[[nm]])
  }

  y <- as.numeric(df$pi0)
  X <- as.matrix(df[, keep, drop = FALSE])

  # Drop (near-)zero variance columns
  vars <- apply(X, 2, function(col) stats::var(col, na.rm = TRUE))
  keep_cols <- names(vars)[is.finite(vars) & vars > 1e-10]
  X <- X[, keep_cols, drop = FALSE]

  # QR decomposition to remove collinearity
  QR <- qr(X)
  rank <- QR$rank
  piv <- QR$pivot[seq_len(rank)]
  X <- X[, piv, drop = FALSE]
  kept_names <- colnames(X)

  # If still too many predictors, keep the most variable ones
  if (nrow(X) < ncol(X) + 20) {
    v <- apply(X, 2, stats::var)
    ord <- order(v, decreasing = TRUE)
    keep_k <- max(1L, nrow(X) - 20L)
    X <- X[, ord[seq_len(keep_k)], drop = FALSE]
    kept_names <- colnames(X)
  }

  # Standardisation for linear QR
  X_mu <- colMeans(X)
  X_sd <- apply(X, 2, stats::sd)
  X_sd[X_sd < 1e-8] <- 1
  Xs <- scale(X, center = X_mu, scale = X_sd)

  if (model == "qr") {
    if (!isTRUE(have_qr)) {
      stop("quantreg not installed; cannot fit 'qr' model.")
    }
    df_fit <- data.frame(y = y, Xs)
    vars <- colnames(Xs)
    form <- as.formula(paste("y ~", paste(vars, collapse = " + ")))

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

    # Fall through to forest or lm if QR fails
    if (isTRUE(have_qrf)) {
      qrf <- quantregForest::quantregForest(X, y, ntree = 500)
      return(list(
        type = "qrf",
        qrf = qrf,
        feature_names = kept_names,
        tau = tau
      ))
    } else {
      lmfit <- stats::lm(y ~ ., data = data.frame(y = y, X))
      return(list(
        type = "lm",
        lmfit = lmfit,
        feature_names = kept_names,
        mu = rep(0, length(kept_names)),
        sd = rep(1, length(kept_names)),
        tau = tau
      ))
    }
  } else {
    # model == "qrf"
    if (!isTRUE(have_qrf)) {
      stop("quantregForest not installed; cannot fit 'qrf' model.")
    }
    qrf <- quantregForest::quantregForest(X, y, ntree = 500)
    return(list(
      type = "qrf",
      qrf = qrf,
      feature_names = kept_names,
      tau = tau
    ))
  }
}

# Internal helper: predict quantile for a single feature vector
predict_qtau <- function(cqr_fit, x_vec) {
  x <- as.numeric(x_vec[cqr_fit$feature_names])

  if (cqr_fit$type == "qr") {
    # Standardise using training means/sds
    xs <- (x - cqr_fit$mu) / cqr_fit$sd
    newd <- as.data.frame(matrix(xs, nrow = 1))
    colnames(newd) <- cqr_fit$vars
    as.numeric(stats::predict(cqr_fit$rqfit, newd))

  } else if (cqr_fit$type == "qrf") {
    X_new <- matrix(x, nrow = 1)
    colnames(X_new) <- cqr_fit$feature_names
    as.numeric(stats::predict(cqr_fit$qrf, X_new, what = cqr_fit$tau))

  } else if (cqr_fit$type == "lm") {
    newd <- as.data.frame(matrix(x, nrow = 1))
    colnames(newd) <- cqr_fit$feature_names
    as.numeric(stats::predict(cqr_fit$lmfit, newd))

  } else {
    stop("Unknown cqr_fit type.")
  }
}

#' Conformal upper bound for pi0 using CQR and density-ratio stabilisation
#'
#' Given a fitted conditional quantile regression model for \code{pi0},
#' this function computes a conformal upper bound at a target feature
#' vector \code{x0}. It uses residuals on calibration data, reweighted by
#' an estimated density ratio between calibration and target feature
#' distributions, with optional clipping.
#'
#' @importFrom stats as.formula approx
#' @param cqr_fit Fitted object from [train_cqr_pi0()].
#' @param cal_df Calibration data frame containing \code{pi0} and features.
#' @param feature_names Character vector of feature names used in the CQR.
#' @param x0 Numeric feature vector at which to evaluate the upper bound.
#'   Must be named and at least contain \code{feature_names}.
#' @param X_target_all Matrix or data frame of target feature rows,
#'   used for density-ratio (logistic) fitting.
#' @param alpha Miscoverage level for the conformal correction.
#' @param weight_clip Numeric length-2 vector giving lower and upper
#'   clipping bounds for the importance weights.
#'
#' @return A scalar upper bound on \code{pi0} in \eqn{[0,1]}.
#' @export
cqr_upper_bound <- function(cqr_fit, cal_df, feature_names, x0,
                            X_target_all, alpha = 0.05,
                            weight_clip = c(0.2, 2)) {
  # Extract calibration and target feature matrices
  X_cal_full <- as.matrix(cal_df[, feature_names, drop = FALSE])
  X_tar_full <- as.matrix(X_target_all[, feature_names, drop = FALSE])

  n_c <- nrow(X_cal_full)
  n_t <- nrow(X_tar_full)
  if (n_c == 0 || n_t == 0) {
    return(1)
  }

  # Match sample sizes for logistic regression
  if (n_c > n_t) {
    sel_c <- sample.int(n_c, n_t)
    X_cal <- X_cal_full[sel_c, , drop = FALSE]
    X_tar <- X_tar_full
  } else if (n_t > n_c) {
    sel_c <- seq_len(n_c)
    X_cal <- X_cal_full
    X_tar <- X_tar_full[sample.int(n_t, n_c), , drop = FALSE]
  } else {
    sel_c <- seq_len(n_c)
    X_cal <- X_cal_full
    X_tar <- X_tar_full
  }

  # Standardise features for logistic regression only
  X_all <- rbind(X_cal, X_tar)
  mu <- colMeans(X_all, na.rm = TRUE)
  sdv <- apply(X_all, 2, stats::sd)
  sdv[sdv < 1e-8] <- 1

  X_cal_s <- scale(X_cal, center = mu, scale = sdv)
  X_tar_s <- scale(X_tar, center = mu, scale = sdv)

  z_cal <- rep(0L, nrow(X_cal_s))
  z_tar <- rep(1L, nrow(X_tar_s))

  df_lr <- data.frame(
    z = c(z_cal, z_tar),
    rbind(X_cal_s, X_tar_s)
  )
  colnames(df_lr) <- c("z", feature_names)

  fit_ok <- TRUE
  p_cal <- NULL

  # Fit logistic regression for density ratio
  suppressWarnings({
    gl <- try(
      stats::glm(z ~ ., data = df_lr, family = stats::binomial(), control = list(maxit = 200)),
      silent = TRUE
    )
    if (inherits(gl, "try-error")) {
      fit_ok <- FALSE
    } else {
      newd <- as.data.frame(X_cal_s)
      colnames(newd) <- feature_names
      p_cal <- try(stats::predict(gl, newdata = newd, type = "response"), silent = TRUE)
      if (inherits(p_cal, "try-error") || any(!is.finite(p_cal))) {
        fit_ok <- FALSE
      }
    }
  })

  # Check for extreme probabilities
  if (fit_ok &&
      (min(p_cal, na.rm = TRUE) < 1e-6 ||
       max(p_cal, na.rm = TRUE) > 1 - 1e-6)) {
    fit_ok <- FALSE
  }

  # Compute clipped importance weights
  if (fit_ok) {
    w_raw <- p_cal / pmax(1 - p_cal, 1e-12)
    w <- pmin(pmax(w_raw, weight_clip[1]), weight_clip[2])
  } else {
    w <- rep(1, nrow(X_cal_s))
  }

  # Residuals on RAW X_cal for upper bound
  qhat_cal <- apply(X_cal, 1, function(xx) {
    predict_qtau(cqr_fit, stats::setNames(as.numeric(xx), feature_names))
  })
  y_bal <- cal_df$pi0[sel_c]
  r <- y_bal - qhat_cal

  # Weighted (1 - alpha) quantile of residuals
  ord <- order(r)
  r_ord <- r[ord]
  w_ord <- w[ord]
  cw <- cumsum(w_ord) / sum(w_ord)

  c_hat <- approx(
    x = c(0, cw),
    y = c(r_ord[1], r_ord),
    xout = 1 - alpha,
    ties = "ordered",
    rule = 2
  )$y
  c_hat <- max(c_hat, 0)

  # Upper bound at x0
  qhat_x0 <- predict_qtau(cqr_fit, stats::setNames(as.numeric(x0), feature_names))
  U <- qhat_x0 + c_hat
  clip01(U, 1e-8)
}
