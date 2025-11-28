prepare_conformal_cqr <- function(
    cqr_fit,
    cal_df,
    feature_names,
    X_target_all,
    alpha = 0.05,
    weight_clip = c(0.2, 2)
) {

  # ---------------------------
  # 0. Extract matrices
  # ---------------------------
  X_cal_full <- as.matrix(cal_df[, feature_names, drop = FALSE])
  X_tar_full <- as.matrix(X_target_all[, feature_names, drop = FALSE])
  y_cal_full <- cal_df$pi0

  n_c <- nrow(X_cal_full)
  n_t <- nrow(X_tar_full)
  if (n_c == 0 || n_t == 0)
    return(list(c_hat = 1, mu = NULL, sdv = NULL))

  # ---------------------------
  # 1. Balance sample sizes
  # ---------------------------
  if (n_c > n_t) {
    idx_cal <- sample.int(n_c, n_t)
    X_cal <- X_cal_full[idx_cal, , drop = FALSE]
    y_cal <- y_cal_full[idx_cal]
    X_tar <- X_tar_full
  } else {
    idx_cal <- seq_len(n_c)
    X_cal <- X_cal_full
    y_cal <- y_cal_full
    X_tar <- X_tar_full[sample.int(n_t, n_c), , drop = FALSE]
  }

  # ---------------------------
  # 2. Standardize (once)
  # ---------------------------
  X_all <- rbind(X_cal, X_tar)
  mu <- colMeans(X_all)
  sdv <- apply(X_all, 2, sd)
  sdv[sdv < 1e-8] <- 1

  X_cal_s <- scale(X_cal, center = mu, scale = sdv)
  X_tar_s <- scale(X_tar, center = mu, scale = sdv)

  # ---------------------------
  # 3. Fit logistic regression (once)
  # ---------------------------
  df_lr <- data.frame(
    z = c(rep(0L, nrow(X_cal_s)), rep(1L, nrow(X_tar_s))),
    rbind(X_cal_s, X_tar_s)
  )
  colnames(df_lr) <- c("z", feature_names)

  fit_ok <- TRUE
  suppressWarnings({
    gl <- try(glm(z ~ ., data = df_lr, family = binomial(),
                  control = list(maxit = 200)),
              silent = TRUE)
    if (inherits(gl, "try-error"))
      fit_ok <- FALSE
  })

  # ---------------------------
  # 4. Density-ratio weights (once)
  # ---------------------------
  if (fit_ok) {
    p_cal <- try(predict(gl, newdata = as.data.frame(X_cal_s),
                         type = "response"), silent = TRUE)

    if (inherits(p_cal, "try-error") ||
        any(!is.finite(p_cal)) ||
        min(p_cal) < 1e-6 ||
        max(p_cal) > 1 - 1e-6) {
      fit_ok <- FALSE
    }
  }

  if (!fit_ok) {
    w <- rep(1, nrow(X_cal_s))
  } else {
    w_raw <- p_cal / pmax(1 - p_cal, 1e-12)
    w <- pmin(pmax(w_raw, weight_clip[1]), weight_clip[2])
  }

  # ---------------------------
  # 5. Calibration residuals (once)
  # ---------------------------
  qhat_cal <- apply(X_cal, 1, function(xx) {
    scoutFDR:::predict_qtau(cqr_fit, setNames(as.numeric(xx), feature_names))
  })

  r <- y_cal - qhat_cal
  ok <- is.finite(r) & is.finite(w)
  r <- r[ok]
  w <- w[ok]

  if (length(r) < 2)
    return(list(c_hat = 1, mu = mu, sdv = sdv))

  # ---------------------------
  # 6. Weighted calibration quantile (once)
  # ---------------------------
  ord <- order(r)
  r_ord <- r[ord]
  w_ord <- w[ord]

  cw <- cumsum(w_ord) / sum(w_ord)
  x_vec <- c(0, cw)
  y_vec <- c(r_ord[1], r_ord)

  ok2 <- is.finite(x_vec) & is.finite(y_vec)
  x_vec <- x_vec[ok2]
  y_vec <- y_vec[ok2]

  if (length(x_vec) < 2)
    return(list(c_hat = 1, mu = mu, sdv = sdv))

  c_hat <- stats::approx(
    x_vec, y_vec,
    xout = 1 - alpha,
    ties = "ordered",
    rule = 2
  )$y
  c_hat <- max(c_hat, 0)

  # ---------------------------
  # Return preparation object
  # ---------------------------
  list(
    c_hat = c_hat,
    mu = mu,
    sdv = sdv,
    feature_names = feature_names,
    cqr_fit = cqr_fit
  )
}


predict_conformal_cqr <- function(prep, x0) {

  q_x0 <- scoutFDR:::predict_qtau(
    prep$cqr_fit,
    setNames(as.numeric(x0), prep$feature_names)
  )

  out <- q_x0 + prep$c_hat
  pmin(pmax(out, 0), 1)
}
