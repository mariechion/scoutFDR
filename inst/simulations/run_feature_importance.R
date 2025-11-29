# ============================================================
# SCOUT-FDR — Feature Importance (Global + Scenario-Weighted)
# ============================================================

library(scoutFDR)
library(quantreg)
library(ggplot2)
library(dplyr)
library(purrr)

# ----- Load calibration + model -----
cal_df <- readRDS("inst/calibration/cal_df.rds")
cqr_fit <- readRDS("inst/calibration/cqr_fit.rds")
feature_names <- cqr_fit$feature_names

# ----- Pinball loss -----
pinball_loss <- function(y, yhat, tau, weights=NULL) {
  u <- y - yhat
  loss <- ifelse(u >= 0, tau*u, (1-tau)*(-u))
  if (is.null(weights)) mean(loss) else mean(loss * weights)
}

# ----- Permutation importance -----
compute_feature_importance_perm <- function(cqr_fit, cal_df, weights=NULL,
                                            repeats=100, seed=123) {

  feat_names <- cqr_fit$feature_names
  keep_rows <- complete.cases(cal_df[, c("pi0", feat_names)])
  df <- cal_df[keep_rows, ]
  y <- df$pi0
  X <- as.matrix(df[, feat_names])

  mu <- cqr_fit$mu %||% rep(0, length(feat_names))
  sd <- cqr_fit$sd %||% rep(1, length(feat_names))
  sd[sd < 1e-8] <- 1

  Xs <- scale(X, center=mu, scale=sd)
  if (!is.null(cqr_fit$vars)) colnames(Xs) <- cqr_fit$vars

  tau <- cqr_fit$tau %||% 0.5

  # Identify model
  model_type <- case_when(
    !is.null(cqr_fit$rqfit) ~ "rq",
    !is.null(cqr_fit$qrf)   ~ "qrf",
    !is.null(cqr_fit$lmfit) ~ "lm",
    TRUE                    ~ "custom"
  )

  baseline_pred <- switch(model_type,
    rq  = { b <- coef(cqr_fit$rqfit); as.numeric(b[1] + Xs %*% b[-1]) },
    qrf = predict(cqr_fit$qrf, X, what=tau),
    lm  = predict(cqr_fit$lmfit, as.data.frame(X)),
    custom = apply(X, 1, function(xx) predict_qtau(cqr_fit, setNames(xx, feat_names)))
  )

  L0 <- pinball_loss(y, baseline_pred, tau, weights)

  set.seed(seed)
  importance <- map_dbl(seq_along(feat_names), function(j) {
    deltas <- map_dbl(seq_len(repeats), function(r) {
      Xs_perm <- Xs; Xs_perm[, j] <- sample(Xs_perm[, j])

      yhat_p <- switch(model_type,
        rq  = { b <- coef(cqr_fit$rqfit); as.numeric(b[1] + Xs_perm %*% b[-1]) },
        qrf = {
          Xp <- X; Xp[, j] <- sample(Xp[, j])
          predict(cqr_fit$qrf, Xp, what=tau)
        },
        lm  = {
          Xp <- X; Xp[, j] <- sample(Xp[, j])
          predict(cqr_fit$lmfit, as.data.frame(Xp))
        },
        custom = {
          Xp <- X; Xp[, j] <- sample(Xp[, j])
          apply(Xp, 1, function(xx) predict_qtau(cqr_fit, setNames(xx, feat_names)))
        }
      )
      pinball_loss(y, yhat_p, tau, weights) - L0
    })
    mean(deltas)
  })

  data.frame(
    feature = feat_names,
    importance = importance,
    type = if (is.null(weights)) "global" else "scenario_weighted",
    stringsAsFactors = FALSE
  )
}

# ----- Density-ratio weights -----
estimate_dr_weights <- function(cal_df, feat_for_weight, X_target_all,
                                clip=c(1/3, 3)) {

  feat_for_weight <- intersect(feat_for_weight,
                               intersect(colnames(cal_df), colnames(X_target_all)))
  if (length(feat_for_weight)==0) return(rep(1, nrow(cal_df)))

  X_cal <- as.matrix(cal_df[, feat_for_weight, drop=FALSE])
  X_tar <- as.matrix(X_target_all[, feat_for_weight, drop=FALSE])

  n_c <- nrow(X_cal); n_t <- nrow(X_tar)
  if (n_c==0 || n_t==0) return(rep(1, nrow(cal_df)))

  if (n_c > n_t) X_cal <- X_cal[sample(n_c, n_t), ]
  if (n_t > n_c) X_tar <- X_tar[sample(n_t, n_c), ]

  X_all <- rbind(X_cal, X_tar)
  mu <- colMeans(X_all)
  sd <- apply(X_all, 2, sd); sd[sd < 1e-8] <- 1

  Xc <- scale(X_cal, center=mu, scale=sd)
  Xt <- scale(X_tar, center=mu, scale=sd)

  df_lr <- data.frame(z=c(rep(0, nrow(Xc)), rep(1, nrow(Xt))),
                      rbind(Xc, Xt))
  colnames(df_lr) <- c("z", feat_for_weight)

  gl <- try(glm(z ~ ., data=df_lr, family=binomial()), silent=TRUE)
  if (inherits(gl, "try-error")) return(rep(1, nrow(cal_df)))

  Xc_full <- scale(as.matrix(cal_df[, feat_for_weight]), center=mu, scale=sd)
  p_cal <- try(predict(gl, newdata=as.data.frame(Xc_full), type="response"), silent=TRUE)
  if (inherits(p_cal, "try-error")) return(rep(1, nrow(cal_df)))

  w_raw <- p_cal / pmax(1 - p_cal, 1e-12)
  pmin(pmax(w_raw, clip[1]), clip[2])
}

# ----- Target feature extraction -----
build_target_features <- function(p, K, feature_names) {
  fold_id <- sample(rep(seq_len(K), length.out=length(p)))
  map_dfr(seq_len(K), function(k) {
    fk <- try(features_ecdf(p[fold_id == k]), silent=TRUE)
    row <- rep(NA_real_, length(feature_names)); names(row) <- feature_names
    if (!inherits(fk,"try-error")) {
      nm <- intersect(names(fk), feature_names)
      row[nm] <- as.numeric(fk[nm])
    }
    row
  })
}

# ============================================================
# 1. Global feature importance
# ============================================================

imp_global <- compute_feature_importance_perm(
  cqr_fit=cqr_fit,
  cal_df=cal_df,
  weights=NULL,
  repeats=100
)

write.csv(imp_global,
          "inst/simulations/feature_importance_cont_perm_global.csv",
          row.names=FALSE)

# Plot
ig <- imp_global[order(-imp_global$importance), ]
topK <- min(30, nrow(ig))
ig_top <- ig[1:topK, ]
ig_top$feature <- factor(ig_top$feature, levels=rev(ig_top$feature))

p <- ggplot(ig_top, aes(x=feature, y=importance)) +
  geom_col(fill="#2c7fb8") +
  coord_flip() +
  theme_minimal(base_size=12)

ggsave("inst/simulations/feature_importance_cont_perm_global.png",
       p, width=8, height=10, dpi=300)

# ============================================================
# 2. Scenario-weighted importance
# ============================================================

scenarios_extra <- data.frame(
  scen_id = paste0("L",1:15),
  n    = c(3,3,3, 5,5,5, 5,5,5, 8,8,8, 10,10,10),
  pi0  = c(0.2,0.4,0.7, 0.2,0.4,0.7, 0.2,0.4,0.7, 0.2,0.4,0.7, 0.2,0.4,0.7),
  power= c(0.20,0.20,0.20, 0.20,0.20,0.20, 0.35,0.35,0.35, 0.35,0.35,0.35, 0.50,0.50,0.50),
  stringsAsFactors=FALSE
)

for (i in seq_len(nrow(scenarios_extra))) {
  scen <- scenarios_extra[i, ]

  p_tar <- scoutFDR:::simulate_pvalues_ttest(
    m=5000, n=scen$n, pi0=scen$pi0, target_power=scen$power
  )

  X_tar <- build_target_features(p_tar, K=2, feature_names=feature_names)

  keep_cols <- vapply(X_tar, function(v) all(is.finite(v)), logical(1))
  feat_ok <- intersect(feature_names, names(X_tar)[keep_cols])

  if (length(feat_ok) < 5) next

  w <- estimate_dr_weights(
    cal_df=cal_df,
    feat_for_weight=feat_ok,
    X_target_all=X_tar[, feat_ok, drop=FALSE]
  )

  imp_w <- compute_feature_importance_perm(
    cqr_fit=cqr_fit, cal_df=cal_df, weights=w, repeats=100
  )

  outfile <- sprintf(
    "inst/simulations/feature_importance_cont_perm_weighted_%s.csv",
    scen$scen_id
  )
  write.csv(imp_w, outfile, row.names=FALSE)
}

cat("\nFeature importance analyses complete.\n")
