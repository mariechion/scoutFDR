# --- Build calibration data --- #

start_time_gen <- Sys.time()
cal_df <- generate_calibration_data(n_sim = 10000, seed = 17)
end_time_gen <- Sys.time()
saveRDS(cal_df, file = "inst/calibration/cal_df.rds")

message("Calibration data generated in: ", end_time_gen - start_time_gen)

# --- Fit CQR model --- #

feat_names <- c(
  "A",
  "x9","x10","x9_05","x9_10","x9_20","x9_curv",
  "A_lam60","A_lam70","A_lam80","A_lam90","A_slope","A_range",
  "F050","ER050",
  "R0","t0","Simes_min",
  "KS","HC",
  "lnp_med","lnp_mean","lnp_p90",
  "m","uniq_ratio","boundary_mass"
)

feat_names <- intersect(feat_names, colnames(cal_df))

start_time_train <- Sys.time()
cqr_fit <- train_cqr_pi0(
  cal_df        = cal_df,
  feature_names = feat_names,
  tau           = 0.8,
  model         = "qr"
)
end_time_train <- Sys.time()

saveRDS(cqr_fit, file = "inst/calibration/cqr_fit.rds")
message("Saved calibration objects to inst/calibration/")
message("CQR fitting time: ", end_time_train - start_time_train)

# --- Optional diagnostics --- #
x_vec <- as.numeric(cal_df[1, feat_names])
names(x_vec) <- feat_names
pred <- predict_qtau(cqr_fit, x_vec)
print(pred)

resids <- sapply(seq_len(nrow(cal_df)), function(i) {
  xv <- as.numeric(cal_df[i, feat_names])
  names(xv) <- feat_names
    qhat <- predict_qtau(cqr_fit, xv)
    cal_df$pi0[i] - qhat
})
print(summary(resids))

base_pi0    <- make_base_pi0("ecdf_storey_median")
feature_fun <- function(pp) as.list(features_ecdf(pp))
p     <- runif(2000)
A     <- base_pi0(p)
feats <- feature_fun(p)
x0    <- c(A = A, unlist(feats))

U <- cqr_upper_bound(
  cqr_fit      = cqr_fit,
  cal_df       = cal_df,
  feature_names = feat_names,
  x0           = x0,
  X_target_all = as.data.frame(t(x0)),
  alpha        = 0.05
)
print(U)

