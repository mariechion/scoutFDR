library(scoutFDR)

# --- Build calibration data --- #

start_time_gen <- Sys.time()
cal_df <- generate_calibration_data(n_sim = 20000, seed = 17)
end_time_gen <- Sys.time()

saveRDS(cal_df, file = "inst/calibration/cal_df.rds")
message("Calibration data generated in: ", end_time_gen - start_time_gen)


# --- Fit CQR model --- #

start_time_train <- Sys.time()
cqr_fit <- train_cqr_pi0(
  cal_df        = cal_df,
  feature_names = "auto",
  tau           = 0.7,
  model         = "qr"
)
end_time_train <- Sys.time()


saveRDS(cqr_fit, file = "inst/calibration/cqr_fit.rds")
message("Saved calibration objects to inst/calibration/")
message("CQR fitting time: ", end_time_train - start_time_train)


# --- Optional diagnostics --- #
feat_names <- cqr_fit$feature_names
x_vec <- as.numeric(cal_df[1, feat_names])
names(x_vec) <- feat_names
pred <- scoutFDR:::predict_qtau(cqr_fit, x_vec)
print(pred)

resids <- apply(cal_df, 1, function(row) {
  xv <- as.numeric(row[feat_names])
  names(xv) <- feat_names
  qhat <- scoutFDR:::predict_qtau(cqr_fit, xv)
  row["pi0"] - qhat
})
print(summary(resids))

base_pi0    <- make_base_pi0("ecdf_storey_median")
feature_fun <- function(pp) as.list(features_ecdf(pp))
p     <- runif(2000)
A     <- base_pi0(p)
feats <- feature_fun(p)

x0 <- unlist(c(A = A, feats))
x0 <- x0[feat_names]

U <- cqr_upper_bound(
  cqr_fit       = cqr_fit,
  cal_df        = cal_df,
  feature_names = feat_names,
  x0            = x0,
  X_target_all  = as.data.frame(t(x0)),
  alpha         = 0.05
)
print(U)
