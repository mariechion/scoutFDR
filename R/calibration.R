#' Simulate two-sample t-test p-values
#'
#' Simulates p-values from independent two-sample t-tests under a mixture
#' model with a proportion \code{pi0} of null hypotheses and a target
#' power for the alternatives.
#'
#' @param m Integer, number of hypotheses.
#' @param n Integer, sample size per group for each t-test.
#' @param pi0 Numeric in \eqn{[0,1]}, proportion of null hypotheses.
#' @param target_power Numeric target power for non-null effects
#'   (used to calibrate the distribution of effect sizes).
#'
#' @return A numeric vector of length \code{m} of p-values.
#' @export
simulate_pvalues_ttest <- function(m, n, pi0, target_power) {
  m0 <- round(m * pi0)
  m1 <- m - m0

  # Calibrate mean effect size to hit target power
  z_power <- stats::qnorm(pmin(pmax(target_power, 1e-4), 1 - 1e-4))
  z_alpha <- stats::qnorm(0.975)
  delta_mean <- (z_alpha + z_power) / sqrt(n / 2)
  delta_sd <- abs(delta_mean) * 0.3

  p_values <- numeric(m)

  # Nulls
  if (m0 > 0) {
    for (i in 1:m0) {
      x1 <- stats::rnorm(n, 0, 1)
      x2 <- stats::rnorm(n, 0, 1)
      p_values[i] <- stats::t.test(x1, x2)$p.value
    }
  }

  # Alternatives
  if (m1 > 0) {
    for (i in (m0 + 1):m) {
      delta <- stats::rnorm(1, delta_mean, delta_sd)
      x1 <- stats::rnorm(n, 0, 1)
      x2 <- stats::rnorm(n, delta, 1)
      p_values[i] <- stats::t.test(x1, x2)$p.value
    }
  }

  p_values
}

#' Generate calibration data for pi0 regression
#'
#' Repeatedly simulates p-values from a t-test model over a grid of
#' \code{(pi0, m, n, target_power)} settings, computes a base pi0 estimate
#' and user-specified features, and returns a data frame. This can be used
#' to train conditional models for pi0, e.g. via quantile regression.
#'
#' @param n_sim Integer, number of simulated datasets.
#' @param grid Optional data frame specifying columns \code{pi0}, \code{m},
#'   \code{n}, and \code{target_power}. If \code{NULL}, a default grid is used.
#' @param base_pi0 A function factory as returned by [make_base_pi0()],
#'   used to compute a scalar base estimate \code{A} from each p-vector.
#' @param feature_fun Function mapping a numeric p-value vector to a
#'   named list or vector of features. Defaults to [features_ecdf()].
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A data frame with one row per simulation, containing columns
#'   \code{pi0}, \code{A}, and the chosen features. All columns are numeric,
#'   making this object easy to use with tidyverse tools.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
generate_calibration_data <- function(n_sim = 2000,
                                      grid = NULL,
                                      base_pi0 = make_base_pi0("ecdf_storey_median"),
                                      feature_fun = function(p) as.list(features_ecdf(p)),
                                      seed = 1) {
  if (is.null(grid)) {
    grid <- expand.grid(
      pi0 = seq(0.1, 0.7, by = 0.05),
      m = seq(1000, 15000, by = 100),
      n = 3:15,
      target_power = seq(0.05, 0.7, by = 0.05)
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  out <- vector("list", n_sim)
  nG <- nrow(grid)

  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  on.exit(close(pb), add = TRUE)

  for (i in seq_len(n_sim)) {
    # Randomly pick a grid row
    g <- grid[sample.int(nG, 1), ]

    # Simulate p-values
    p <- simulate_pvalues_ttest(
      m = g$m,
      n = g$n,
      pi0 = g$pi0,
      target_power = g$target_power
    )

    # Base estimator and features
    A <- base_pi0(p)

    feats <- feature_fun(p)
    row <- c(pi0 = g$pi0, A = A, unlist(feats))

    storage.mode(row) <- "double"
    out[[i]] <- row

    setTxtProgressBar(pb, i)
  }

  df <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)

  # Ensure all columns are numeric (tidyverse-friendly)
  for (nm in colnames(df)) {
    df[[nm]] <- as.numeric(df[[nm]])
  }

  df
}
