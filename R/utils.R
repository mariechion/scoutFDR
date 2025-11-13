suppressPackageStartupMessages({
  have_cp4p <- requireNamespace("cp4p", quietly = TRUE)
  have_qr   <- requireNamespace("quantreg", quietly = TRUE)
  have_qrf  <- requireNamespace("quantregForest", quietly = TRUE)
  have_MASS <- requireNamespace("MASS", quietly = TRUE)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

clip01 <- function(x, eps = 1e-8) {
  pmin(pmax(as.numeric(x), eps), 1)
}