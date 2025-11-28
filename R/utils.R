
# Internal function
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Optional package availability 

suppressPackageStartupMessages({
  have_cp4p  <- requireNamespace("cp4p", quietly = TRUE)
  have_qr    <- requireNamespace("quantreg", quietly = TRUE)
  have_qrf   <- requireNamespace("quantregForest", quietly = TRUE)
  have_MASS  <- requireNamespace("MASS", quietly = TRUE)
})

# Needed because some functions refer to these names non-locally
utils::globalVariables(c(
  "have_cp4p", "have_qr", "have_qrf", "have_MASS"
))
