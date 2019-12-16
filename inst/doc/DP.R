## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)


## ----setup--------------------------------------------------------------------
library(PowerTOST) # attach the library


## -----------------------------------------------------------------------------
mod.fibo <- function(lowest, levels) {
  # modified Fibonacci dose-escalation
  fib      <- c(2, 1 + 2/3, 1.5, 1 + 1/3)
  doses    <- numeric(levels)
  doses[1] <- lowest
  level    <- 2
  repeat {
    if (level <= 4) {
      doses[level] <- doses[level-1] * fib[level-1]
    } else {  # ratio 1.33 for all higher doses
      doses[level] <- doses[level-1] * fib[4]
    }
    level <- level + 1
    if (level > levels) {
      break
    }
  }
  return(signif(doses, 3))
}
lowest <- 10
levels <- 3
doses  <- mod.fibo(lowest, levels)
sampleN.dp(CV = 0.20, doses = doses, beta0 = 1.02)


## -----------------------------------------------------------------------------
levels <- 4
doses  <- mod.fibo(lowest, levels)
x <- sampleN.dp(CV = 0.20, doses = doses, beta0 = 1.02) # we need the data.frame later


## -----------------------------------------------------------------------------
res <- data.frame(n = seq(x[["Sample size"]], 12, -1),
                  power = NA)
for (i in 1:nrow(res)) {
  res$power[i] <- signif(suppressMessages(
                           power.dp(CV = 0.20, doses = doses,
                                    beta0 = 1.02,
                                    n = res$n[i])), 5)
}
res <- res[res$power >= 0.80, ]
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
levels <- 5
doses  <- mod.fibo(lowest, levels)
per    <- 3
block  <- levels*(levels-1)/(per-1)
dm     <- crossdes::balmin.RMD(levels, block, per)
x      <- sampleN.dp(CV = 0.20, doses = doses, beta0 = 1.02,
                     design = "IBD", dm = dm)


## -----------------------------------------------------------------------------
res <- data.frame(n = seq(x[["Sample size"]], nrow(dm), -1),
                  power = NA)
for (i in 1:nrow(res)) {
  res$power[i] <- signif(suppressMessages(
                           power.dp(CV = 0.20, doses = doses,
                                    beta0 = 1.02, design = "IBD",
                                    dm = dm, n = res$n[i])), 5)
}
res <- res[res$power >= 0.80, ]
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
doses  <- 2^(seq(0, 8, 2))
levels <- length(doses)
sampleN.dp(CV = 0.30, doses = doses, beta0 = 1.02,
           design = "crossover")


## -----------------------------------------------------------------------------
sampleN.dp(CV = 0.30, doses = doses, beta0 = 1.02,
           design = "crossover", theta1 = 0.75)
