## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)


## ----setup--------------------------------------------------------------------
library(PowerTOST) # attach the library


## -----------------------------------------------------------------------------
sampleN.TOST(CV = 0.30)


## -----------------------------------------------------------------------------
sampleN.TOST(CV = 0.30, details = FALSE, print = FALSE)[["Sample size"]]


## -----------------------------------------------------------------------------
power.TOST(CV = 0.30, n = 39)


## -----------------------------------------------------------------------------
designs <- c("2x2x2", "2x2x3", "2x3x3", "2x2x4")
# data.frame of results
res <- data.frame(design = designs, n = NA, power = NA, n.do = NA,
                  power.do = NA, stringsAsFactors = FALSE)
for (i in 1:4) {
  # print = FALSE and details = FALSE suppress output to the console
  # we are only interested in columns 7-8
  # let's also calculate power for one dropout
  res[i, 2:3] <- signif(
                   sampleN.TOST(CV = 0.30, design = res$design[i],
                                print = FALSE, details = FALSE)[7:8], 6)
  res[i, 4]   <-  res[i, 2] - 1
  res[i, 5]   <- suppressMessages(
                   signif(
                     power.TOST(CV = 0.30, design = res$design[i],
                                n = res[i, 4]), 6))
}
print(res, row.names = FALSE)


## ---- echo = FALSE------------------------------------------------------------
designs <- c("2x2x2", "2x2x3", "2x3x3", "2x2x4")
res <- data.frame(design = rep(NA, 4), name = NA, n = NA, formula = NA,
                  df = NA, t.value = NA, stringsAsFactors = FALSE)
res[, c(2:1, 4)] <- known.designs()[which(known.designs()[, 2] %in% designs),
                                    c(9, 2:3)]
for (i in 1:4) {
  res$n[i]  <- sampleN.TOST(CV = 0.30, design = res$design[i],
                            print = FALSE, details = FALSE)[["Sample size"]]
  e         <- parse(text=res[i, 4], srcfile=NULL)
  n         <- res$n[i]
  res[i, 5] <- eval(e)
  res$t.value[i] <- signif(qt(1-0.05, df = res[i, 5]), 4)
}
res <- res[with(res, order(-n, design)), ]
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
ngrp <- function(capacity, n) {
  # split sample size into >=2 groups based on capacity
  if (n <= capacity) { # make equal groups
    ngrp <- rep(ceiling(n/2), 2)
  } else {             # at least one = capacity
    ngrp    <- rep(0, ceiling(n / capacity))
    grps    <- length(ngrp)
    ngrp[1] <- capacity
    for (j in 2:grps) {
      n.tot <- sum(ngrp) # what we have so far
      if (n.tot + capacity <= n) {
        ngrp[j] <- capacity
      } else {
        ngrp[j] <- n - n.tot
      }
    }
  }
  return(ngrp = list(grps = length(ngrp), ngrp = ngrp))
}
CV        <- 0.30
capacity  <- 24 # clinical capacity
res       <- data.frame(n = NA, grps = NA, pwr.1 = NA, pwr.2 = NA)
x         <- sampleN.TOST(CV = CV, print = FALSE, details = FALSE)
res$n     <- x[["Sample size"]]
res$pwr.1 <- x[["Achieved power"]]
x         <- ngrp(capacity = capacity, n = res$n)
res$grps  <- x[["grps"]]
ngrp      <- x[["ngrp"]]
res$pwr.2 <- power.TOST.sds(CV = CV, n = res$n, grps = res$grps,
                            ngrp = ngrp, gmodel = 2, progress = FALSE)
res$loss  <- signif(100*(res$pwr.2 - res$pwr.1)/res$pwr.1, 5)
print(signif(res, 6), row.names = FALSE)


## -----------------------------------------------------------------------------
sampleN.RatioF(CV = 0.20, CVb = 0.40)


## -----------------------------------------------------------------------------
sampleN.TOST(CV = 0.20, theta0 = 0.92)


## -----------------------------------------------------------------------------
df <- 16 - 2 # degrees of freedom of the 2x2x2 crossover pilot
CVCL(CV = 0.20, df = df, side = "upper", alpha = 0.20)[["upper CL"]]


## -----------------------------------------------------------------------------
CL.upper <- CVCL(CV = 0.20, df = 16 - 2, side = "upper",
                 alpha = 0.20)[["upper CL"]]
res <- sampleN.TOST(CV = CL.upper, theta0 = 0.92, print = FALSE)
print(res[7:8], row.names = FALSE)


## -----------------------------------------------------------------------------
CL.upper <- CVCL(CV = 0.20, df = 16 - 2, side = "upper",
                 alpha = 0.20)[["upper CL"]]
power.TOST(CV = CL.upper, theta0 = 0.92, n = 28)


## -----------------------------------------------------------------------------
power.TOST(CV = 0.22, theta0 = 0.90, n = 40)


## -----------------------------------------------------------------------------
expsampleN.TOST(CV = 0.20, theta0 = 0.92, prior.type = "CV",
                prior.parm = list(m = 16, design = "2x2x2"))


## -----------------------------------------------------------------------------
expsampleN.TOST(CV = 0.20, theta0 = 0.92, prior.type = "theta0",
                prior.parm = list(m = 16, design = "2x2x2"))


## -----------------------------------------------------------------------------
expsampleN.TOST(CV = 0.20, theta0 = 0.92, prior.type = "both",
                prior.parm = list(m = 16, design = "2x2x2"),
                details = FALSE)


## -----------------------------------------------------------------------------
CV  <- 0.20
res <- data.frame(design = c("3x6x3", "2x2x2"), n = NA, power = NA,
                  stringsAsFactors = FALSE)
for (i in 1:2) {
  res[i, 2:3] <- signif(
                   sampleN.TOST(CV = CV, design = res$design[i],
                                details = FALSE, print = FALSE)[7:8], 5)
}
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
CV  <- 0.20
res <- data.frame(design = c("4x4", "2x2x2"), n = NA, power = NA,
                  stringsAsFactors = FALSE)
for (i in 1:2) {
  res[i, 2:3] <- signif(
                  sampleN.TOST(CV = CV, design = res$design[i],
                               details = FALSE, print = FALSE)[7:8], 5)
}
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
n     <- 28
n.seq <- rep(n/2, 2) - c(1, 2)
round(100*CI.BE(pe = 0.90, CV = 0.25, n = n.seq), 2)


## -----------------------------------------------------------------------------
power.TOST(CV = 0.25, theta0 = 0.90, n = c(13, 12))


## -----------------------------------------------------------------------------
power.TOST(CV = 0.25, theta0 = 1, n = c(13, 12))


## -----------------------------------------------------------------------------
CVs <- ("  CV |  n | design | study
         0.20 | 16 |  2x2x2 | pilot
         0.25 | 25 |  2x2x2 | pivotal")
txtcon <- textConnection(CVs)
data   <- read.table(txtcon, header = TRUE, sep = "|",
                     strip.white = TRUE, as.is = TRUE)
close(txtcon)
print(CVpooled(data, alpha = 0.20), digits = 4, verbose = TRUE)


## -----------------------------------------------------------------------------
pa.ABE(CV = 0.20, theta0 = 0.92)


## -----------------------------------------------------------------------------
balance <- function(x, seqs) {
  x <- ceiling(x) + ceiling(x) %% seqs
  return(x)
}
do.rate   <- 0.10
seqs      <- 3
n         <- seq(12, 120, 12)
res       <- data.frame(n = n,
                        ad1 = balance(n * (1 + do.rate), seqs),
                        el1 = NA, diff1 = NA,
                        ad2 = balance(n / (1 - do.rate), seqs),
                        el2 = NA, diff2 = NA)
res$el1   <- floor(res$ad1 * (1 - do.rate))
res$diff1 <- sprintf("%+i", res$el1 - n)
res$el2   <- floor(res$ad2 * (1 - do.rate))
res$diff2 <- sprintf("%+i", res$el2 - n)
invisible(
  ifelse(res$el2 - n >=0, res$opt <- res$el2, res$opt <- res$el1))
res$diff  <- sprintf("%+i", res$opt - n)
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
CVfromCI(lower = 0.8323, upper = 1.0392,
         design = "2x2x4", n = 26)


## -----------------------------------------------------------------------------
n     <- 26
n1    <- balance(seq(n, 12, -1), 2) / 2
n2    <- n - n1
nseqs <- unique(data.frame(n1 = n1, n2 = n2, n = n))
res   <- data.frame(n1 = nseqs$n1, n2 = nseqs$n2, CV = NA)
for (i in 1:nrow(res)) {
  res$CV[i] <- CVfromCI(lower = 0.8323, upper = 1.0392,
                        design = "2x2x4",
                        n = c(res$n1[i], res$n2[i]))
}
print(res, row.names = FALSE)


## -----------------------------------------------------------------------------
CV   <- 0.21
d    <- 0.05 # delta 5%, direction unknown
n    <- sampleN.TOST(CV = CV, theta0 = 1 - d, print = FALSE,
                     details = FALSE)[["Sample size"]]
res1 <- data.frame(CV = CV, theta0 = c(1 - d, 1 / (1 - d)),
                   n = n, power = NA)
for (i in 1:nrow(res1)) {
  res1$power[i] <- power.TOST(CV = CV, theta0 = res1$theta0[i], n = n)
}
n    <- sampleN.TOST(CV = CV, theta0 = 1 + d, print = FALSE,
                     details = FALSE)[["Sample size"]]
res2 <- data.frame(CV = CV, theta0 = c(1 + d, 1 / (1 + d)),
                   n = n, power = NA)
for (i in 1:nrow(res1)) {
  res2$power[i] <- power.TOST(CV = CV, theta0 = res2$theta0[i], n = n)
}
res <- rbind(res1, res2)
print(signif(res[order(res$n, res$theta0), ], 4), row.names = FALSE)
